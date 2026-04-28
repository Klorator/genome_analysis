#!/bin/bash

set -euo pipefail

die() {
  echo "Error: $*" >&2
  exit 1
}

require_commands() {
  local cmd
  for cmd in "$@"; do
    command -v "$cmd" >/dev/null 2>&1 || die "Required command not found: $cmd"
  done
}

check_blast_flag_support() {
  local blast_cmd="$1"
  shift

  local help_text
  help_text="$($blast_cmd -help 2>&1)" || die "Failed to inspect $blast_cmd help output"

  local flag
  for flag in "$@"; do
    if ! grep -q -- "$flag" <<<"$help_text"; then
      die "$blast_cmd does not appear to support $flag"
    fi
  done
}

validate_positive_int() {
  local name="$1"
  local value="$2"

  [[ "$value" =~ ^[1-9][0-9]*$ ]] || die "$name must be a positive integer; got '$value'"
}

validate_positive_number() {
  local name="$1"
  local value="$2"

  [[ "$value" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$|^[0-9]+[eE][+-]?[0-9]+$ ]] || die "$name must be a positive number; got '$value'"
}

preflight_checks() {
  local blast_db="$1"
  local threads="$2"
  local top_n="$3"
  local max_target_seqs="$4"
  local min_identity="$5"
  local min_qcov="$6"
  local evalue="$7"

  require_commands blastn blastdbcmd awk find sort tar cksum dirname basename

  check_blast_flag_support blastn \
    -task \
    -query \
    -db \
    -num_threads \
    -max_target_seqs \
    -max_hsps \
    -perc_identity \
    -qcov_hsp_perc \
    -evalue \
    -outfmt \
    -out \
    -taxids \
    -negative_taxids

  blastdbcmd -db "$blast_db" -info >/dev/null 2>&1 || die "BLAST database '$blast_db' is unavailable"

  validate_positive_int THREADS "$threads"
  validate_positive_int TOP_CONTIGS "$top_n"
  validate_positive_int MAX_TARGET_SEQS "$max_target_seqs"
  validate_positive_number BLAST_PERC_IDENTITY "$min_identity"
  validate_positive_number BLAST_QCOV_HSP_PERC "$min_qcov"
  validate_positive_number BLAST_EVALUE "$evalue"

sync_human_readable_outputs() {
  local source_dir="$1"
  local target_subdir="$2"
  local target_dir="$RESULT_DIR/$target_subdir"

  [ -d "$source_dir" ] || return 0

  mkdir -p "$target_dir"
  while IFS= read -r -d '' file; do
    local rel_path dest_path
    rel_path="${file#${source_dir}/}"
    dest_path="$target_dir/$rel_path"
    mkdir -p "$(dirname "$dest_path")"
    cp -f "$file" "$dest_path"
  done < <(
    find "$source_dir" -type f -size -20M \
      \( -iname "*.html" -o -iname "*.pdf" -o -iname "*.txt" -o -iname "*.tsv" -o -iname "*.csv" -o -iname "*.log" -o -iname "*.json" -o -iname "*.yaml" -o -iname "*.yml" -o -iname "*.md" -o -iname "*.png" -o -iname "*.svg" \) \
      -print0
  )
}

build_top_n_query_fasta() {
  local output_fasta="$1"
  local top_n="$2"
  shift 2

  local -a fasta_files=("$@")
  [ "${#fasta_files[@]}" -gt 0 ] || die "No FASTA files were provided for query reduction"

  awk -v top_n="$top_n" '
    function flush_record(   original, first_token, rest) {
      if (header == "") return

      original = substr(header, 2)
      first_token = original
      rest = ""
      if (match(original, /[ \t]/)) {
        first_token = substr(original, 1, RSTART - 1)
        rest = substr(original, RSTART)
      }

      record_count++
      headers[record_count] = sprintf(">%s|r%d|%s%s", source_tag, file_record_index, first_token, rest)
      lengths[record_count] = seq_len
      seqs[record_count] = seq
    }

    FNR == 1 {
      if (NR > 1) {
        flush_record()
      }
      header = ""
      seq = ""
      seq_len = 0
      file_record_index = 0
      source_tag = FILENAME
      sub(/^.*[\\/]/, "", source_tag)
      gsub(/[^[:alnum:]_.-]/, "_", source_tag)
    }

    /^>/ {
      flush_record()
      header = $0
      seq = ""
      seq_len = 0
      file_record_index++
      next
    }

    {
      gsub(/[ \t\r\n]/, "", $0)
      seq = seq $0
      seq_len += length($0)
    }

    END {
      flush_record()

      if (record_count < 1) {
        exit 1
      }

      for (i = 1; i <= record_count; i++) idx[i] = i
      for (i = 1; i <= record_count; i++) {
        for (j = i + 1; j <= record_count; j++) {
          if (lengths[idx[j]] > lengths[idx[i]]) {
            tmp = idx[i]
            idx[i] = idx[j]
            idx[j] = tmp
          }
        }
      }

      limit = (top_n < record_count) ? top_n : record_count
      for (i = 1; i <= limit; i++) {
        k = idx[i]
        print headers[k]
        print seqs[k]
      }
    }
  ' "${fasta_files[@]}" > "$output_fasta" || die "Failed to build reduced query FASTA: $output_fasta"

  [ -s "$output_fasta" ] || die "Failed to build reduced query FASTA: $output_fasta"
}

log_run_parameters() {
  local script_name="$1"
  local assembly_mode="$2"
  local assembly_run_dir="$3"
  local query_fasta="$4"
  local query_sources="$5"
  local blast_db="$6"
  local top_n="$7"
  local threads="${8}"
  local max_target_seqs="${9}"
  local min_identity="${10}"
  local min_qcov="${11}"
  local evalue="${12}"

  cat <<EOF
=== BLAST related-species run ===
Script           : $script_name
Assembly mode    : $assembly_mode
Assembly run dir : $assembly_run_dir
Query FASTA      : $query_fasta
Query sources    : $query_sources
BLAST database   : $blast_db
Top sequences    : $top_n
Threads          : $threads
Max targets      : $max_target_seqs
Max HSPs         : 1
Min identity     : $min_identity
Min qcov         : $min_qcov
E-value          : $evalue
EOF
}

summarize_blast_hits() {
  local results_tsv="$1"
  local summary_tsv="$2"

  awk -F'\t' '
    {
      species = $8
      if (species == "" || species == "N/A") species = $9
      if (species == "" || species == "N/A") species = $2
      if (species == "" || species == "N/A") species = "Unknown"
      count[species]++
      if (!(species in best_bitscore) || $6 > best_bitscore[species]) best_bitscore[species] = $6
    }
    END {
      print "species\thit_count\tbest_bitscore"
      for (s in count) print s "\t" count[s] "\t" best_bitscore[s]
    }
  ' "$results_tsv" | sort -t$'\t' -k2,2nr -k3,3nr > "$summary_tsv"
}
