#!/bin/bash -l

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
BATCH_ROOT="$SCRIPT_DIR/batch_scripts"
LOG_DIR="$SCRIPT_DIR/output_slurm/submissions"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/submit_$(date +%Y%m%d_%H%M%S).log"

DRY_RUN=0
SUBMIT_CMD="sbatch"
DRY_COUNTER=1000

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    --submit-cmd)
      SUBMIT_CMD="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

log_submission() {
  local job_label="$1"
  local job_id="$2"
  local dep_spec="$3"
  local script_path="$4"
  printf '%s\t%s\t%s\t%s\t%s\n' "$(date +%Y-%m-%dT%H:%M:%S)" "$job_label" "$job_id" "$dep_spec" "$script_path" >> "$LOG_FILE"
}

format_dependency_arg() {
  local dep_ids=()
  for dep in "$@"; do
    if [[ -n "$dep" ]]; then
      dep_ids+=("$dep")
    fi
  done

  if [[ ${#dep_ids[@]} -eq 0 ]]; then
    return 0
  fi

  local joined
  joined=$(IFS=, ; echo "${dep_ids[*]}")
  printf '%s' "--dependency=afterok:${joined}"
}

submit_job() {
  local job_label="$1"
  local script_rel="$2"
  shift 2

  local script_path="$BATCH_ROOT/$script_rel"
  if [[ ! -f "$script_path" ]]; then
    echo "Missing batch script: $script_path" >&2
    exit 1
  fi

  local dep_arg
  dep_arg=$(format_dependency_arg "$@")

  local sbatch_output
  local job_id

  if [[ "$DRY_RUN" -eq 1 ]]; then
    DRY_COUNTER=$((DRY_COUNTER + 1))
    job_id="$DRY_COUNTER"
    if [[ -n "$dep_arg" ]]; then
      echo "DRY-RUN: $SUBMIT_CMD $dep_arg $script_path"
    else
      echo "DRY-RUN: $SUBMIT_CMD $script_path"
    fi
  else
    if [[ -n "$dep_arg" ]]; then
      sbatch_output=$($SUBMIT_CMD "$dep_arg" "$script_path")
    else
      sbatch_output=$($SUBMIT_CMD "$script_path")
    fi
    job_id=$(sed -n 's/^Submitted batch job \([0-9]\+\)$/\1/p' <<< "$sbatch_output")

    if [[ -z "$job_id" ]]; then
      echo "Could not parse job ID for $job_label" >&2
      echo "$sbatch_output" >&2
      exit 1
    fi
  fi

  log_submission "$job_label" "$job_id" "${dep_arg:-none}" "$script_path"
  echo "$job_id"
}

echo "Logging submissions to $LOG_FILE"
printf 'timestamp\tjob_label\tjob_id\tdependencies\tscript_path\n' > "$LOG_FILE"

job_preprocess_illumina=$(submit_job "preprocess_illumina" "01_preprocessing/preprocess_Illumina.batch")
job_rna_fastqc_trim=$(submit_job "rna_fastqc_trim" "01_preprocessing/rna_fastqc_trim.batch")
job_assembly_canu=$(submit_job "assembly_canu" "02_assembly/genome_assembly_PacBio.batch")
job_assembly_spades=$(submit_job "assembly_spades" "02_assembly/genome_assembly_SPAdes.batch" "$job_preprocess_illumina")

job_quast_spades=$(submit_job "quast_spades" "03_evaluation/genome_evaluation_QUAST_SPAdes.batch" "$job_assembly_spades")
job_busco_spades=$(submit_job "busco_spades" "03_evaluation/genome_evaluation_BUSCO_SPAdes.batch" "$job_assembly_spades")
job_prokka_spades=$(submit_job "prokka_spades" "04_annotation/genome_annotation_PROKKA_SPAdes.batch" "$job_assembly_spades")

job_quast_canu=$(submit_job "quast_canu" "03_evaluation/genome_evaluation_QUAST_Canu.batch" "$job_assembly_canu")
job_busco_canu=$(submit_job "busco_canu" "03_evaluation/genome_evaluation_BUSCO_Canu.batch" "$job_assembly_canu")
job_prokka_canu=$(submit_job "prokka_canu" "04_annotation/genome_annotation_PROKKA_Canu.batch" "$job_assembly_canu")

job_blast_spades=$(submit_job "blast_spades" "05_comparative_genomics/genome_related_species_BLAST_SPAdes.batch" "$job_prokka_spades")
job_mummer_spades=$(submit_job "mummer_spades" "05_comparative_genomics/genome_synteny_MUMmerplot_SPAdes.batch" "$job_blast_spades")

job_blast_canu=$(submit_job "blast_canu" "05_comparative_genomics/genome_related_species_BLAST_Canu.batch" "$job_prokka_canu")
job_mummer_canu=$(submit_job "mummer_canu" "05_comparative_genomics/genome_synteny_MUMmerplot_Canu.batch" "$job_blast_canu")

job_rna_align_spades=$(submit_job "rna_align_spades" "06_rnaseq/rna_align_spades_conditions.batch" "$job_rna_fastqc_trim" "$job_assembly_spades")
job_rna_align_canu=$(submit_job "rna_align_canu" "06_rnaseq/rna_align_canu_conditions.batch" "$job_rna_fastqc_trim" "$job_assembly_canu")
job_rna_counts_htseq=$(submit_job "rna_counts_htseq" "06_rnaseq/rna_read_counts_htseq.batch" "$job_rna_align_spades" "$job_rna_align_canu" "$job_prokka_spades" "$job_prokka_canu")

cat <<EOF
Submitted pipeline jobs:
  preprocess_illumina: $job_preprocess_illumina
  rna_fastqc_trim:     $job_rna_fastqc_trim
  assembly_canu:       $job_assembly_canu
  assembly_spades:     $job_assembly_spades
  quast_spades:        $job_quast_spades
  busco_spades:        $job_busco_spades
  prokka_spades:       $job_prokka_spades
  quast_canu:          $job_quast_canu
  busco_canu:          $job_busco_canu
  prokka_canu:         $job_prokka_canu
  blast_spades:        $job_blast_spades
  mummer_spades:       $job_mummer_spades
  blast_canu:          $job_blast_canu
  mummer_canu:         $job_mummer_canu
  rna_align_spades:    $job_rna_align_spades
  rna_align_canu:      $job_rna_align_canu
  rna_counts_htseq:    $job_rna_counts_htseq

Submission log:
  $LOG_FILE
EOF
