#!/usr/bin/env bash
set -eo pipefail

usage() {
  echo "Usage: $0 -i <input_dir> -o <output_dir>" >&2
  exit 1
}

IN=""
OUT=""
while getopts "i:o:" opt; do
  case "$opt" in
    i) IN="$OPTARG" ;;
    o) OUT="$OPTARG" ;;
    *) usage ;;
  esac
done
[[ -z "${IN}" || -z "${OUT}" ]] && usage

FASTAD="${IN}/fasta"
GENESD="${IN}/genes"
FA="${FASTAD}/genome.fa"
GTFGZ="${GENESD}/genes.gtf.gz"

[[ -d "${FASTAD}" ]] || { echo "ERROR: Missing directory: ${FASTAD}" >&2; exit 2; }
[[ -d "${GENESD}" ]] || { echo "ERROR: Missing directory: ${GENESD}" >&2; exit 2; }
[[ -f "${FA}" ]]     || { echo "ERROR: Missing file: ${FA}" >&2; exit 2; }
[[ -f "${GTFGZ}" ]]  || { echo "ERROR: Missing file: ${GTFGZ}" >&2; exit 2; }

mkdir -p "${OUT}"
SIZES_OUT="${OUT}/genome.chrom.sizes"
GTF_OUT="${OUT}/genes.gtf"

if command -v samtools >/dev/null 2>&1; then
  TMPDIR="${OUT}/.tmp_make_ref_$$"
  mkdir -p "${TMPDIR}"
  trap 'rm -rf "${TMPDIR}"' EXIT

  cp "${FA}" "${TMPDIR}/genome.fa"
  samtools faidx "${TMPDIR}/genome.fa"
  cut -f1,2 "${TMPDIR}/genome.fa.fai" > "${SIZES_OUT}"
else
  awk '
    BEGIN{hdr="";len=0}
    /^>/{
      if(hdr!=""){print hdr "\t" len}
      hdr=substr($0,2); sub(/ .*/,"",hdr); len=0; next
    }
    {
      gsub(/[ \t\r\n]/,"");
      len+=length($0)
    }
    END{ if(hdr!=""){print hdr "\t" len} }
  ' "${FA}" > "${SIZES_OUT}"
fi

gzip -dc "${GTFGZ}" > "${GTF_OUT}"

echo "Wrote: ${SIZES_OUT}"
echo "Wrote: ${GTF_OUT}"
