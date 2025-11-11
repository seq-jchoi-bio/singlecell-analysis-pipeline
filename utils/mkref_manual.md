# Building Reference Genomes for Non-Human/Mouse Species

This guide explains how to prepare reference genomes for species other than human and mouse using `mkref`.  
**Required files:**  
- Genome FASTA (`.fa`) file  
- Gene annotation GTF (`.gtf`) file  

**Important:**  
The chromosome names in both `.fa` and `.gtf` **must match exactly**. This document includes steps to verify and modify them.

For detailed command-line options, refer to:

```bash
mkref --help
# or
mkref -h
```

1. Check and prepare the FASTA file

- Use the command below to list chromosome headers in the FASTA file:
```bash
$ head rice.fa
$ grep "^>" rice.fa | sed 's/^>//' | cut -d' ' -f1 | sort
Chr1
Chr10
Chr11
Chr12
Chr2
Chr3
Chr4
Chr5
Chr6
Chr7
Chr8
Chr9
ChrSy
ChrUn
```

- Select only the desired chromosomes and adjust naming styles if needed.
This typically includes removing unplaced contigs (`ChrUn` or `ChrSy`), mitochondrial (`ChrM`), and chloroplast (`ChrPt`) sequences.

- Ensure the final chromosome naming convention (e.g., chr1, chr2, ..., chr12) matches that of the GTF file.

```bash
# Change Chr to chr
$ sed 's/^>Chr/>chr/' rice.fa > rice_chrstyle.fa

# Check the changed file
$ grep "^>" rice_chrstyle.fa | sed 's/^>//' | cut -d' ' -f1 | sort # check
```

2. Check and prepare the GTF file

- Use the following command to list chromosome names used in the GTF:
```bash
$ cut -f1 RAPDB_annot.gtf | sort | uniq
chr01
chr02
chr03
chr04
chr05
chr06
chr07
chr08
chr09
chr10
chr11
chr12
```

- **The GTF file must contain exactly 9 tab-separated fields per line.**
Extra tabs inside the 9th field (attributes) must be removed or modified:
```bash
$ cut -f9 RAPDB_annot.gtf | head
gene_id
gene_id
gene_id
gene_id
gene_id
gene_id
gene_id
gene_id
gene_id
gene_id

# Change chromosome style and modify the 9th field
awk 'BEGIN{OFS="\t"} {
  if ($1 ~ /^chr0[1-9]$/) sub(/^chr0/, "chr", $1)
  attr = ""
  for (i=9; i<=NF; i++) attr = attr $i " "
  $9 = attr
  NF = 9
  print
}' RAPDB_annot.gtf > RAPDB_annot_mkrep.gtf

# Check the result
$ cut -f9 RAPDB_annot_mkrep.gtf | head
gene_id "Os01g0100100"; transcript_id "Os01t0100100-01"; 
gene_id "Os01g0100100"; transcript_id "Os01t0100100-01"; 
gene_id "Os01g0100100"; transcript_id "Os01t0100100-01"; 
gene_id "Os01g0100100"; transcript_id "Os01t0100100-01"; 
gene_id "Os01g0100100"; transcript_id "Os01t0100100-01"; 
gene_id "Os01g0100100"; transcript_id "Os01t0100100-01"; 
gene_id "Os01g0100100"; transcript_id "Os01t0100100-01"; 
gene_id "Os01g0100100"; transcript_id "Os01t0100100-01"; 
gene_id "Os01g0100100"; transcript_id "Os01t0100100-01"; 
gene_id "Os01g0100100"; transcript_id "Os01t0100100-01";

$ cut -f1 RAPDB_annot_mkrep.gtf | sort | uniq
chr1
chr10
chr11
chr12
chr2
chr3
chr4
chr5
chr6
chr7
chr8
chr9
```

3. Create the configuration YAML
- Prepare a config.yaml file required by `mkref` using text editors (`nano`, `vim`):

```yaml
organism: "Oryza sativa"
genome: ["IRGSP_1_0"]
input_fasta: ["/data/shared/SativaGenome/rice_chrstyle.fa"] # Fasta absolute path
input_gtf: ["/data/shared/gtf/os_IRGSP/RAPDB_annot_mkrep.gtf"] # gtf absolute path
non_nuclear_contigs: ["chrSy", "chrUn"] # OPTIONAL!
```

* Note: Genome names must consist of letters, digits, hyphens (-), or underscores (_).
Dots (.) are not allowed.

4. Run mkref
- **No Conda environment activation is required.**
- Run mkref via the official cellranger-atac binary:

```bash
$ cd /data/programs/cellranger/atac-2.2.0
$ ./cellranger-atac mkref --config=IRGSP.config.yaml

# output example
$ ./cellranger-atac mkref --config=IRGSP.config.yaml
>>> Creating reference for IRGSP_1_0 <<<

Creating new reference folder at /data/programs/cellranger/atac-2.2.0/IRGSP_1_0
...done

Writing genome FASTA file into reference folder...
...done

Indexing genome FASTA file...
...done

Writing genes GTF file into reference folder...
...done

Writing genome metadata JSON file into reference folder...
Computing hash of genome FASTA file...
...done

Computing hash of genes GTF file...
...done

...done

Generating bwa index (may take over an hour for a 3Gb genome)...
[bwa_index] Pack FASTA... 1.46 sec
[bwa_index] Construct BWT for the packed sequence...
[BWTIncCreate] textLength=748942480, availableWord=64697916
[BWTIncConstructFromPacked] 10 iterations done. 99538208 characters processed.
[BWTIncConstructFromPacked] 20 iterations done. 190781872 characters processed.
[BWTIncConstructFromPacked] 30 iterations done. 271875488 characters processed.
[BWTIncConstructFromPacked] 40 iterations done. 343947696 characters processed.
[BWTIncConstructFromPacked] 50 iterations done. 408001632 characters processed.
[BWTIncConstructFromPacked] 60 iterations done. 464928944 characters processed.
[BWTIncConstructFromPacked] 70 iterations done. 515522048 characters processed.
[BWTIncConstructFromPacked] 80 iterations done. 560485312 characters processed.
[BWTIncConstructFromPacked] 90 iterations done. 600444752 characters processed.
[BWTIncConstructFromPacked] 100 iterations done. 635956784 characters processed.
[BWTIncConstructFromPacked] 110 iterations done. 667515936 characters processed.
[BWTIncConstructFromPacked] 120 iterations done. 695561744 characters processed.
[BWTIncConstructFromPacked] 130 iterations done. 720484912 characters processed.
[BWTIncConstructFromPacked] 140 iterations done. 742632656 characters processed.
[bwt_gen] Finished constructing BWT in 144 iterations.
[bwa_index] 163.74 seconds elapse.
[bwa_index] Update BWT... 1.46 sec
[bwa_index] Pack forward-only FASTA... 0.97 sec
[bwa_index] Construct SA from BWT and Occ... 103.99 sec
[main] Version: 0.7.17-r1198-dirty
[main] CMD: bwa index /data/programs/cellranger/atac-2.2.0/IRGSP_1_0/fasta/genome.fa
[main] Real time: 275.570 sec; CPU: 271.627 sec
done

Writing TSS and transcripts bed file...
...done

Writing genome metadata JSON file into reference folder...
Computing hash of genome FASTA file...
...done

Computing hash of genes GTF file...
...done

...done

>>> Reference successfully created at IRGSP_1_0 <<<
```