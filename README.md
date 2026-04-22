## About This Fork

This repository is a fork of [CIRI3](https://github.com/gyjames/CIRI3) that adds four decoupled sub-commands — **SCAN1**, **BUILD_UNIVERSE**, **SCAN2**, and **FINALIZE** — designed for high-throughput processing of large cohorts. In the original pipeline, all samples must be processed together in a single JVM invocation. The decoupled sub-commands split the pipeline so that per-sample stages can run independently (and in parallel across a cluster), while the two joint stages run once after all per-sample jobs complete.

---

## High-Throughput Processing

```
[SCAN1]          per sample — first-pass BSJ detection
[BUILD_UNIVERSE] joint      — build shared circRNA universe from all SCAN1 outputs
[SCAN2]          per sample — second-pass FSJ counting against the universe
[FINALIZE]       joint      — merge results and write expression matrices
```

### Sub-command reference

#### `SCAN1` — per-sample first scan

```
java -jar CIRI3_decoupled.jar SCAN1 \
    -I  sample.sam          \  # SAM/BAM; for STAR: Chimeric.out.junction,Aligned.out.sam,bwa.sam
    -O  outdir/sample       \  # output prefix; BSJ files and .scan1_meta written here
    -F  ref.fa              \
   [-A  ref.gtf]            \  # annotation (optional; required for -It 1)
   [-T  8]                  \  # threads (default: 1)
   [-Ma 0]                  \  # 0 = BWA-MEM (default), 1 = STAR
   [-S  0]                     # strigency filter applied during scan (0/1/2)
```

Outputs written to `{output_prefix}.*`:
- `{output_prefix}.scan1_meta` — checkpoint metadata read by downstream stages
- `{output_prefix}BSJ1` … `{output_prefix}BSJ{N}` — per-thread BSJ fragment files

> **Note**: BSJ files are written next to the output prefix, **not** next to the input SAM file.
> SCAN2 and FINALIZE locate them via `bsjPrefix=` in the `.scan1_meta` file.

#### `BUILD_UNIVERSE` — joint universe construction

```
java -jar CIRI3_decoupled.jar BUILD_UNIVERSE \
    -I  samples_scan1.tsv   \  # two-column TSV: samFile<TAB>scan1_meta_path
    -F  ref.fa              \
    -O  universe/cohort        # output prefix; writes {prefix}.universe
```

`samples_scan1.tsv` format (one row per sample):
```
/abs/path/sample1/bwa.sam   /abs/path/sample1/sample1.scan1_meta
/abs/path/sample2/bwa.sam   /abs/path/sample2/sample2.scan1_meta
```

#### `SCAN2` — per-sample second scan

```
java -jar CIRI3_decoupled.jar SCAN2 \
    -I  sample.sam          \  # same input triple as SCAN1
    -CU cohort.universe     \  # universe file from BUILD_UNIVERSE
    -SM sample.scan1_meta   \  # scan1_meta for this sample (locates BSJ files)
    -O  outdir/sample       \  # output prefix; writes {prefix}.fsj_counts
    -F  ref.fa              \
   [-T  8]                  \  # threads (should match SCAN1 thread count)
   [-Ma 0]                     # must match -Ma used in SCAN1
```

> **`-SM`** is strongly recommended. It tells SCAN2 where the BSJ files from SCAN1 are
> (via `bsjPrefix=` in the meta file). If omitted, SCAN2 falls back to looking for BSJ
> files adjacent to the input SAM file (backwards-compatible behaviour).

#### `FINALIZE` — joint merge and matrix writing

```
java -jar CIRI3_decoupled.jar FINALIZE \
    -I  finalize_samples.tsv \  # five-column TSV (see below)
    -CU cohort.universe      \  # universe file from BUILD_UNIVERSE (recommended)
    -F  ref.fa               \
    -O  finalize/result      \  # output prefix; writes .BSJ_Matrix and .FSJ_Matrix
   [-A  ref.gtf]             \  # annotation (optional)
   [-S  0]                   \  # strigency (0/1/2)
   [-E  0]                   \  # rel_exp threshold (default: 0)
   [-It 0]                      # intronic circRNAs (0/1, default: 0)
```

`finalize_samples.tsv` format — **five columns**, one row per sample:
```
/abs/bwa.sam  /abs/sample1.fsj_counts  8  sample1  /abs/sample1/sample1
/abs/bwa.sam  /abs/sample2.fsj_counts  8  sample2  /abs/sample2/sample2
```

| Column | Content |
|--------|---------|
| 1 | Absolute path to the input SAM file |
| 2 | Absolute path to the `.fsj_counts` file from SCAN2 |
| 3 | `fileSplitNum` value (read from `.scan1_meta`) |
| 4 | Sample name (used as column header in output matrices) |
| 5 | BSJ prefix — same value as `-O` passed to SCAN1/SCAN2 (locates BSJ files) |

### Intermediate file formats

**`.scan1_meta`** — written by SCAN1, read by BUILD_UNIVERSE, SCAN2 (`-SM`), and FINALIZE
```
samFile=/abs/path/to/bwa.sam
readLen=151
readNum=10503264
fileSplitNum=8
bsjPrefix=/abs/path/sample1/sample1
```

**`.universe`** — written by BUILD_UNIVERSE, read by SCAN2 (`-CU`) and FINALIZE (`-CU`)
```
seqLen=139
chr1	1379083	1414681
chr1	1384101	1397537
```
First line: `seqLen=<int>` (global max readLen − 12). Subsequent lines: one circRNA candidate per line (`chr start end`).

**`.fsj_counts`** — written by SCAN2, read by FINALIZE
```
chr1	1379083	1414681	3
chr1	1384101	1397537	0
```
Tab-separated: `chr start end fsjCount` — one line per circRNA in the universe.

---

### Complete example (STAR-based cohort)

Adjust paths at the top then run each block as a separate job (or sequentially).

```bash
# ── configuration ────────────────────────────────────────────────────────────
THREADS=64
REF_FASTA=/path/to/GRCh38.fa
GTF_FILE=/path/to/gencode.v32.annotation.gtf
STARDIR=/path/to/star_outputs          # contains STAR_output_<SAMPLE>/ sub-dirs
DECOUPLED_JAR=/path/to/CIRI3_decoupled.jar
WORKDIR=/path/to/workdir

SAMPLES=(
  "Div_100_S91"
  "Div_101_S92"
  "PARDOS_1_S1"
  "PARDOS_2_S2"
)

# ── SCAN1  (one job per sample) ───────────────────────────────────────────────
> "${WORKDIR}/samples_scan1.tsv"

for SAMPLE in "${SAMPLES[@]}"; do
  CHIMERIC="${STARDIR}/STAR_output_${SAMPLE}/Chimeric.out.junction"
  ALIGNED="${STARDIR}/STAR_output_${SAMPLE}/Aligned.out.sam"
  BWA_SAM="${STARDIR}/STAR_output_${SAMPLE}/bwa.sam"

  mkdir -p "${WORKDIR}/${SAMPLE}"
  OUT_PREFIX="${WORKDIR}/${SAMPLE}/${SAMPLE}"
  META="${OUT_PREFIX}.scan1_meta"

  java -jar "${DECOUPLED_JAR}" SCAN1 \
    -I "${CHIMERIC},${ALIGNED},${BWA_SAM}" \
    -O "${OUT_PREFIX}" \
    -F "${REF_FASTA}" \
    -A "${GTF_FILE}" \
    -Ma 1 -S 2 -T "${THREADS}"

  echo -e "${BWA_SAM}\t${META}" >> "${WORKDIR}/samples_scan1.tsv"
done

# ── BUILD_UNIVERSE  (single joint job) ───────────────────────────────────────
mkdir -p "${WORKDIR}/universe"
java -jar "${DECOUPLED_JAR}" BUILD_UNIVERSE \
  -I "${WORKDIR}/samples_scan1.tsv" \
  -F "${REF_FASTA}" \
  -O "${WORKDIR}/universe/cohort"

UNIVERSE="${WORKDIR}/universe/cohort.universe"

# ── SCAN2  (one job per sample) ───────────────────────────────────────────────
> "${WORKDIR}/finalize_samples.tsv"

for SAMPLE in "${SAMPLES[@]}"; do
  CHIMERIC="${STARDIR}/STAR_output_${SAMPLE}/Chimeric.out.junction"
  ALIGNED="${STARDIR}/STAR_output_${SAMPLE}/Aligned.out.sam"
  BWA_SAM="${STARDIR}/STAR_output_${SAMPLE}/bwa.sam"
  OUT_PREFIX="${WORKDIR}/${SAMPLE}/${SAMPLE}"
  META="${OUT_PREFIX}.scan1_meta"
  SPLIT_NUM=$(grep "^fileSplitNum=" "${META}" | cut -d= -f2)

  java -jar "${DECOUPLED_JAR}" SCAN2 \
    -I "${CHIMERIC},${ALIGNED},${BWA_SAM}" \
    -CU "${UNIVERSE}" \
    -SM "${META}" \
    -O "${OUT_PREFIX}" \
    -F "${REF_FASTA}" \
    -Ma 1 -T "${THREADS}"

  echo -e "${BWA_SAM}\t${OUT_PREFIX}.fsj_counts\t${SPLIT_NUM}\t${SAMPLE}\t${OUT_PREFIX}" \
    >> "${WORKDIR}/finalize_samples.tsv"
done

# ── FINALIZE  (single joint job) ─────────────────────────────────────────────
java -jar "${DECOUPLED_JAR}" FINALIZE \
  -I  "${WORKDIR}/finalize_samples.tsv" \
  -CU "${UNIVERSE}" \
  -F  "${REF_FASTA}" \
  -O  "${WORKDIR}/result" \
  -A  "${GTF_FILE}" \
  -S  2
```

**Output directory layout after a successful run:**

```
workdir/
  Div_100_S91/
    Div_100_S91.scan1_meta      ← SCAN1 checkpoint
    Div_100_S91BSJ1 … BSJ{N}   ← BSJ fragment files
    Div_100_S91.fsj_counts      ← SCAN2 output
  Div_101_S92/  …
  universe/
    cohort.universe             ← shared circRNA universe
  samples_scan1.tsv
  finalize_samples.tsv
  result.BSJ_Matrix             ← column headers = sample names
  result.FSJ_Matrix
  result.txt                    ← per-circRNA annotation
```

---

## Validation

Three regression scripts verify the decoupled pipeline produces matrices
equivalent to the joint `-W 1` pipeline:

- `scripts/test_decoupled_pipeline.sh` — small BWA test data bundled with
  the repo (no external data needed).
- `scripts/test_decoupled_pipeline_bwa_fullsize.sh` — full-size BWA data
  (edit `DATA_DIR` and `SAMPLES` at the top of the script).
- `scripts/test_decoupled_pipeline_star_fullsize.sh` — full-size STAR
  data (same pattern).

All three support `--threads N`, `--keep` (preserve the output
directory), and — on the full-size scripts — `--intron` (enable intron
mode `-It 1`) and `--use-current-joint` (use this repo's jar for the
joint run instead of the original `CIRI3_Java_*.jar`, for apples-to-apples
equivalence testing).

After running one of the full-size scripts with `--keep`, render the
validation report:

```bash
Rscript -e "rmarkdown::render('scripts/validation_report.Rmd')"
# or point at a different test output directory:
Rscript -e "rmarkdown::render('scripts/validation_report.Rmd', \
    params=list(out_root='/path/to/decoupled_comparison'))"
```

That produces
[`scripts/validation_report.html`](https://rushalz.github.io/CIRI3/scripts/validation_report.html)
— per-sample BSJ/FSJ scatter plots, cell-level agreement tables,
Pearson/Spearman correlations, the largest remaining disagreements, and
the per-stage benchmark (wall time, CPU%, peak RAM) collected by
`scripts/_bench.sh`.

---

## Building from Source

### Prerequisites

- JDK 8 or later (`javac`, `jar`)
- `lib/htsjdk-3.0.4.jar` (included in the repository)

### Quick build (recommended)

```bash
bash scripts/build_jar.sh
```

This compiles `src/` with `-source 8 -target 8`, bundles
`lib/htsjdk-3.0.4.jar`, and writes a runnable
`CIRI3_decoupled.jar` at the repository root. The jar is also checked
in — you only need to rebuild it if you change Java sources.

### Manual compile

```bash
mkdir -p bin
find src -name "*.java" > sources.txt
javac -source 8 -target 8 -cp "lib/htsjdk-3.0.4.jar" -d bin @sources.txt
rm sources.txt
```

### Manual distributable JAR

Create a manifest file that sets the main class and bundles the htsjdk dependency:

```bash
# Extract htsjdk into the bin directory so it is included in the JAR
cd bin
jar xf ../lib/htsjdk-3.0.4.jar
cd ..

# Write the manifest
cat > MANIFEST.MF <<'EOF'
Manifest-Version: 1.0
Main-Class: com.zx.test.TestParameters
EOF

# Build the JAR (name it after your Java version, e.g. CIRI3_Java_1.8.0.jar)
JAVA_VERSION=$(java -version 2>&1 | head -1 | awk -F'"' '{print $2}')
jar cfm "CIRI3_Java_${JAVA_VERSION}.jar" MANIFEST.MF -C bin .
rm MANIFEST.MF
```

The resulting `CIRI3_Java_<version>.jar` (or `CIRI3_decoupled.jar` from
the helper script) is self-contained and can be used as:

```bash
java -jar CIRI3_decoupled.jar -I sample.sam -O result -F ref.fa
java -jar CIRI3_decoupled.jar SCAN1 -I sample.sam -O sample -F ref.fa -T 8
```

---

<p align="center">
  <img src="./data/CIRI3_logo.png" width="60%">
</p>


## About

CIRI3 is a comprehensive analysis package for the detection and quantification of circRNAs in RNA-Seq data, while providing differential expression analysis of circRNAs at multiple levels.

## Installation

The CIRI3 was constructed based on java. The `environment.yaml` was provided and the dependencies can be installed as the follow:
```
git clone https://github.com/gyjames/CIRI3.git
cd CIRI3
conda env create -n CIRI3 -f ./environment.yaml
conda activate CIRI3
```

## Usage

### Step1. Mapping

#### 1)BWA
Recommended protocols for running BWA-MEM:
```
bwa index -a bwtsw ref.fa
bwa mem –T 19 ref.fa reads.fq > sample.sam (single-end reads)
bwa mem –T 19 ref.fa read1.fq read2.fq > sample.sam (paired-end reads)
```
* ref.fa: FASTA file of all reference sequences

#### 2)STAR
```
STAR --runThreadN 10 \
--genomeDir star_index \
--outSAMtype SAM \
--readFilesIn reads1.fq read2.fq \
--outFileNamePrefix output_dir/ \
--outReadsUnmapped Fastx \
--outSJfilterOverhangMin 15 12 12 12 \
--alignSJoverhangMin 15 \
--alignSJDBoverhangMin 15 \
--outFilterMultimapNmax 20 \
--outFilterScoreMin 1 \
--outFilterMatchNmin 1 \
--outFilterMismatchNmax 2 \
--chimSegmentMin 15 \
--chimScoreMin 15 \
--chimJunctionOverhangMin 15

cd output_dir/
bwa mem -T 19 ref.fa Unmapped.out.mate1 Unmapped.out.mate2 > bwa.sam
```

### Step2. CircRNA detection and quantification

CIRI3 provides multiple input options for the identification and quantification of circRNAs, including single-sample input, multiple-sample input, multiple-sample files containing RNase R treated information. In addition, users have the option to input a collection of circRNAs of interest, allowing CIRI3 to quantify the BSJ and FSJ of these circRNAs within the samples.

A small test data set is provided in [data/](data/). The unpacked files are:

#### 1)Single-SAM/BAM files as input
CIRI3 will give information about circRNAs in the sample.
```
# without the annotation gtf
java -jar CIRI3.jar -I ./data/circRNA/Single/sample.sam -O ./data/circRNA/Single/result.txt -F ./data/circRNA/ref.fa (BWA-MEM)
java -jar CIRI3.jar -I /path/to/Chimeric.out.junction,/path/to/Aligned.out.sam,/path/to/bwa.sam -O output.ciri3 -F ref.fa -Ma 1 (STAR)

# with the annotation gtf
Java -jar CIRI3.jar -I sample.sam -O result.txt -F ref.fa -A ref.gtf (BWA-MEM)
java -jar CIRI3.jar -I /path/to/Chimeric.out.junction,/path/to/Aligned.out.sam,/path/to/bwa.sam -O output.ciri3 -F ref.fa -A ref.gtf -Ma 1 (STAR)

```
* sample.sam and bwa.sam: SAM file generated by BWA-MEM
* Aligned.out.sam and Chimeric.out.junction: output files generated by STAR
* ref.gtf: GTF/GFF3 formatted annotation file

#### 2)Multiple-SAM/BAM files as input
CIRI3 will give information about circRNA in the samples, circRNA BSJ expression matrix and circRNA FSJ expression matrix.
```
# without the annotation gtf
java -jar CIRI3.jar -I ./data/circRNA/Mutiple/samples.tsv -O ./data/circRNA/Mutiple/result.txt -F ./data/circRNA/ref.fa -W 1 (BWA-MEM)
java -jar CIRI3.jar -I samples.tsv -O output.ciri3 -F ref.fa -Ma 1 -W 1 (STAR)

# with the annotation gtf
Java -jar CIRI3.jar -I samples.tsv -O result.txt -F ref.fa -A ref.gtf -W 1 (BWA-MEM)
java -jar CIRI3.jar -I samples.tsv -O output.ciri3 -F ref.fa -A ref.gtf -Ma 1 -W 1 (STAR)
```
* samples.tsv is a tab separated file like (BWA-MEM):
  + Column 1: absolute path to the sam/bam file
```
/path/to/sample1.sam
/path/to/sample2.sam
/path/to/sample3.sam
/path/to/sample4.sam
```
* samples.tsv is a tab separated file like (STAR):
  + Column 1: absolute path to the sam/bam file
```
/path/to/sample1/Chimeric.out.junction,/path/to/sample1/Aligned.out.sam,/path/to/sample1/bwa.sam
/path/to/sample2/Chimeric.out.junction,/path/to/sample2/Aligned.out.sam,/path/to/sample2/bwa.sam
/path/to/sample3/Chimeric.out.junction,/path/to/sample3/Aligned.out.sam,/path/to/sample3/bwa.sam
/path/to/sample4/Chimeric.out.junction,/path/to/sample4/Aligned.out.sam,/path/to/sample4/bwa.sam
```
#### 3)Multiple-SAM/BAM files containing RNase R treated information as input
CIRI3 gives information on circRNAs in the samples, circRNA BSJ expression matrix, circRNA FSJ expression matrix and circRNA enrichment ratios.
```
# without the annotation gtf
java -jar CIRI3.jar -I ./data/circRNA/RNase/sample_infor.tsv -O ./data/circRNA/RNase/result.txt -F ./data/circRNA/ref.fa -W 2

# with the annotation gtf
Java -jar CIRI3.jar -I sample_infor.tsv -O result.txt -F ref.fa -A ref.gtf -W 2
```
* samples.tsv is a tab separated file like:
  + Column 1: absolute path to the sam/bam file
  + Column 2: RNase R treated information. 0 (without RNase R treated)    1 (with RNase R treated)   2 (RNase R treatment unknown)
```
/path/to/sample1.sam  0
/path/to/sample2.sam  0
/path/to/sample3.sam  1
/path/to/sample4.sam  1
```

#### 4)CircRNA collections of interest and sam/bam files as inputs
CIRI3 will give the BSJ expression matrix and FSJ expression matrix of the input circRNA collection in the sample.
```
java -jar CIRI3.jar -I ./data/circRNA/User/sample.sam -UC ./data/circRNA/User/circRNA_List.txt -O ./data/circRNA/User/result.txt -F ./data/circRNA/ref.fa
```

* circRNA_List.txt is a tab separated file like:
  + #Column 1: chromosome
  + #Column 2: start loci of a circRNA on the chromosome(1-based)
  + #Column 3: end loci of a circRNA on the chromosome
```
chr1	1379083	1414681
chr1	1384101	1397537
```

### Step3. Differential Expression Analysis

CIRI3 provides three levels of differential expression analysis algorithms, including the expression of a single circRNA, the junction ratio of a single circRNA, and the relative expression of multiple circRNAs derived from a single gene. The entire differential expression analysis script incorporates parts of the code from the edgeR package in R or rMATS to calculate the statistical significance of differential expression (p-values).

#### 1)the expression of a single circRNA

##### Study without biological replicate(Input in two ways)

```
java -jar CIRI3.jar DE_BSJ -I ./data/DE/BSJ/Without_RE/infor.tsv -O ./data/DE/BSJ/Without_RE/result.txt

java -jar CIRI3.jar DE_BSJ -I ./data/DE/BSJ/Without_RE/infor.tsv -M ./data/DE/BSJ/Without_RE/BSJ_Matrix.txt -O ./data/DE/BSJ/Without_RE/result.txt
```

* infor.tsv is a tab separated file like:
  + Column 1: sample name
  + Column 2: absolute path to CIRI3 output result files
  + Column 3: sample type
  + Column 4: number of mapped reads (It can be obtained from the log file output by CIRI3)
```
Sample	Path	Class	Map_reads
Case1	path/to/case1_CIRI3.txt	Case	1050326
Control1	path/to/control1_CIRI3.txt	Control	1055952
```

* BSJ_Matrix.txt is a tab separated file like:
```
circRNA_ID	Case1	Control1
chr9:33886881|33900271	5	6
chr9:83311887|83343287	4	0
chr9:33318731|33338589	0	6
```

##### Study with biological replicate(Input in two ways)

```
java -jar CIRI3.jar DE_BSJ -I ./data/DE/BSJ/With_RE/infor.tsv -G ./data/DE/BSJ/With_RE/Gene_Expression.txt -O ./data/DE/BSJ/With_RE/result.txt

java -jar CIRI3.jar DE_BSJ -I ./data/DE/BSJ/With_RE/infor.tsv -G ./data/DE/BSJ/With_RE/Gene_Expression.txt -M ./data/DE/BSJ/With_RE/BSJ_Matrix.txt -O ./data/DE/BSJ/With_RE/result.txt
```

* infor.tsv is a tab separated file like:
  + Column 1: sample name
  + Column 2: absolute path to CIRI3 output result files
  + Column 3: sample type
  + Column 4: subject (optional, only for paired samples)
```
Sample	Path	Class	Num
Case1	path/to/case1_CIRI3.txt	Case	1
Case2	path/to/case2_CIRI3.txt	Case	2
Control1	path/to/control1_CIRI3.txt	Control	1
Control2	path/to/control2_CIRI3.txt	Control	2
```

* Gene_Expression.txt is a tab separated file like Gene_Expression.txt is a tab separated file like (which can be obtained by running featureCounts on a BWA-aligned and sorted BAM files):
```
Geneid	Case1	Case2	Control1	Control2
ENSG00000236875.3	14	4	13	4
ENSG00000181404.18	59	21	37	10
```

#### 2)the junction ratio of a single circRNA

```
java -jar CIRI3.jar DE_Ratio -I ./data/DE/Ratio/infor.tsv -BM ./data/DE/Ratio/BSJ_Matrix.txt -FM ./data/DE/Ratio/FSJ_Matrix.txt -O ./data/DE/Ratio/result.txt
```

* infor.tsv is a tab separated file like:
  + Column 1: sample name
  + Column 2: sample type
```
Sample	Class
Case1	Case
Case2	Case
Control1	Control
Control2	Control
```

* FSJ_Matrix.txt is a tab separated file like:
```
circRNA_ID	Case1	Case2	Control1	Control2
chr5:171214|173067	15	0	12	0
chr1:18854213|18857631	5	6	0	0
```

#### 3)the relative expression of multiple circRNAs derived from a single gene(Input in two ways)

```
java -jar CIRI3.jar DE_Relative -I ./data/DE/Relative_Expression/infor.tsv -O ./data/DE/Relative_Expression/result.txt

java -jar CIRI3.jar DE_Relative -I ./data/DE/Relative_Expression/infor.tsv -M ./data/DE/Relative_Expression/BSJ_Matrix.txt -GC ./data/DE/Relative_Expression/circ_Gene.txt -O ./data/DE/Relative_Expression/result.txt
```

* infor.tsv is a tab separated file like:
  + Column 1: sample name
  + Column 2: absolute path to CIRI3 output result files
  + Column 3: sample type
```
Sample	Path	Class
Case1	path/to/case1_CIRI3.txt	Case
Case2	path/to/case2_CIRI3.txt	Case
Control1	path/to/control1_CIRI3.txt	Control
Control2	path/to/control2_CIRI3.txt	Control
```

* circ_Gene is a tab separated file like:
```
circRNA_ID	gene_id
chr9:465962|467219	ENSG00000227155.7
chr9:846960|894195	ENSG00000137090.12
chr9:2161686|2181676	ENSG00000080503.24
```

### All Arguments

```
java -jar CIRI3.jar -H

Arguments:

    -I, --in
          Input SAM file name or SAM files list(required; generated by BWA-MEM)
    -O, --out
          Output circRNA file name(required)
    -F, --ref_file
          FASTA file of all reference sequences. Please make sure this file is
          the same one provided to BWA-MEM (required).
    -A, --anno
          input GTF/GFF3 formatted annotation file name (optional)
    -G, --log
          output log file name (optional)
    -H, --help
          show this help information
    -Max, --max_span
          max spanning distance of circRNAs (default: 200000)
    -Min, --min_span
          min spanning distance of circRNAs (default: 140)
    -S, --strigency
          2: only output circRNAs supported by more than 2 distinct PCC signals (default)
          1: only output circRNAs supported by more than 2 junction reads
          0: output all circRNAs regardless junction read or PCC signal counts
    -U, --mapq_uni
          set threshold for mappqing quality of each segment of junction reads
          (default: 10; should be within [0,30])
    -E, --rel_exp
          set threshold for relative expression calculated based on counts of
          junction reads and non-junction reads (optional: e.g. 0.1)
    -Mc, --mitochondria
          0: Skip the recognition of  mitochondria circRNA (default)
          1: Perform the recognition of mitochondria circRNA (the ID of mitochondrion in reference file is required)
    -M, --chrM
          tell CIRI3 the ID of mitochondrion in reference file(s) (default:
          chrM)
    -T, --thread_num
          set number of threads for parallel running (default: 1)
    -It, --intron
          0: Skip the recognition of  intronic self-ligated circRNA (default)
          1: Perform the recognition of intronic self-ligated circRNA (formatted annotation file is required)
    -Sp, --splicing_signals
	    0: Only canonical GT-AG splice signals were considered (default)
	    1: Both canonical GT-AG and non-canonical splice signals were considered
    -Ma, --mapper\r\n"
	    0: The SAM file generated by the BWA-MEM (default)
	    1: the SAM file generated by STAR
    -W, --way
          0: Input a single SAM file to identify circRNA (default)
          1: Input a SAM files list (including the absolute path to SAM files
             and sample name) to identify circRNA
          2: Input a SAM files list (including the absolute path to SAM files,
             RNase information, and sample names) to identify circRNA, this can generate
             circRNA confidence score
    -UC, --user_circ
          file of circRNA collections of interest to users

```

```
java -jar CIRI3.jar DE_BSJ -H

Arguments:

    -I, --in
          the file of sample list(required)
    -O, --out
          output differential expression result(required)
    -M, --matrix
          circRNA BSJ expression matrix (optional)
    -G, --gene
          gene expression matrix
    -P, --pvalue
          p value threshold for DE score calculation (default: 0.05)
    -H, --help
          show this help information
```

```
java -jar CIRI3.jar DE_Ratio -H

Arguments:

    -I, --in
          the file of sample list(required)
    -O, --out
          output differential expression result(required)
    -BM, --bsj_matrix
          circRNA BSJ expression matrix (required)
    -FM, --fsj_matrix
          circRNA FSJ expression matrix (required)
    -T, --thread_num
          set number of threads for parallel running (default: 1)
    -H, --help
          show this help information
```

```
java -jar CIRI3.jar DE_Relative -H

Arguments:

    -I, --in
          the file of sample list(required)
    -O, --out
          output differential expression result(required)
    -M, --matrix
          circRNA BSJ expression matrix (optional)
    -GC, --circ_gene
          the gene information file corresponding to circRNAs (optional)
    -T, --thread_num
          set number of threads for parallel running (default: 1)
    -H, --help
          show this help information
```

## Output

The three main output files are:
* `result.txt` is an information file for predicted circRNA.
  + Column 1: ID of a predicted circRNA in the pattern of "chr:start|end";
  + Column 2: chromosome of a predicted circRNA
  + Column 3: start loci of a predicted circRNA on the chromosome
  + Column 4: end loci of a predicted circRNA on the chromosome
  + Column 5: circular junction read (also called as back-spliced junction read) count of a predicted circRNA 
  + Column 6: unique CIGAR types of a predicted circRNA. For example, a circRNAs have three junction reads: read A (80M20S, 80S20M), read B (80M20S, 80S20M), read 
              C (40M60S, 40S30M30S, 70S30M), then its has two SM types (80S20M, 70S30M), two MS types (80M20S, 70M30S) and one SMS type (40S30M30S). Thus its 
              SM_MS_SMS should be 2_2_1.
  + Column 7: non-junction read count of a predicted circRNA that mapped across the circular junction but consistent with linear RNA instead of being back-spliced
  + Column 8: ratio of circular junction reads calculated by 2*#junction_reads/(2*#junction_reads+#non_junction_reads). 
  + Column 9: type of a circRNA according to positions of its two ends on chromosome (exon, intron or intergenic_region; only available when annotation file is 
              provided)
  + Column 10: ID of the gene(s) where an exonic or intronic circRNA locates
  + Column 11: strand info of a predicted circRNAs (new in CIRI2)
  + Column 12: all of the circular junction read IDs (split by ",")
  + Column 13: the count of high-quality BSJ reads from the first scan
* `result.txt.BSJ_Matrix` is a BSJ expression matrix file for predicted circRNA.
  + Each detected circRNA is reported on a separate line.
* `result.txt.FSJ_Matrix` is a FSJ expression matrix file for predicted circRNA.
* `result.txt.Score` is a enrichment ratio information file for predicted circRNA.

---

## References

1. Li, Heng. "Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM." arXiv preprint arXiv:1303.3997 (2013).
2. Shen S, Park JW, Lu ZX, et al. rMATS: robust and flexible detection of differential alternative splicing from replicate RNA-Seq data. Proc Natl Acad Sci U S A. 2014;111(51):E5593-E5601. doi:10.1073/pnas.1419161111
3. Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics. 2010;26(1):139-140. doi:10.1093/bioinformatics/btp616
4. Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics. 2014;30(7):923-930. doi:10.1093/bioinformatics/btt656