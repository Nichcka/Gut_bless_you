# Gut_bless_you

A pipeline for RNA-Seq quality control, mapping, result aggregation, and basic statistical analysis.
___

## Quality control

### FastQC

Run FastQC for initial quality assessment of raw FASTQ files.

``` cd qc/fastqc ```

``` fastqc *.fastq.gz -o /path/to/results_folder ```

### MultiQC

Aggregate FastQC reports into a single summary with MultiQC.

``` cd ../multiqc ```

``` multiqc /path/to/result_fastqc_folder ```
___

## Mapping

### Salmon

Build a Salmon index and quantify transcript abundance for paired-end reads.

``` cd ../../salmon ```

for one assemble:

``` salmon index -t ./path/to/reference.fasta -i salmon_index -k 21 ```

``` salmon quant -i salmon_index/ -l A -1 /path/to/sample_R1.fastq.gz -2 /path/to/sample_R2.fastq.gz -o /path/to/result_folder --meta --minScoreFraction 0.25 --consensusSlack 0.2 --maxRecoverReadOcc 1000 ```

**or** use the helper script to process multiple samples **automatically**:

``` ./run_salmon_hc.sh /path/to/folder_with_*fastq.gz /path/to/results/ ```
___

## Making result tables

Aggregate Salmon quantification results with metadata into summary tables.

``` ./make_result_table.sh /path/to/salmon_results/ /path/to/SQ_loci_meta.csv /path/to/SraRunTable.tsv ```

___

## Statistics

### Wilcoxon test

Perform the Wilcoxon statistical test on TPM tables for comparative analysis.

``` python wilcoxon.py /path/to/all_tpm_ibs.tsv /path/to/all_tpm_hc.tsv ```




