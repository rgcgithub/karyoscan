# KaryoScan: abnormal karyotype detection from high-throughput whole-exome sequence data

## What KaryoScan is for

KaryoScan is a method for analyzing whole-exome sequence data to detect chromosomal anomalies including aneuploidies, partial chromosome events, mosaic events, and copy neutral events. KaryoScan detects abnormal karyotypes through analysis of read depths and SNP allele balance distributions.

KaryoScan is designed for use on large sample cohorts (e.g. at least in the hundreds or thousands) as it requires a model training step. For smaller sample sets, portions of the code may work (e.g. allele balance analysis) but with lower accuracy.

While KaryoScan can detect some cancer-related chromosomal anomalies, it is not intended for use in tumor-normal sequencing. Additionally, KaryoScan was not designed for whole-genome sequencing data; it would conceptually be applicable, but has not been tested.

Lastly, KaryoScan was developed as a complementary approach to our exome-based CNV detection method [CLAMMS](https://github.com/rgcgithub/clamms). The first part of this README mimics the relevant pre-processing steps from the CLAMMS tutorial required to generate the input read coverage tracks used by KaryoScan. At a high level, this includes computing read coverage over pre-defined exon "windows" (exons or 500-1000bp segments of large exons) along with annotation of windows for GC content, mappability statistics, and known problematic regions. Notably, KaryoScan analyzes the raw window-based coverage values, NOT the CLAMMS normalized coverage values (i.e. runs `annotate_windows.sh` & `sam_gatk_coverage_to_bed`, but not `normalize_coverage`).

## Related Publications

* KaryoScan methods paper: Maxwell EK, Gonzaga-Jauregui C, O'Dushlaine C, et al. (2017) KaryoScan: abnormal karyotype detection from whole-exome sequence data. bioRxiv.
* CLAMMS (CNV) papers:
    * method: Packer JS, Maxwell EK, O'Dushlaine C, et al. (2015) CLAMMS: a scalable algorithm for calling common and rare copy number variants from exome sequencing data. Bioinformatics 32 (1): 133-135.) [link](http://bioinformatics.oxfordjournals.org/content/32/1/133) describes the methods of CLAMMS, as well as the results of validation experiments we used to evaluate its performance in comparison to previous tools.
    * application to 50K DiscovEHR Study samples: Maxwell EK, Packer JS, O'Dushlaine C, McCarthy SE, Hare-Harris A, Gonzaga-Jauregui C, et al. (2017) Profiling copy number variation and disease associations from 50,726 DiscovEHR Study exomes. bioRxiv. [http://biorxiv.org/content/early/2017/03/22/119461](http://biorxiv.org/content/early/2017/03/22/119461)

## Software Requirements

KaryoScan is designed for use within a Linux environment. The following software packages and libraries are required:

* bash
* Python
* R
    * R plotting libraries (optional, but highly recommended for diagnostic plotting)
        * ggplot2
        * grid
        * gridBase
        * gridExtra
        * stringr
        * reshape2
* bcftools >= v1.3 (requires BCFtools/RoH module)
    * https://github.com/samtools/bcftools
* bedtools2 (requires subtract module)
    * https://github.com/arq5x/bedtools2
* datamash
    * https://www.gnu.org/software/datamash/
* CLAMMS pre-processing scripts

## Read coverage pre-processing with CLAMMS (from [CLAMMS GitHub](https://github.com/rgcgithub/clamms))

First, clone the CLAMMS Github repository and compile the code:

    git clone https://github.com/rgcgithub/clamms.git
    cd clamms
    make

Set the environment variable `CLAMMS_DIR` to the appropriate path using the `export` command

    export $CLAMMS_DIR=/path/to/clamms

Now you will need to generate a file `windows.bed`. This file will list windows along the exome for which CLAMMS will estimate copy numbers, along with metadata for those windows. Most windows will simply be exons from your exome capture design, but large exons (>= 1000 bp) will be split up into equally-sized calling windows of ~500 bp.

To generate `windows.bed`, you will need four input files:

1. targets.bed — a BED file listing your exome capture regions.
1. genome.fa — an indexed FASTA file for the reference genome you are using.
1. mappability.bed — a BED file listing mappability scores across the genome. More detail on this below.
1. clamms_special_regions.bed — provided in the data/ directory with the code distribution (hg19 coordinates).

The chromosome names in the BED files and in the genome index should not have "chr" preceding the number/letter (i.e. "1" instead of "chr1"). The BED files must be sorted using either `bedtools sort` or `sort -k1,1 -k2,2n`.

The FASTA index should be generated from the raw FASTA file using [BWA](http://bio-bwa.sourceforge.net/bwa.shtml): `bwa index genome.fa`.

The mappability score for a given base is one divided by the number of locations in the genome that the k-mer starting at that base aligns to (k = the length of your reads), with up to two mismatches allowed (see [here](http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability) for more details). You can download mappability tracks for 75-mers or 100-mers on the GRCh37 human reference genome from the link above and convert them to CLAMMS-ready BED files (requires `bigWigToWig` tool from [UCSC](http://genome.ucsc.edu/goldenpath/help/bigWig.html)):

    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign75mer.bigWig
    bigWigToWig wgEncodeCrgMapabilityAlign75mer.bigWig wgEncodeCrgMapabilityAlign75mer.wig
    grep -v '^#' wgEncodeCrgMapabilityAlign75mer.wig | sed 's/^chr//g' > mappability.bed

Once you have the input files ready, you can generate `windows.bed` with the following commands. This will take ~5 minutes. Note that the preprocessing script `annotate_windows.sh` requires [Bedtools](http://bedtools.readthedocs.org) to be installed and in your system PATH.

    export INSERT_SIZE=200
    chmod +x $CLAMMS_DIR/annotate_windows.sh
    $CLAMMS_DIR/annotate_windows.sh targets.bed genome.fa mappability.bed $INSERT_SIZE $CLAMMS_DIR/data/clamms_special_regions.bed >windows.bed

The `INSERT_SIZE` variable should be set to a value that is a little bit larger than the average insert size for your sequencing process (so that most reads will come from inserts of size <= this value). For example, we use `INSERT_SIZE = 200` when our mean insert size is ~150 bp. If a window is smaller than `INSERT_SIZE`, it is extended to the length of `INSERT_SIZE` for purposes of calculating it's GC content. This is because according to [Benjamini and Speed (2012)](http://www.ncbi.nlm.nih.gov/pubmed/22323520), GC coverage bias is best estimated based on the GC content of the insert, not the individual reads.

#### Troubleshooting windows.bed file generation
If you have trouble generating the windows.bed file, your input files are likely improperly formatted or have inconsistencies. A few things you should check:

* All BED files are sorted properly, using `sort -k1,1 -k2,2n`. This sorts by chromosome name (string sort) then by start position (numeric sort). You should re-sort all files in the event that your system locale settings differ from those that were used to sort the externally sourced input files.
* Chromosome naming consistency: Make sure that all chromosomes are named consistently (i.e. chromosome 1 is "1", not "chr1"). This must be the case in all input files, including the genome.fa input file.

Here is a simple test you can run on your input files to make sure they are consistent and sorted properly:

    cut -f 1 targets.bed | uniq
    cut -f 1 mappability.bed | uniq
    cut -f 1 clamms_special_regions.bed | uniq

All of these should return the same chromosome names and sort order:

    1
    10
    ...
    19
    2
    20
    21
    22
    3
    4
    ...
    9
    X
    Y

Also check the chromosome names in the genome FASTA file (sort order is not important):

    grep '^>' -m 24 genome.fa
    >1
    >2
    ...
    >22
    >X
    >Y

#### Computing depths of coverage

You will need a BED file for each of your samples listng the mean depth of coverage at each of the exact windows listed in `windows.bed`. The coverage files must be named in the following format: `sample_name.coverage.bed`

The sample depth-of-coverage files can be generated from BAM files using [Samtools](http://www.htslib.org/):

    # 30 = minimum mapping quality for a read to be counted
    samtools bedcov -Q 30 windows.bed sample.bam \
    | awk '{ printf "%s\t%d\t%d\t%.6g\n", $1, $2, $3, $NF/($3-$2); }' \
    >sample.coverage.bed

They can also be generated from GATK DepthOfCoverage output files:

    $CLAMMS_DIR/sam_gatk_coverage_to_bed sample.gatk_readDepth_1x_q30.out windows.bed >sample.coverage.bed

This step will almost certainly take longer than the CNV calling process itself. To speed it up, consider processing several samples in parallel using `xargs`:

    cat list.of.samples.txt | xargs -P $NUM_PROCESSES --max-args 1 ./compute_coverage.sh

Where `list.of.samples.txt` lists each sample name (one per line) and `compute_coverage.sh` is a shell script that takes a sample name as its argument and generates its coverage file using one of the two methods shown above.

## Getting started with KaryoScan

Clone the repository

    git clone https://github.com/rgcgithub/karyoscan

Generate a 'windows.bed' file (from CLAMMS, see section Read coverage pre-processing with CLAMMS). This file contains exon target windows annotated with GC-content and mappability, which will be used to select high-confidence windows (mid-range GC content, high mappability) to pool read coverage from.

Generate sample coverage statistics across pre-computed windows (from CLAMMS, see section Read coverage pre-processing with CLAMMS).

Prep single-sample VCF files. Each should have a Tabix index, or alternatively, 'tabix' should be installed and in your system PATH to generate one on the fly. Each VCF file should contain DP (total read depth) and AD (allele depth) flags, or equivalent.

Generate/customize the config file.

## Config file specifications

The KaryoScan config file sets several universal parameters for use within the analysis pipeline and is passed as the first argument to the main driver script. The parameters defined in the config file include file/binary paths and variables used within the KaryoScan pipeline.

The config file is a tab-separated text file with at least 2 columns; columns 3+ and lines beginning with "#" are ignored (comments). Columns one and two are key/value pairs respectively. See an example in "example/config.txt".

Config file parameters:

#### Input file paths
1. `rdcov_input_tsv` sample name:coverage filepath pairs file, tab-separated
1. `vcf_input_tsv` sample name:VCF filepath pairs file, tab-separated
1. `output_dir` Output files directory - existing files will be overwritten
1. `windows_bed` Path to window definitions BED file (see CLAMMS preprocessing steps)
1. `targets_bed` Path to target regions BED file
1. `chr_boundaries_bed` Path to pre-computed chromosomal boundary definitions BED file
1. `roh_afs_gz` Path to SNP allele frequencies file formatted for BCFtools/RoH --AF-file parameter (see bcftools manual)
1. `qc_metadata` Meta-data for samples, including QC metrics and reported gender
1. `karyo_sex` KaryoScan sex assignments file (optional)

#### Binary paths
1. `bcftools_dir` Path to bcftools executable w/ BCFtools/RoH compatability (e.g. /usr/local/bin/bcftools-1.x)
1. `bedtools_dir` Path to bedtools executable parent directory (e.g. /usr/local/bin/bedtools2/bin)
1. `datamash_dir` Path to datamash executable parent directory (e.g. /usr/local/bin)

#### Parameters
1. `ad_flag` Allele depth VCF flag (required) - e.g. 'AD'
1. `dp_flag` Total depth VCF flag (required) - e.g. 'DP'
1. `rd_gc_min` Minimum window GC content for read depth (RD) calculations - e.g. 0.45
1. `rd_gc_max` Maximum " - e.g. 0.55
1. `rd_mappability` Minimum window mappability for RD calculations - e.g. 0.95
1. `ab_mappability` Minimum window mappability for inclusion in allele balance (AB) analysis - e.g. 0.9
1. `homAB_threshold` Minimum lesser allele balance to consider as homozygous SNP (anything above this threshold is considered a putative het) - e.g. 0.02
`ab_mappability` Minimum window mappability for inclusion in allele balance (AB) analysis - e.g. 0.9
1. `bcftools_view_args` Command line arguments for bcftools view query/filter string - e.g. -i 'QD > 5 && DP >= 20' -M2 -v snps -f 'PASS'
1. `bcftools_roh_args` Command line arguments for bcftools roh - e.g. -G30 -I -a 6.6e-09 -H 5e-10
1. `chrY_dup_factor` chrY duplications are naively assigned to male samples having gamma_chrY >= `chrY_dup_factor` * median(gamma_chrY for all male samples). Alternatively, setting this to -1 enables chrY to be modeled like the autosomes and chrX.
1. `chromab_rd_covar` Sample read depth QC metric to use for chromosome-wide allele balance regression covariate - e.g. Picard metric PCT_TARGET_BASES_50X
1. `make_plots` TRUE/FALSE, flag for whether or not to create diagnostic plots (read depth and allele balance plots)


The `example/` directory contains example file formats for each file referenced in the config. These files are not meant to produce an example KaryoScan run, only to show what the expected file formats are. Note that the code uses the string "SampleID" in several places, be sure that any input files have headers named accordingly.

### Chr boudaries file
This BED file contains exon-capture boundaries for the chromosome. Specifically, it contains the start-stop coordinates of the centromeres buffered to the nearest exon boundaries, as well as the chormosomal end coordinate truncated to the last exon. If a chromosome is acrocentric (or otherwise contains no captured exons on a chromosomal arm), the centromere start-end coordinates are extended to encompass the full chromosomal arm that is absent of sequence coverage.

Note: these boundaries should be set for your specific sequencing capture regions. 'bedtools closest' is useful for this computation.

Columns:

1. Chromosome
1. Centromere start coordinate (buffered to nearest upstream exon)
1. Centromere end coordinate (buffered to nearest downstream exon)
1. Chromosomal end coordinate (buffered to nearest exon)
1. Expected Chromosome-wide Median Allele Balance (for visualization purposes only, e.g. 0.45)

## Running KaryoScan

The main driver program is `ks_driver.sh`.

    Usage: ./ks_driver.sh <config_file> [<run_mode>]
        run_mode options: 0 = run all steps (skip manual sex assignment, default)
                          1 = process read depths & allele balance files, stop for manual sex assignment
                          2 = proceed post-sex assignemnt (post run_mode=1 call)

The full KaryoScan analysis pipeline performs the following high-level steps:

1. Compute the `gamma_i` matrix based on the CLAMMS raw coverage files referenced in `rdcov_input_tsv`.
1. Scan the single sample VCF files referenced in `vcf_input_tsv`, outputting signals for chromosome-wide allele balance (AB), candidate local AB anomalies, and ROH segments. These output files are named `all_chrom_ab.txt`, `all_local_ab.txt`, and `all_roh.txt`, respectively.
1. [OPTIONAL] Generally, it is suggested that sex assignments be inspected and corrected here. This step can be skipped if a sex assignment is included in the QC metrics file and inspection is not desired. To inspect and correct sex assignments entails observing the `gamma_i` matrix distributions and determining a genetic sex assignment. A process to do this in R is provided in `assign_ks_sex.R`. Note that this is not an Rscript, but an interactive process - the comments should be read carefully. Alternative methods for inspecting and reassigning sex are fine, but should output a two-column file with SampleID:KaryoSex pairs.
1. Perform cohort-level regressions on `gamma_i` matrix and `all_chrom_ab.txt` files to identified anomalous values of each metric.
1. Aggregation of signals from raw ROH and local AB anomaly outputs as well as regressed `gamma_i` and chromosome-wide AB statistics.
1. Apply "decision tree" style logic to aggregated signals, scoring and ranking individual anamolies and generating a sample-level ranking and report for identified karyotypic anomalies.
1. Outputs of KaryoScan are the `karyoscan_anomalies.txt` summary file, and if `make_plots` was set, an output sub-directory `plots/` will contain the KaryoScan diagnostic plots.

### Output summary file karyoscan_anomalies.txt
This file has aggregated chromosomal anomaly calls at the individual sample and chromosome levels. Anomalies are scored into "tier" ratings, where Tier 1 anomalies are the most confident and Tier 3 are the least. Notably, tier rankings are assigned for individual anomaly signals (e.g., an ROH signal has a tier rating), and aggregated over a single chromosome for a sample to generate a chromosome-level tier rating (e.g. in a trisomy 21 sample, there may be a Tier 1 coverage anomaly and a Tier 1 AB anaomly, they are aggregated at the chromosome level to make a tier 1 chromosome 21 anomaly call).

In addition to ranking and scoring events, KaryoScan makes a call at each chromosome for each sample attempting to determine if the anomaly is:

* dosage-altering (gain, loss, or neutral)
* mosaic or not
* whole or partial chromosome scale

This final determination is made using the aggregate of anomalous signals identified for a chromosome. Thus each individual signal is flagged as a "Primary", "Secondary", or "Tertiary" signal with respect to its contribution to the determination. Using the same trisomy 21 example described above, both the coverage anomaly and chromosome-wide allele balance anomalies would be considered "Primary". However, if an ROH signal or local AB signal were also present, they may determined "Secondary" or "Tertiary" given that they do not drive the call determination, but may provide supporting (or unrelated) information about the predicted event.

The format of the final `karyoscan_anomalies.txt` output file is a tab-separated file with the following column definitions, where each row represents a single anomalous signal:

1. SampleID
1. Chromosome
1. Start position of anomaly ("1" if chromosome-wide event)
1. End position of anomaly (chromosome end coordinate if chromosome-wide event)
1. Event type (ROH, Coverage, ChromMedAB, or LocalMedAB)
1. Aggregate chromosomal fraction estimate from read coverage (1 = one chromosomal copy, aggregated over all events on this chromosome)
1. Aggregate chromosomal fraction estimate from allele balance (1 = one chromosomal copy, aggregated over all events on this chromosome))
1. KaryoScan anomaly dosage determination (gain, loss, or neutral)
1. KaryoScan anomaly scale determination (whole or partial chromosome)
1. KaryoScan anomaly mosaic determination (yes or no)
1. Raw event metric for estimating chromosomal fraction (context specific, may be an allele balance or coverage-based metric)
1. Event chromosomal fraction estimate from read depth or allele balance (context specific)
1. Sample sex
1. Tier rating of this specific event signal
1. Combined tier rating for the chromosome (aggregated over all raw event signals)

### Optimizing for your sequencing pipeline
The KaryoScan code here has been trained to fit our internal sequencing pipeline, with high-level parameters being configurable via the config file. However, it is likely that your sequencing pipeline may be better fit by adjusting internal parameters that are currently hard-coded. As described in the method paper (see reference above), many of the score thresholds (e.g. determination of tier ratings) require manual inspection of the distributions of KaryoScan metrics across a large subset of your samples to be properly tuned for your pipeline. This is not a process that should be automated; the user needs to inspect their data to determine what thresholds are applicable.

* To make modifications to the Local AB and ROH event thresholds, parameters are defined in `aggregate_signals.sh`
* To make modifications to the read coverage and chromosome-wide AB regression analyses, see `regress_rd_ab.R`
* To make modifications to the final event aggregation and scoring logic, see `call_anomalies.py`
