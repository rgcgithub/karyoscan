#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: ./ks_driver.sh <config_file> [<run_mode>]"
    echo -e "\t run_mode options: 0 = run all steps (skip manual sex assignment)"
    echo -e "\t                   1 = process read depths & allele balance files, stop for manual sex assignment"
    echo -e "\t                   2 = proceed post-sex assignemnt (post run_mode=1 call)"
    exit 0
fi

fconfig="$1"
run_mode="0"
if [ $# -gt 1 ]; then
    run_mode="$2"
fi

if [[ "$run_mode" != "0" && "$run_mode" != "1" && "$run_mode" != "2" && "$run_mode" != "3" && "$run_mode" != "4" ]]; then
    echo "Invalid run_mode option: $run_mode"
    echo "Usage: ./ks_driver.sh <config_file> [<run_mode>]"
    echo -e "\t run_mode options: 0 = run all steps (skip manual sex assignment)"
    echo -e "\t                   1 = process read depths & allele balance files, stop for manual sex assignment"
    echo -e "\t                   2 = proceed post-sex assignemnt (post run_mode=1 call)"
    exit 0
fi

KS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Read config file params
echo "Reading configuration file...."
declare -A cfgvars
while read key val; do
    cfgvars["$key"]="$val"
done < <(awk -F'\t' '$1 !~ /^$/ && $1 !~ /^#/{print $1"\t"$2;}' $fconfig)

mkdir -p ${cfgvars[output_dir]} ${cfgvars[output_dir]}/plots

if [[ "$run_mode" == "0" || "$run_mode" == "1" ]]; then
    # Analyze read-coverage data:
    echo "Calculating r_i, gamma_i (read coverage chromosomal summary) values for all samples in ${cfgvars[rdcov_input_tsv]} file..."
    ${KS_DIR}/calculate_gammai.py ${cfgvars[windows_bed]} \
                          ${cfgvars[rd_mappability]} \
                          ${cfgvars[rd_gc_min]} \
                          ${cfgvars[rd_gc_max]} \
                          ${cfgvars[rdcov_input_tsv]} \
                          "${cfgvars[output_dir]}/gammai_matrix.txt"

    # Analyze allele-balance data from VCFs
    echo "Parsing VCF files to identify het allele balance metrics and ROH segments..."
    ${KS_DIR}/scan_vcf_abstats.sh ${cfgvars[vcf_input_tsv]} \
                          ${cfgvars[bcftools_dir]} \
                          "${cfgvars[bcftools_view_args]}" \
                          "${cfgvars[bcftools_roh_args]}" \
                          ${cfgvars[ad_flag]} \
                          ${cfgvars[dp_flag]} \
                          ${cfgvars[ab_mappability]} \
                          ${cfgvars[homAB_threshold]} \
                          ${cfgvars[chr_boundaries_bed]} \
                          ${cfgvars[roh_afs_gz]} \
                          ${cfgvars[windows_bed]}\
                          ${cfgvars[targets_bed]} \
                          ${cfgvars[make_plots]} \
                          ${cfgvars[output_dir]}
    
    echo "Completed stage 1 (gammai matrix calculation and VCF parsing)..."
fi

if [[ "$run_mode" == "1" ]]; then
    echo "Follow assign_ks_sex.R code tutorial to assign sex and chrY duplication threshold. Then proceed with run_mode=2 (./ks_driver $fconfig 2)."
    echo "Done"
    exit 0
fi

if [[ "$run_mode" == "0" || "$run_mode" == "2" || "$run_mode" == "3" ]]; then
    echo "Performing cohort-wide regression modelling on read-coverage (gamma_i) and chromosomal het allele balance (ChromHetAB)..."
    Rscript ${KS_DIR}/regress_rd_ab.R "${cfgvars[output_dir]}/gammai_matrix.txt" \
                            "${cfgvars[output_dir]}/all_chrom_ab.txt" \
                            ${cfgvars[qc_metadata]} \
                            ${cfgvars[windows_bed]} \
                            ${cfgvars[chrY_dup_factor]} \
                            ${cfgvars[chromab_rd_covar]} \
                            ${cfgvars[make_plots]} \
                            ${cfgvars[output_dir]} \
                            ${cfgvars[karyo_sex]}

    echo "Completed stage 2 (gamma_i and ChromHetAB regression)..."
fi

if [[ "$run_mode" == "0" || "$run_mode" == "2" || "$run_mode" == "4" ]]; then
    fsex_arg="${cfgvars[qc_metadata]}"
    if [ ${cfgvars[karyo_sex]} ]; then
        fsex_arg="${cfgvars[karyo_sex]}"
    fi

    echo "Aggregating signals from all metrics..."
    ${KS_DIR}/aggregate_signals.sh "${cfgvars[output_dir]}/coverage_anomalies.txt" \
                           "${cfgvars[output_dir]}/chrom_ab_anomalies.txt" \
                           "${cfgvars[output_dir]}/all_local_ab.txt" \
                           "${cfgvars[output_dir]}/all_roh.txt" \
                           ${cfgvars[chr_boundaries_bed]} \
                           $fsex_arg \
                           "${cfgvars[output_dir]}/all_signals.txt" \
                           ${cfgvars[bedtools_dir]} \
                           ${cfgvars[datamash_dir]}

    echo "Scoring anomalies..."
    ${KS_DIR}/call_anomalies.py "${cfgvars[output_dir]}/all_signals.txt" \
                        $fsex_arg \
                        ${cfgvars[chr_boundaries_bed]} \
    > ${cfgvars[output_dir]}/karyoscan_anomalies.txt

    echo "Completed stage 3 (signal aggregation and event scoring)..."
fi

echo "Done, all output placed in ${cfgvars[output_dir]}"
