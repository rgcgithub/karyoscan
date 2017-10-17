#!/bin/bash

set -e -o pipefail

if [ $# -lt 14 ]; then
    echo "Usage: ./scan_vcf_abstats.sh <sample:vcf input list TSV file> <bcftools_dir> <bcftools_view_args> <bcftools_roh_args> <allele depth VCF flag (e.g. AD)> <read depth VCF flag (e.g. DP)> <Min mappability> <homozygous AB threshold> <chromosome boundaries BED> <bcftools roh --AF-file> <windows BED> <targets BED> <make_plots_flag> <output stem>"
    exit
fi

fvcfmap="$1"
export bcftools_dir="$2"
export bcftools_view_args="$3"
export bcftools_roh_args="$4"
export ad_flag="$5"
export dp_flag="$6"
export ab_map="$7"
export homAB_threshold="$8"
export chr_boundaries="$9"
export afs_gz="${10}"
export windows="${11}"
export targets_bed="${12}"
export plot_flag="${13}"
export output_dir="${14}"

process_one() {
    in_args=(${1// / })

    sample="${in_args[0]}"
    vcf="${in_args[1]}"

    echo "Processing $sample $vcf"
    # Check for TBI file named ${vcf}.tbi
    if [ ! -e "${vcf}.tbi" ]; then
        command -v tabix >/dev/null 2>&1 || { echo >&2 "No index exists for ${vcf} and tabix is not in you PATH to create one. Aborting"; exit 1; }
        tabix -p vcf ${vcf}
    fi

    eval $bcftools_dir/bcftools view "${bcftools_view_args}" \
      -T <(awk -v abmap=${ab_map} '$NF == 3 && $(NF-1) >= abmap {print $1"\t"$2"\t"$3;}' $windows) $vcf -Ou \
    | $bcftools_dir/bcftools query -f "%CHROM\t%POS\t[%${ad_flag}{0}\t%${dp_flag}]\n%CHROM\t%POS\t[%${ad_flag}{1}\t%${dp_flag}]\n" \
    | awk '
$4 > 0 && ( ($3/$4) < 0.5 || (($3/$4) == 0.5 && (NR % 2) == 0) ){
    $3=$3/$4;
    print $1"\t"$2"\t"$3;
}' > ${output_dir}/${sample}.allele_balance.txt

    eval $bcftools_dir/bcftools roh "$bcftools_roh_args" --AF-file $afs_gz -R ${targets_bed} $vcf \
    | awk -v samp="$sample" '
BEGIN {
    start = -1;
    cur = -1;
    end = -1;
    cnt = 0;
    max_conf = 0;
    curchr = "";
}
$1 !~ /^#/ { 
    if (curchr == "")
        curchr = $1;
    cur = $2*$3;
    chr = $1;
    conf = $4;
    if (curchr != chr || cur == 0) {
        if (start >= 0) {
            print curchr"\t"start"\t"end"\t"end-start"\t"cnt"\t"max_conf"\t"samp;
            start = -1;
            end = -1;
            cnt = 0;
            max_conf = 0;
        }
    }
    if (cur > 0) {
        if (start < 0) {
            start = cur;
            cnt++;
            end = cur;
            max_conf = conf;
        }
        else {
            cnt++;
            end = cur;
            if (conf > max_conf)
                max_conf = conf;
        }
    }
    curchr = chr;
}
END {
    if (start >= 0)
        print curchr"\t"start"\t"end"\t"end-start"\t"cnt"\t"max_conf"\t"samp;
}' > ${output_dir}/${sample}.roh.txt

    Rscript analyze_plot_allele_balance.R ${output_dir}/${sample}.allele_balance.txt ${output_dir}/${sample}.roh.txt $chr_boundaries $homAB_threshold $plot_flag
}

export -f process_one

while read sample vcf; do
    echo "$sample $vcf";
done < $fvcfmap \
| xargs -P $(nproc) -I{} -n 1 bash -c 'process_one "$@"' _ {}

cat ${output_dir}/*.roh.txt > ${output_dir}/all_roh.txt
cat ${output_dir}/*.eventstats.txt > ${output_dir}/all_local_ab.txt
cat ${output_dir}/*.chrmedians.txt > ${output_dir}/all_chrom_ab.txt

rm ${output_dir}/*.roh.txt
rm ${output_dir}/*.eventstats.txt
rm ${output_dir}/*.chrmedians.txt
rm ${output_dir}/*.allele_balance.txt

mv ${output_dir}/*.allele_balance.png ${output_dir}/plots/
