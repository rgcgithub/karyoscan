#!/bin/bash

usage_str="<coverage_anomalies> <chromAB_anomalies> <localAB_events> <roh_events> <chromomosome_boundardies_bed> <sex_qc_file> <output_file> [<bedtools_dir> <datamash_dir>]"

if [ $# -lt 7 ]; then
    echo "Usage: $0 $usage_str"
    exit 1
fi

coverage="$1"
chromab="$2"
localab="$3"
roh="$4"
centromeres="$5"
sex_qc="$6"
fout="$7"

bedtools_dir=""
datamash_dir=""

if [ $# -gt 7 ]; then
    bedtools_dir="${8}/"
    datamash_dir="${9}/"
fi

locab_area=0.02     # min event area cutoff for local running median anomalies
locab_slope=-3000   # linear cutoff slope for local median anomalies
locab_y0_tier1=230  # linear cutoff intercept, tier 1
locab_y0_tier2=170  # linear cutoff intercept, tier 2
locab_y0_tier3=110  # linear cutoff intercept, tier 3
roh_minbp=5000000   # min size (b.p.) of an ROH event to consider
roh_minsnp=40	    # min # of SNPs in an ROH event to consider
roh_agg=20000000    # min total roh events on a chromosome (sum of events passing roh_event)

echo "Parsing $coverage file..."
awk 'BEGIN{OFS="\t";} NR>1{
sampleid=$1; chr=substr($2,4); fraction=$3; p=$4; q=$5; sex=$7; tier=$8;
print sampleid, chr, 0, 0, "Coverage", tier, fraction, p, q, sex;
}' $coverage \
> ${fout}.tmp

echo "Parsing $chromab file..."
awk 'BEGIN{OFS="\t";} {
if (NR == FNR) { #Parsing centromere coordinates file
    chr=$1; chrend=$4;
    chrends[chr] = chrend;
} else { #Parsing the ChromMedAB anomalies file
    sampleid=$1; chr=$2; ab=$3; sex=$5; fit=$6; residual=$7; p=$8; q=$9; tier=$10;
    if (tier <= 3) {
        print sampleid, chr, 1, chrends[chr], "ChromMedAB", tier, ab, p, q, fit;
    }
} }' $centromeres $chromab \
>> ${fout}.tmp

echo "Parsing $localab file..."
# Parse local_median_events file (Anomalies in the running median allele balance (sub-chromosome scale))
awk '{print $7"\t"$0;}' $localab \
| sort -k1,1 \
| join -t$'\t' - <(sort -k1,1 $sex_qc) \
| cut -f 2- \
| awk -v min_area="${locab_area}" -v slope="${locab_slope}" \
      -v y0_1="${locab_y0_tier1}" -v y0_2="${locab_y0_tier2}" -v y0_3="${locab_y0_tier3}" \
'BEGIN {OFS="\t";}
$1 !~ /^Y/ && ($1 !~ /^X/ || $NF ~ /^f/ || $NF ~ /^F/ || $NF == 2) && $5 > min_area {
    chr=$1; start=$2; end=$3; nsnps=$4; area=$5; minab=$6; sampleid=$7;
    tier = 0;
    test_nsnps = nsnps - (area * slope); # Given area, minimum SNPs required to pass tiers
    if (test_nsnps >= y0_1) {
        tier = 1;
    } else if (test_nsnps >= y0_2) {
        tier = 2;
    } else if (test_nsnps >= y0_3) {
        tier = 3;
    }
    
    if (tier > 0) {
       print sampleid, chr, start, end, "LocalMedAB", tier, minab, area, nsnps, ".";
    }
}' \
>> ${fout}.tmp

echo "Parsing $roh file (this will take the longest)..."
# Parse roh_events file, (Large runs-of-homozygosity anomalies)
awk '{print $7"\t"$0;}' $roh \
| sort -k1,1 \
| join -t$'\t' - <(sort -k1,1 $sex_qc) \
| cut -f 2- \
| awk -v min_bp="${roh_minbp}" -v min_snp="${roh_minsnp}" 'BEGIN{OFS="\t";}
$1 !~ /^Y/ && ($1 != "X" || $NF ~ /^F/ || $NF ~ /^f/ || $NF == 2) && $4 > min_bp && $5 > min_snp {
    $4=$4"\t"$1":"$2":"$3":"$5":"$6;
    print;
}' \
| sort -k1,1 -k2,2n \
| ${bedtools_dir}bedtools subtract -a - -b $centromeres \
| awk -v min_bp="${roh_minbp}" 'BEGIN {OFS="\t";} 
($3-$2) > min_bp {
    $4=$3-$2;
    print $8"\t"$0;
}' \
| sort -k1,1 \
| ${datamash_dir}datamash -g 1 sum 5 unique 6 \
| awk -v min_agg_bp="${roh_agg}" 'BEGIN{OFS="\t";} {
tier = 2;
if ($2 > min_agg_bp) {
    tier = 1;
}
sample=$1; tot_rohbp=$2;
n=split($3,events,",");
for (i=1; i<=n; i++) {
    split(events[i],tokens,":");
    print sample, tokens[1], tokens[2], tokens[3], "ROH", tier, 0, tokens[4], tokens[5], tot_rohbp;
} }' \
>> ${fout}.tmp

echo "Parsed all files, merging output"
sort -k1,1 -k2,2 -k3,3n ${fout}.tmp > ${fout}
rm "${fout}.tmp"
echo "Done"
