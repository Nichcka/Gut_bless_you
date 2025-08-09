MAPPING_DIR=$1
SQ_LOCI_META=$2
SRARUNTABLE=$3

./all_tpm.sh $MAPPING_DIR all_tpm.tsv

(echo -e "assemble\tlocus\ttpm"; cat all_tpm.tsv) > 1.tsv

awk 'BEGIN{FS=OFS="\t"} 
    NR==1 {print $1, $2, $3, "log_tpm"; next} 
    {printf "%s\t%s\t%.6f\t%.6f\n", $1, $2, $3, log($3+1)}' 1.tsv > 2.tsv

python make_table.py $SQ_LOCI_META

python add.py $SRARUNTABLE

rm 1.tsv 2.tsv 3.tsv
