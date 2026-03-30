# Run TensorQTL from the BroadInstitute: https://github.com/broadinstitute/tensorqtl?tab=readme-ov-file

## load arguments
PLINK_DIR=$1
EXPRESSION_BED=$2
COVARIATE_FILE=$3
OUTPUT_PREFIX=$4

## map cis-QTLs for all variant-phenotype pairs
echo "Mapping cis-QTLS..."
python3 -m tensorqtl --covariates $COVARIATE_FILE --mode cis_nominal $PLINK_DIR $EXPRESSION_BED $OUTPUT_PREFIX 2> $4.cisqtl.error

## convert cis-QTL outputs into tsv format
echo "Converting cis-QTL outputs to CSV..."
for i in $4*.parquet; do
    parquet-tools csv $i >> $4.cis_qtl_pairs.txt
done

echo "Compressing cis-QTL output..."
gzip $4.cis_qtl_pairs.txt

## map trans-QTLS
echo "Mapping trans-QTLs..."
python3 -m tensorqtl --covariates $COVARIATE_FILE --mode trans --output_text $PLINK_DIR $EXPRESSION_BED $OUTPUT_PREFIX 2> $4.transqtl.error

## send message that running is complete
echo "Mapping Done."