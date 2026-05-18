## filter VCF based on chromosome

# first variable is the vcf
VCF_IN=$1
CHROMOSOME_LIST=$2
VCF_PREFIX=$3

cat ${CHROMOSOME_LIST} | while read i
do

echo "Filtering for chromosome ${i}..."
bcftools view ${VCF_IN} --regions Chr_${i} -o ${VCF_PREFIX}_chr${i}.vcf.gz -Oz
echo "Finished filtering for chromosome ${i}."

echo "Create PCA for chromosome ${i}..."
plink --vcf ${VCF_PREFIX}_chr${i}.vcf.gz --double-id --allow-extra-chr --make-bed --mind --pca --out ${VCF_PREFIX}_chr${i}
echo "Finished PCA for chromosome ${i}."

done