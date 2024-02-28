#1/bin/bash -e

echo -e "donor\tnvars" > num_variants.txt

for donor in p0{07,08,09,13,14,16,20,21,26,35}; do 
        nvars=$(bcftools view -v snps ../../zenodo_repository/WGS/${donor}_filtered.vcf.gz | grep -v "^#" | wc -l)
        echo -e "${donor}\t${nvars}" >> num_variants.txt
done
