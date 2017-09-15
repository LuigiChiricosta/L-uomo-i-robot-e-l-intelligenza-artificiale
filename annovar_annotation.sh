#!/bin/bash
vcfin=$1
vcfout=$2
export ANNODIR=/export/um99/bin/annovar
export VCFTOOLS=/export/um99/bin/vcftools_0.1.13/bin
$ANNODIR/table_annovar.pl $vcfin $ANNODIR/humandb/ -buildver hg19 -out $vcfout -remove -protocol refGene,gnomad_genome -operation g,f -nastring . -vcfinput
rm $vcfout.avinput
rm $vcfout.hg19_multianno.txt
$VCFTOOLS/vcf-sort -c  $vcfout.hg19_multianno.vcf |gzip -c > $vcfout.hg19_multianno.vcf.gz
