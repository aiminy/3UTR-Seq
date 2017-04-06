$annovar/annotate_variation.pl -buildver mm10 -downdb -webfrom annovar refGene mm10db/
$annovar/annotate_variation.pl --buildver mm10 --downdb seq mm10db/mm10_seq
$annovar/retrieve_seq_from_fasta.pl mm10db/mm10_refGene.txt -seqdir mm10db/mm10_seq -format refGene -outfile mm10db/mm10_refGeneMrna.fa
$annovar/annotate_variation.pl -buildver mm10 -downdb cytoBand mm10db/
$annovar/annotate_variation.pl -buildver mm10 -downdb genomicSuperDups humandb/ 
$annovar/annotate_variation.pl -buildver mm10 -downdb -webfrom annovar esp6500siv2_all mm10db/
$annovar/annotate_variation.pl -buildver mm10 -downdb -webfrom annovar 1000g2014oct mm10db/
$annovar/annotate_variation.pl -buildver mm10 -downdb -webfrom annovar snp138 mm10db/ 
$annovar/annotate_variation.pl -buildver mm10 -downdb -webfrom annovar ljb26_all mm10db/

#$annovar/table_annovar.pl example/ex1.avinput mm10db/ -buildver mm10 -out myanno -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring . -csvout

$annovar/table_annovar.pl myAnalysis/results/passed.somatic.snvs.vcf mm10db/ -buildver mm10 -out myanno -remove -protocol refGene,cytoBand,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring . -vcfinput

$annovar/table_annovar.pl myAnalysis/results/passed.somatic.snvs.vcf mm10db/ -buildver mm10 -out myanno -remove -protocol refGene,cytoBand  -operation g,r -nastring . -vcfinput