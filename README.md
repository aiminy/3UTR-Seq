# 3UTR-Seq
An R package for processing and analyzing 3UTR,DOGs and 5UTR data

# To Install

```{r eval=TRUE}
#In R console
library(devtools)

# This version of PathwaySplice can be installed if your R version is >= 3.4.0
install_github("aiminy/3UTR-Seq",ref = 'DoGs')

#In pegasus terminal 
R -e 'library(devtools);install_github("aiminy/3UTR-Seq",ref = "DoGs")'

```
# Usage

Rscript ~/R/lib64/R/library/ThreeUTR/bin/rscript/GetCountUseBash.r

Do you want to get count?

Yes

please define the input Bam file directory:

/nethome/zxg161/ZG-BBC-project/David_Ho/Alignment

please define the output file directory:

Test

please define the annotation file:

/nethome/axy148/3UTR_GFF2.bed

please define the upstream base of each transcript :

0

please define the downstream base of each transcript :

0

please define the sample that will not used for count :

No

Do you want to use bigmem to run your job:
Yes


Do you want to get count?

No

I already have counts, so I need to perform analysis.

Please define the following setting parameters: 

input.count.file.dir (ex: UTR/Counts)

input.file.pattern (ex: count.txt)

output.anlysis.dir (ex:UTR/Analysis)

out.file.pattern.interested(ex:DoGs_adjust_by_batch_interested_gene) 

out.file.pattern.positive.gene(ex:DoGs_adjust_by_batch_positive) 

out.file.pattern.negative.gene(ex:DoGs_adjust_by_batch_negative) 
 
out.file.pattern.all(ex:DoGs_adjust_by_batch_all) 

dir.name.gene.list(ex:/scratch/projects/exempt/tr/azhang/for_bioinfo_core/RNA_seq) 

pattern.4.gene.list(ex:final_list.csv) 

adjust_by_batch(ex:Yes) 
