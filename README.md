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
