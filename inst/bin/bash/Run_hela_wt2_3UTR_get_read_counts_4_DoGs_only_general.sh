
#LSBATCH: User input
#!/bin/bash
#BSUB -P bbc
#BSUB -J DogFT
#BSUB -o %J.DogFT.log
#BSUB -e %J.DogFT.err
#BSUB -W 72:00
#BSUB -n 8
#BSUB -q general
#BSUB -u aimin.yan@med.miami.edu
module load bedtools/2.17.0

#For each gene, get count of the reads that fall in the 4500bp region of downstream of each gene

#+-
awk -F '\t' '$6=="+"&&$12=="-"' "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/""hela_wt2".DoGs.only.strand.based.hg19.total2.bed | awk '$8<$2&&$9>=$2'| awk '{print $4}' | sort | uniq -c | sort -nr > "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/Counts/""hela_wt2".plus.gene.minus.read.below.DoGs.count.2.txt

awk -F '\t' '$6=="+"&&$12=="-"' "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/""hela_wt2".DoGs.only.strand.based.hg19.total2.bed | awk '$8>=$2&&$9<=$3'| awk '{print $4}' | sort | uniq -c | sort -nr > "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/Counts/""hela_wt2".plus.gene.minus.read.DoGs.count.2.txt

awk -F '\t' '$6=="+"&&$12=="-"' "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/""hela_wt2".DoGs.only.strand.based.hg19.total2.bed | awk '$8<=$3&&$9>$3' | awk '{print $4}' | sort | uniq -c | sort -nr > "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/Counts/""hela_wt2".plus.gene.minus.read.over.DoGs.count.2.txt


#++
awk -F '\t' '$6=="+"&&$12=="+"' "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/""hela_wt2".DoGs.only.strand.based.hg19.total2.bed | awk '$8<$2&&$9>=2'| awk '{print $4}' | sort | uniq -c | sort -nr > "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/Counts/""hela_wt2".plus.gene.plus.read.below.DoGs.count.2.txt

awk -F '\t' '$6=="+"&&$12=="+"' "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/""hela_wt2".DoGs.only.strand.based.hg19.total2.bed | awk '$8>=$2&&$9<=$3'| awk '{print $4}' | sort | uniq -c | sort -nr > "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/Counts/""hela_wt2".plus.gene.plus.read.DoGs.count.2.txt

awk -F '\t' '$6=="+"&&$12=="+"' "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/""hela_wt2".DoGs.only.strand.based.hg19.total2.bed | awk '$8<=$3&&$9>$3' | awk '{print $4}' | sort | uniq -c | sort -nr > "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/Counts/""hela_wt2".plus.gene.plus.read.over.DoGs.count.2.txt

#-+
awk -F '\t' '$6=="-"&&$12=="+"' "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/""hela_wt2".DoGs.only.strand.based.hg19.total2.bed | awk '$8<$2&&$9>=$2' | awk '{print $4}' | sort | uniq -c | sort -nr > "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/Counts/""hela_wt2".minus.gene.plus.read.over.DoGs.count.2.txt

awk -F '\t' '$6=="-"&&$12=="+"' "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/""hela_wt2".DoGs.only.strand.based.hg19.total2.bed | awk '$8>=$2&&$9<=$3'| awk '{print $4}' | sort | uniq -c | sort -nr > "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/Counts/""hela_wt2".minus.gene.plus.read.DoGs.count.2.txt

awk -F '\t' '$6=="-"&&$12=="+"' "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/""hela_wt2".DoGs.only.strand.based.hg19.total2.bed | awk '$8<=$3&&$9>$3' | awk '{print $4}' | sort | uniq -c | sort -nr > "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/Counts/""hela_wt2".minus.gene.plus.read.below.DoGs.count.2.txt


#--
awk -F '\t' '$6=="-"&&$12=="-"' "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/""hela_wt2".DoGs.only.strand.based.hg19.total2.bed | awk '$8<$2&&$9>=$2'| awk '{print $4}' | sort | uniq -c | sort -nr > "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/Counts/""hela_wt2".minus.gene.minus.read.over.DoGs.count.2.txt

awk -F '\t' '$6=="-"&&$12=="-"' "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/""hela_wt2".DoGs.only.strand.based.hg19.total2.bed | awk '$8>=$2&&$9<=$3'| awk '{print $4}' | sort | uniq -c | sort -nr > "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/Counts/""hela_wt2".minus.gene.minus.read.DoGs.count.2.txt

awk -F '\t' '$6=="-"&&$12=="-"' "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/""hela_wt2".DoGs.only.strand.based.hg19.total2.bed | awk '$8<=$3&&$9>$3' | awk '{print $4}' | sort | uniq -c | sort -nr > "/scratch/projects/bbc/Ramin_azhang/for_bioinfo_core/Results4DoGsOnly/Counts/""hela_wt2".minus.gene.minus.read.below.DoGs.count.2.txt

