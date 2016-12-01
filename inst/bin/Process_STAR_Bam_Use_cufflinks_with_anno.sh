#!/bin/bash

#Usage: sh Process_STAR_Bam_Use_cufflinks_with_anno.sh  input_files_dir
#Example: sh ~/Script_bash/Process_STAR_Bam_Use_cufflinks_with_anno.sh /projects/scratch/bbc/Project2016/GOSJ/Alignment/

DIR="$1*__STAR"

results_dir="$1Output_cufflinks_with_anno/"

if [ -e $results_dir ]; then echo "already exists";
else
    mkdir $results_dir

fi


for file in $DIR

do

f=`echo "$file"`

echo "$f"

bam_file="$f/STAR_out.sorted.bam"

dir_name=$(dirname "$bam_file")
file_name=$(basename "$dir_name")

echo "$dir_name"
echo "$file_name"

sample_name=$(echo "$file_name" | awk -F"__" '{print $1}')

echo "$sample_name"

cat > ~/Script_bash/Run_"$sample_name"_cufflinks_anno.sh <<EOF

# LSBATCH: User input
#!/bin/bash
#BSUB -P bbc
#BSUB -J CufflinksAnno
#BSUB -o %J.CufflinksAnno.log
#BSUB -e %J.CufflinksAnno.err
#BSUB -W 48:00
#BSUB -n 16
#BSUB -q bigmem
#BSUB -R 'rusage[mem=36864] span[ptile=8]'
#BSUB -u aimin.yan@med.miami.edu
module load cufflinks/2.2.1

cufflinks $1$file_name/STAR_out.sorted.bam -o $1$file_name/cufflinks."$sample_name".with.anno -p 8 -G /projects/scratch/bbc/Homo_REF/Homo_sapiens.GRCh38.84.processed.sorted.gtf

mv $1$file_name/cufflinks."$sample_name".with.anno/transcripts.gtf $1$file_name/cufflinks."$sample_name".with.anno/"$sample_name".transcripts.gtf

echo $1$file_name/cufflinks."$sample_name".with.anno/"$sample_name".transcripts.gtf >> "$results_dir"assemblies_2.txt

EOF

bsub -K -P bbc < ~/Script_bash/Run_"$sample_name"_cufflinks_anno.sh& 
wait 

done

cat > ~/Script_bash/Run_cuffmerge_anno.sh <<EOF

# LSBATCH: User input
#!/bin/bash
#BSUB -P bbc
#BSUB -J CuffmergeAnno
#BSUB -o %J.CuffmergeAnno.log
#BSUB -e %J.CuffmergeAnno.err
#BSUB -W 48:00
#BSUB -n 16
#BSUB -q bigmem
#BSUB -R 'rusage[mem=36864] span[ptile=8]'
#BSUB -u aimin.yan@med.miami.edu
module load cufflinks/2.2.1

cuffmerge -o "$results_dir" -g /projects/scratch/bbc/Homo_REF/Homo_sapiens.GRCh38.84.processed.sorted.gtf -s /nethome/yxb173/Genome_Ref/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa -p 8 "$results_dir"assemblies_2.txt

EOF

cat > ~/Script_bash/Run_cuffdiff_merge_anno.sh <<EOF

# LSBATCH: User input
#!/bin/bash
#BSUB -P bbc
#BSUB -J CuffdiffAnno
#BSUB -o %J.CuffdiffAnno.log
#BSUB -e %J.CuffdiffAnno.err
#BSUB -W 48:00
#BSUB -n 16
#BSUB -q bigmem
#BSUB -R 'rusage[mem=36864] span[ptile=8]'
#BSUB -u aimin.yan@med.miami.edu
module load cufflinks/2.2.1

cuffdiff -o "$results_dir"diff_out -b /nethome/yxb173/Genome_Ref/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa -p 8 -L Mut,HC -u /nethome/axy148/merged_asm/merged.gtf "$1"SRR1660308__STAR/STAR_out.sorted.bam,"$1"SRR1660311__STAR/STAR_out.sorted.bam,"$1"SRR1660313__STAR/STAR_out.sorted.bam "$1"SRR1660320__STAR/STAR_out.sorted.bam,"$1"SRR1660321__STAR/STAR_out.sorted.bam,"$1"SRR1660322__STAR/STAR_out.sorted.bam,"$1"SRR1660323__STAR/STAR_out.sorted.bam,"$1"SRR1660324__STAR/STAR_out.sorted.bam

EOF

bsub -K -P bbc < ~/Script_bash/Run_cuffmerge_anno.sh&
wait

bsub -P bbc < ~/Script_bash/Run_cuffdiff_merge_anno.sh