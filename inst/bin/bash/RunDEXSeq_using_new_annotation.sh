#python /home/aiminyan/R/x86_64-pc-linux-gnu-library/3.2/DEXSeq/python_scripts/dexseq_count.py -p no -s no -a 10 /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.DEXSeq.gtf /media/H_driver/2015/Nimer_Cheng/164.alignments.sam /media/H_driver/2015/Nimer_Cheng/164.alignments.count

#convert bam files to sam files
#samtools view -h -o out.sam in.bam

#samtools sort -n /media/H_driver/2015/Nimer_Cheng/162.alignments.sam /media/H_driver/2015/Nimer_Cheng/162.alignments.sorted.sam
#samtools sort -n /media/H_driver/2015/Nimer_Cheng/163.alignments.sam /media/H_driver/2015/Nimer_Cheng/163.alignments.sorted.sam
#samtools sort -n /media/H_driver/2015/Nimer_Cheng/164.alignments.sam /media/H_driver/2015/Nimer_Cheng/164.alignments.sorted.sam
#samtools sort -n /media/H_driver/2015/Nimer_Cheng/78.alignments.sam /media/H_driver/2015/Nimer_Cheng/78.alignments.sorted.sam
#samtools sort -n /media/H_driver/2015/Nimer_Cheng/79.alignments.sam /media/H_driver/2015/Nimer_Cheng/79.alignments.sorted.sam
#samtools sort -n /media/H_driver/2015/Nimer_Cheng/80.alignments.sam /media/H_driver/2015/Nimer_Cheng/80.alignments.sorted.sam

#python /home/aiminyan/R/x86_64-pc-linux-gnu-library/3.2/DEXSeq/python_scripts/dexseq_count.py -p no -s no -a 10 /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.DEXSeq.12cs.gtf /media/H_driver/2015/Nimer_Cheng/78.alignments.sam /media/H_driver/2015/Nimer_Cheng/78.alignments.12cs.count

python /home/aiminyan/R/x86_64-pc-linux-gnu-library/3.2/DEXSeq/python_scripts/dexseq_count.py -p no -s no -a 10 /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.DEXSeq.cleaned.sorted.gtf /media/H_driver/2015/Nimer_Cheng/162.alignments.sam /media/H_driver/2015/Nimer_Cheng/162.alignments.new.anno.count

python /home/aiminyan/R/x86_64-pc-linux-gnu-library/3.2/DEXSeq/python_scripts/dexseq_count.py -p no -s no -a 10 /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.DEXSeq.cleaned.sorted.gtf /media/H_driver/2015/Nimer_Cheng/163.alignments.sam /media/H_driver/2015/Nimer_Cheng/163.alignments.new.anno.count

python /home/aiminyan/R/x86_64-pc-linux-gnu-library/3.2/DEXSeq/python_scripts/dexseq_count.py -p no -s no -a 10 /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.DEXSeq.cleaned.sorted.gtf /media/H_driver/2015/Nimer_Cheng/164.alignments.sam /media/H_driver/2015/Nimer_Cheng/164.alignments.new.anno.count

python /home/aiminyan/R/x86_64-pc-linux-gnu-library/3.2/DEXSeq/python_scripts/dexseq_count.py -p no -s no -a 10 /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.DEXSeq.cleaned.sorted.gtf /media/H_driver/2015/Nimer_Cheng/78.alignments.sam /media/H_driver/2015/Nimer_Cheng/78.alignments.new.anno.count

python /home/aiminyan/R/x86_64-pc-linux-gnu-library/3.2/DEXSeq/python_scripts/dexseq_count.py -p no -s no -a 10 /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.DEXSeq.cleaned.sorted.gtf /media/H_driver/2015/Nimer_Cheng/79.alignments.sam /media/H_driver/2015/Nimer_Cheng/79.alignments.new.anno.count

python /home/aiminyan/R/x86_64-pc-linux-gnu-library/3.2/DEXSeq/python_scripts/dexseq_count.py -p no -s no -a 10 /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.DEXSeq.cleaned.sorted.gtf /media/H_driver/2015/Nimer_Cheng/80.alignments.sam /media/H_driver/2015/Nimer_Cheng/80.alignments.new.anno.count
