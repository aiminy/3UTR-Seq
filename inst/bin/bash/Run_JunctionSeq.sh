java -Xmx4000M -jar /home/aiminyan/QoRTs/QoRTsFullExampleData/QoRTsRelease/QoRTs.jar QC --singleEnded --noGzipOutput --maxReadLength 1000 --keepMultiMapped /media/H_driver/2015/Nimer_Cheng/164.alignments.bam /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.processed.sorted.gtf /media/H_driver/2015/Nimer_Cheng/164_junction_seq_new_gtf_7


java -Xmx4000M -jar /home/aiminyan/QoRTs/QoRTsFullExampleData/QoRTsRelease/QoRTs.jar makeFlatGff /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.processed.sorted.gtf /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.JunctionSeq.flat.gff

java -Xmx4000M -jar /home/aiminyan/QoRTs/QoRTsFullExampleData/QoRTsRelease/QoRTs.jar makeFlatGff --stranded --DEXSeqFmt /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.processed.sorted.gtf /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.JunctionSeq.flat.1.gff

java -Xmx4000M -jar /home/aiminyan/QoRTs/QoRTsFullExampleData/QoRTsRelease/QoRTs.jar makeFlatGff --stranded /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.processed.sorted.gtf /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.JunctionSeq.flat.2.gff

java -Xmx4000M -jar /home/aiminyan/QoRTs/QoRTsFullExampleData/QoRTsRelease/QoRTs.jar makeFlatGff --DEXSeqFmt /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.processed.sorted.gtf /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.JunctionSeq.flat.3.gff

#Use strand information and specify to the first strand
java -Xmx4000M -jar /home/aiminyan/QoRTs/QoRTsFullExampleData/QoRTsRelease/QoRTs.jar QC --singleEnded ----stranded --noGzipOutput --maxReadLength 1000 --keepMultiMapped /media/H_driver/2015/Nimer_Cheng/164.alignments.bam /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.JunctionSeq.flat.2.gff /media/H_driver/2015/Nimer_Cheng/164_junction_seq_new_gtf_stranded

#Use starnd inforamtion and specify to the second strand
java -Xmx4000M -jar /home/aiminyan/QoRTs/QoRTsFullExampleData/QoRTsRelease/QoRTs.jar QC --singleEnded --stranded --stranded_fr_secondstrand --noGzipOutput --maxReadLength 1000 --keepMultiMapped /media/H_driver/2015/Nimer_Cheng/164.alignments.bam /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.JunctionSeq.flat.2.gff /media/H_driver/2015/Nimer_Cheng/164_junction_seq_new_gtf_second_stranded
