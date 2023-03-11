#原始数据过滤
import os 
from glob import glob
from threading import Thread
lg_gtf='/gss1/home/dingbq01/reference/LF10g_v2.0/LF10g_v2.0.gtf'
mp_gtf='/gss1/home/dingbq01/reference/Mparishii_genome/Mparg_v2.0.gtf'
trimmomatic='/gss1/home/dingbq01/software/Trimmomatic/bin/trimmomatic'
gatk='/gss1/App_dir/anaconda3/envs/gatk/bin/gatk'
lg_fa='/gss1/home/dingbq01/reference/LF10g_v2.0/LF10g_v2.0.fa'
mp_fa='/gss1/home/dingbq01/reference/Mparishii_genome/Mparg_v2.0.fa'
samtool='/gss1/App_dir/anaconda3/envs/samtools-1.11/bin/samtools'
fastqc='/gss1/App_dir/anaconda3/envs/fastqc/bin/fastqc'
hisat='/gss1/App_dir/anaconda3/envs/hisat2/bin'
STAR='/gss1/home/dingbq01/anaconda3/envs/STAR/bin/STAR'
picard='/gss1/App_dir/anaconda3/envs/picard/bin/picard'
bedtools='/gss1/App_dir/anaconda3/envs/bedtools_v2.30.0/bin/bedtools'

files=glob('LPF1-5mm-huaguan*_R1.fq.gz')
n=len(files[0][0:-9])+3

os.chdir('LPF1_5MM_huaguan')
#bam文件去重
paths = glob('*Aligned.sortedByCoord.out.bam')
for path in paths:
    file = os.path.basename(path)
    os.system('%s MarkDuplicates -I %s -O %s.marked.bam -M %s.dups.txt' %(gatk, file, file[0:n], file[0:n]))
    os.system('%s MarkDuplicates REMOVE_DUPLICATES=true  I= %s.marked.bam O= %s.dup.bam M= %s.out.txt' %(picard,file[0:n],file[0:n],file[0:n]))

os.system('rm -rf *.marked.bam')

#修改BAM文件的Read Group
paths = glob('*.dup.bam')
lg=[]
for path in paths:
    lg.append(os.path.basename(path))
for i in range(0,8):
    os.system('%s AddOrReplaceReadGroups -I %s -O %s.add.bam -LB %s -PL ILLUMINA -PU %s -SM %s' %(gatk, lg[i], lg[i][0:n], lg[i][0:n-3], lg[i][0:n-3], lg[i][0:n-3]))


#建立索引
paths = glob('*.add.bam')
lg=[]
for path in paths:
    lg.append(os.path.basename(path)) 
for i in range(0,8):
    os.system('%s index %s' %(samtool, lg[i]))