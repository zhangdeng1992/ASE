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

#获取短变异
def gvcf_lf(a,c):
    os.system('%s --java-options "-Xmx100g -XX:ParallelGCThreads=2" HaplotypeCaller \
    -I %s \
    -R /gss1/home/dingbq01/reference/LF10g_v2.0/LF10g_v2.0.fa \
    -ERC GVCF \
    -O %s.g.vcf.gz' %(gatk, a,c))

def gvcf_mp(a,c):
    os.system('%s --java-options "-Xmx100g -XX:ParallelGCThreads=2" HaplotypeCaller \
    -I %s \
    -R /gss1/home/dingbq01/reference/Mparishii_genome/Mparg_v2.0.fa \
    -ERC GVCF \
    -O %s.g.vcf.gz' %(gatk, a,c))

paths = glob('*_lg.add.bam')
lg=[]
for path in paths:
    lg.append(os.path.basename(path))
paths = glob('*_mp.add.bam')
mp=[]
for path in paths:
    mp.append(os.path.basename(path))
lg.sort()
mp.sort()   
for i in range(0,4):
    t = Thread(target=gvcf_lf, args=(lg[i],lg[i][0:n]))
    t.start()

for i in range(0,4):
    t = Thread(target=gvcf_mp, args=(mp[i],mp[i][0:n]))
    t.start()
