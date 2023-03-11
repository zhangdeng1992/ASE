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
path=[]
for file in files:
    path.append(file)

os.chdir('LPF1_5MM_huaguan')

def star(deviceid1,deviceid2):
    os.system('%s --twopassMode Basic --runThreadN 2 --genomeDir /gss1/home/dingbq01/zhangdeng/Allele_specific_expression/Cleandata/data/%s --alignIntronMin 20 --alignIntronMax 50000 --outSAMmapqUnique 60 --outSAMtype BAM SortedByCoordinate --sjdbOverhang 149 --outFileNamePrefix %s_%s. --readFilesCommand gunzip -c --readFilesIn /gss1/home/dingbq01/zhangdeng/Allele_specific_expression/Cleandata/data/%s_R1.fq.gz /gss1/home/dingbq01/zhangdeng/Allele_specific_expression/Cleandata/data/%s_R2.fq.gz' %(STAR,deviceid1, deviceid2,deviceid1, deviceid2,deviceid2))


for file in path:
    for b in ['lg','mp']:
        t = Thread(target=star, args=(b,file[0:-9]))
        t.start()