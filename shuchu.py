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

'''
#gVCF转VCF
paths = glob('*_lg.g.vcf.gz')
lg=[]
for path in paths:
    lg.append(os.path.basename(path))
paths = glob('*_mp.g.vcf.gz')
mp=[]
for path in paths:
    mp.append(os.path.basename(path))
lg.sort()
mp.sort()   
for i in range(0,4):
    os.system('%s --java-options "-Xmx100g -XX:ParallelGCThreads=4" GenotypeGVCFs -R %s -V %s -O %s.vcf.gz'%(gatk, lg_fa, lg[i], lg[i][0:n]))
    os.system('%s --java-options "-Xmx100g -XX:ParallelGCThreads=4" GenotypeGVCFs -R %s -V %s -O %s.vcf.gz' %(gatk, mp_fa, mp[i], mp[i][0:n]))
#vcf文件提取SNP
paths = glob('*_lg.vcf.gz')
lg=[]
for path in paths:
    lg.append(os.path.basename(path))
paths = glob('*_mp.vcf.gz')
mp=[]
for path in paths:
    mp.append(os.path.basename(path))
lg.sort()
mp.sort()   
for i in range(0,4):
    os.system('%s SelectVariants -R %s -V %s --select-type-to-include SNP -O %s_snp.vcf'%(gatk, lg_fa, lg[i], lg[i][0:n]))
    os.system('%s SelectVariants -R %s -V %s --select-type-to-include SNP -O %s_snp.vcf' %(gatk, mp_fa, mp[i], mp[i][0:n]))



#VCF文件质控
paths = glob('*_lg_snp.vcf')
lg=[]
for path in paths:
    lg.append(os.path.basename(path))
paths = glob('*_mp_snp.vcf')
mp=[]
for path in paths:
    mp.append(os.path.basename(path))
lg.sort()
mp.sort()   
for i in range(0,4):
    os.system('%s VariantFiltration -R %s -V %s -O %s_Filter_SNP.vcf \
   --filter-name "QUAL30" --filter-expression "QUAL < 30.0" \
   --filter-name "FS60" --filter-expression "FS > 60.0" \
   --filter-name "QD2" --filter-expression "QD < 2.0" \
   --filter-name "SOR3" --filter-expression "SOR > 3.0" \
   --filter-name "MQ40" --filter-expression "MQ < 40.0" \
   --filter-name "ReadPosRankSum-8.0" --filter-expression "ReadPosRankSum < -8.0" \
   --filter-name "MQRankSum-12.5" --filter-expression "MQRankSum < -12.5"'%(gatk, lg_fa, lg[i], lg[i][0:n]))
    os.system('%s VariantFiltration -R %s -V %s -O %s_Filter_SNP.vcf \
   --filter-name "QUAL30" --filter-expression "QUAL < 30.0" \
   --filter-name "FS60" --filter-expression "FS > 60.0" \
   --filter-name "QD2" --filter-expression "QD < 2.0" \
   --filter-name "SOR3" --filter-expression "SOR > 3.0" \
   --filter-name "MQ40" --filter-expression "MQ < 40.0" \
   --filter-name "ReadPosRankSum-8.0" --filter-expression "ReadPosRankSum < -8.0" \
   --filter-name "MQRankSum-12.5" --filter-expression "MQRankSum < -12.5"' %(gatk, mp_fa, mp[i], mp[i][0:n]))
#ASE计数
paths = glob('*_lg_Filter_SNP.vcf')
lg=[]
for path in paths:
    lg.append(os.path.basename(path))
paths = glob('*_mp_Filter_SNP.vcf')
mp=[]
for path in paths:
    mp.append(os.path.basename(path))
lg.sort()
mp.sort()   
for i in range(0,4):
    os.system('%s ASEReadCounter -R %s -I %s.add.bam -V %s -O %s_ASE1.table'%(gatk, lg_fa,lg[i][0:n],lg[i],lg[i][0:n]))
    os.system('%s ASEReadCounter -R %s -I %s.add.bam -V %s -O %s_ASE1.table' %(gatk, mp_fa,mp[i][0:n],mp[i],mp[i][0:n]))

os.system('rm -rf *dup.bam')
os.system('rm -rf *.bai')
'''

import pandas as pd
import sys
def f(df,i,n):
    df=pd.read_table(df)
    df=df[['contig','position','refCount','altCount']]
    df=df.copy()
    df['contig'] = df['contig'].astype('str')
    df['position'] = df['position'].astype('str')
    df.insert(loc=0, column='snp_id', value=df['contig'].str.cat(df['position'],sep='-'))
    del df['contig']
    del df['position']
    df.rename(columns={'refCount': '%s_R%d_refCount' %(i[0:4],n), 'altCount': '%s_R%d_altCount' %(i[0:4],n)}, inplace=True)    
    return df

paths = glob('*_lg_ASE1.table')
lg=[]
for path in paths:
    lg.append(os.path.basename(path))
lg.sort()
df=f(lg[0],lg[0],1)
for i in range(1,len(lg)):
    df=pd.merge(df,f(lg[i],lg[i],i+1),on='snp_id',how='inner')
lg=df
paths = glob('*_mp_ASE1.table')
mp=[]
for path in paths:
    mp.append(os.path.basename(path))
mp.sort()
df=f(mp[0],mp[0],1)
for i in range(1,len(mp)):
    df=pd.merge(df,f(mp[i],mp[i],i+1),on='snp_id',how='inner')
mp=df
lg['a']=lg['snp_id'].apply(lambda x : x.split('-')[0])
lg['b']=lg['snp_id'].apply(lambda x : x.split('-')[1])
lg['c']=lg['snp_id'].apply(lambda x : x.split('-')[1])
mp['a']=mp['snp_id'].apply(lambda x : x.split('-')[0])
mp['b']=mp['snp_id'].apply(lambda x : x.split('-')[1])
mp['c']=mp['snp_id'].apply(lambda x : x.split('-')[1])
lg[['a','b','c']].to_csv('lg_id',sep='\t',header=None,index=None)
mp[['a','b','c']].to_csv('mp_id',sep='\t',header=None,index=None)
del lg['a']
del lg['b']
del lg['c']
del mp['a']
del mp['b']
del mp['c']
lg.to_csv('lg.count',index=None)
mp.to_csv('mp.count',index=None)
lf_gff='/gss1/home/dingbq01/reference/LF10g_v2.0/LF10g_v2.0.gff3'
mp_gff='/gss1/home/dingbq01/reference/Mparishii_genome/Mparg_v2.0.gff3'
os.system('%s sort -chrThenSizeA -i lg_id > LPF1_lg_sort.pos' %bedtools)
os.system('%s sort -chrThenSizeA -i mp_id > LPF1_mp_sort.pos' %bedtools)
os.system('%s intersect -a LPF1_lg_sort.pos -b %s -wb > LPF1_lg_gene.pos' % (bedtools,lf_gff))
os.system('%s intersect -a LPF1_mp_sort.pos -b %s -wb > LPF1_mp_gene.pos' %(bedtools,mp_gff))


df=pd.read_table('LPF1_lg_gene.pos',header=None)
df[0] = df[0].astype('str')
df[1] = df[1].astype('str')
df.insert(loc=0, column='snp_id', value=df[0].str.cat(df[1],sep='-'))
df=df[df[5]=='CDS']
df1=pd.read_table('lg.count',sep=',')
df2=pd.merge(df,df1,on='snp_id',how='inner')
df2=df2[['snp_id',11,'LPF1_R1_refCount','LPF1_R1_altCount','LPF1_R2_refCount','LPF1_R2_altCount','LPF1_R3_refCount','LPF1_R3_altCount','LPF1_R4_refCount','LPF1_R4_altCount']]
df2.columns=[['snp_id','gene_id','LPF1_R1_refCount','LPF1_R1_altCount','LPF1_R2_refCount','LPF1_R2_altCount','LPF1_R3_refCount','LPF1_R3_altCount','LPF1_R4_refCount','LPF1_R4_altCount']]
df2=df2.dropna()
df2=df2[~df2.isin([0])]
df2=df2.dropna()
for a in df2.columns:
    if df2[a].dtypes ==float:
        df2[a]=df2[a].astype('int') 
df2.to_csv('lg_2.count',sep='\t',index=None)

df=pd.read_table('LPF1_mp_gene.pos',header=None)
df[0] = df[0].astype('str')
df[1] = df[1].astype('str')
df.insert(loc=0, column='snp_id', value=df[0].str.cat(df[1],sep='-'))
df=df[df[5]=='CDS']
df1=pd.read_table('mp.count',sep=',')
df2=pd.merge(df,df1,on='snp_id',how='inner')
df2=df2[['snp_id',11,'LPF1_R1_refCount','LPF1_R1_altCount','LPF1_R2_refCount','LPF1_R2_altCount','LPF1_R3_refCount','LPF1_R3_altCount','LPF1_R4_refCount','LPF1_R4_altCount']]
df2.columns=[['snp_id','gene_id','LPF1_R1_refCount','LPF1_R1_altCount','LPF1_R2_refCount','LPF1_R2_altCount','LPF1_R3_refCount','LPF1_R3_altCount','LPF1_R4_refCount','LPF1_R4_altCount']]
df2=df2.dropna()
df2=df2[~df2.isin([0])]
df2=df2.dropna()
for a in df2.columns:
    if df2[a].dtypes ==float:
        df2[a]=df2[a].astype('int') 
df2.to_csv('mp_2.count',sep='\t',index=None)

df=pd.read_table("lg_2.count")
m=df.columns[2:]
for i in m:
    df=df[(df[i]>10)]
df=df.groupby('gene_id').sum()
df[['LPF1_R1_refCount','LPF1_R1_altCount','LPF1_R2_refCount','LPF1_R2_altCount','LPF1_R3_refCount','LPF1_R3_altCount','LPF1_R4_refCount','LPF1_R4_altCount']]=df[['LPF1_R1_refCount','LPF1_R2_refCount','LPF1_R3_refCount','LPF1_R4_refCount','LPF1_R1_altCount','LPF1_R2_altCount','LPF1_R3_altCount','LPF1_R4_altCount']]
df.columns=['LPF1_R1_refCount','LPF1_R2_refCount','LPF1_R3_refCount','LPF1_R4_refCount','LPF1_R1_altCount','LPF1_R2_altCount','LPF1_R3_altCount','LPF1_R4_altCount']
df.to_csv('lg_3.count',sep='\t')

df=pd.read_table("mp_2.count")
m=df.columns[2:]
for i in m:
    df=df[(df[i]>10)]
df=df.groupby('gene_id').sum()
df[['LPF1_R1_refCount','LPF1_R1_altCount','LPF1_R2_refCount','LPF1_R2_altCount','LPF1_R3_refCount','LPF1_R3_altCount','LPF1_R4_refCount','LPF1_R4_altCount']]=df[['LPF1_R1_refCount','LPF1_R2_refCount','LPF1_R3_refCount','LPF1_R4_refCount','LPF1_R1_altCount','LPF1_R2_altCount','LPF1_R3_altCount','LPF1_R4_altCount']]
df.columns=['LPF1_R1_refCount','LPF1_R2_refCount','LPF1_R3_refCount','LPF1_R4_refCount','LPF1_R1_altCount','LPF1_R2_altCount','LPF1_R3_altCount','LPF1_R4_altCount']
df.to_csv('mp_3.count',sep='\t')