#ASE
主要分四部分进行
bidui.py 将RNA-seq的数据比对到参考基因组
guolv.py 使用picard过滤掉测序产生的重复序列
gvcf.py 使用gatk获得snp位点的信息
shuchu.py 使用ASEReadCounter获得每个位点的count数，然后进行数据的过滤，选择每个位点大于 10 个count数的SNP位点
