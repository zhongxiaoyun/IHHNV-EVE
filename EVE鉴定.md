# 二、EVE鉴定

## 1.物种基因组比对到病毒蛋白序列

```bash
#构建蛋白序列索引
#blastx比对  如果蛋白库太大，blastx比对的慢 可以用diamond
query=genome.fa
blastx -query  ${query} -out ${virus}\.blast.out -db ${virus}\.db  -evalue 1e-6 -num_threads 3  -outfmt '6 qseqid qstart qend salltitles evalue qframe pident qcovs sstart send slen'
```

## 2.转换为BDE格式

```bash
python2 /public/home/zhongxiaoyun/software/EVE_software/EVE_finder/Scripts/EVEfinder/EVE_finder/Blast_to_Bed3.py ${virus}\.blast.out

这个脚本里面容易把id中的逗号给换到下一行，所以在脚本里添加一行把逗号替换成没有 name.replace(",","")
根据坐标判断正负 脚本里写错了进行改正，将if unsorted[0]>unsorted[1] :改成if se[0]>se[1]: strand=”-“
```

## 3.BED排序及合并,找到比分最高的（E值）

```bash
#同一个病毒的多个比对文件合并，也可以在一开始合并蛋白序列
cat ./*/$i*EVEs.bed > $i\_EVEs.bed
#排序
bedtools sort -i $i\_EVEs.bed > $i\_EVEs.sorted.bed
#合并
bedtools merge -i $i\_EVEs.sorted.bed -c 4,5,6,7,8,9,10,11 -o collapse,collapse,distinct,collapse,collapse,collapse,collapse,collapse > $i\_EVEs.merge.bed
#转换E值
python change_0_to_e.py $i\_EVEs.merge.bed test
mv test $i\_EVEs.merge.bed
#找到最佳比对
python2 ~/software/EVE_software/EVE_finder/Scripts/EVEfinder/EVE_finder/Top_score_BED2.py $i\_EVEs.merge.bed $i\_EVEs.top.bed
#由于E值的不规范，5.90e-8 会被转换成5.9e-8所以要在Top_score_BED2.py中修改成保留两位数format(sorted_evalues[0],'0.2e)
#但是改成这个脚本之后原文件（merge.bed）中的0.0会被转换成0.00e+00 不修改脚本
#写一个脚本将0.0替换为0.00e+00 


```

```python
vi change_0_to_e.py
#!/usr/bin/python

import sys
f1=open(sys.argv[1],'r')
f2=open(sys.argv[2],'w')

for i in f1:
        i=i.strip()
        l=i.split('\t')
        if '0.0' not in l[4]:
                f2.write(i+"\n")
        else:
                a=l[4].replace('0.0','0.00e+00')
                f2.write("%s\t%s\t%s\n"%('\t'.join(l[0:4]),a,'\t'.join(l[5:])))
```



## 4.获得比对上的物种基因组序列

```bash
bedtools getfasta -s -name -fi genome.fa -bed $i\_EVEs.top.bed -fo $i\_EVEs.top.fasta
```

## 5.比对到NR数据库

```bash
#构建NR数据库索引，将序列比对到NR数据库 要获得staxids 必须在构建索引时加---taxonmap 并对应索引id文件
diamond blastx --db /public/home/zhongxiaoyun/database/02_database/NR/NR/nr_diamond.dmnd  -e 1e-6 --threads 8 -f 6 qseqid qstart qend salltitles evalue qframe pident qcovhsp sstart send slen staxids -q  $i\_EVEs.top.fasta -o $i\_EVEs.top_nr.blastx
##为了使blastx的结果和top bed的内容匹配上（满足脚本的需要），需要对blastx的第一列进行处理，只需要留下Chr:a-b 这一步必需，如果少这一步后续结果没有详细信息
awk -F '::' '{print \$2}' $i\_EVEs.top_nr.blastx |sed  's/(-)//g'  |sed 's/(+)//g' >  $i\_EVEs.top_nr.blastx.out

```

## 6.筛选比对到NR库的结果及分类，去除假阳性（不是已知病毒序列、真核生物序列、低复杂度序列）

```bash
bash /public/home/zhongxiaoyun/software/EVE_software/Refine_EVEs_annotation/Refine_EVE_Annotation.sh \ #要手动将这个脚本里的TOOl=‘blastx‘改成TOOL=’diamond‘ 已经改过了就行了 
                -pipeline_directory /public/home/zhongxiaoyun/software/EVE_software/Refine_EVEs_annotation \
                -tool diamond \ #一般用diamond 
                -VHC_directory /public/home/zhongxiaoyun/software/EVE_software/VHost-Classifier \
                -file_blastx /public/home/zhongxiaoyun/work/04_EVE/01_blastx/$i\_EVEs.top_nr.blastx \
                -file_bed_tophit /public/home/zhongxiaoyun/work/04_EVE/01_blastx/$i\_EVEs.top.bed \
                -output_directory /public/home/zhongxiaoyun/work/04_EVE/01_blastx/Result/$i\_Output_directory \ #如果是重复执行同一个样本，这个文件不会被覆盖，需要删除后重新执行
                -taxonkit_exe taxonkit
  
```
