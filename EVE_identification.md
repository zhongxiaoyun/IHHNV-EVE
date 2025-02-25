# EVE identification

## 1. Species genome alignment to viral protein sequence

```bash
# Construct protein sequence index
#blastx comparison If the protein library is too large and blastx comparison is slow, diamond can be used
query=genome.fa
blastx -query  ${query} -out ${virus}\.blast.out -db ${virus}\.db  -evalue 1e-6 -num_threads 3  -outfmt '6 qseqid qstart qend salltitles evalue qframe pident qcovs sstart send slen'
```

## 2. Convert to 'bed' format

```bash
python2 /public/home/zhongxiaoyun/software/EVE_software/EVE_finder/Scripts/EVEfinder/EVE_finder/Blast_to_Bed3.py ${virus}\.blast.out
```

## 3.Sort 'bed' and merge to find the highest score (E value)

```bash

cat ./*/$i*EVEs.bed > $i\_EVEs.bed
#sort
bedtools sort -i $i\_EVEs.bed > $i\_EVEs.sorted.bed
#merge
bedtools merge -i $i\_EVEs.sorted.bed -c 4,5,6,7,8,9,10,11 -o collapse,collapse,distinct,collapse,collapse,collapse,collapse,collapse > $i\_EVEs.merge.bed
#convert E-value
python change_0_to_e.py $i\_EVEs.merge.bed test
mv test $i\_EVEs.merge.bed
#find best hit
python2 ~/software/EVE_software/EVE_finder/Scripts/EVEfinder/EVE_finder/Top_score_BED2.py $i\_EVEs.merge.bed $i\_EVEs.top.bed
#Due to the irregular E value, 5.90e-8 will be converted to 5.9e-8, so it should be modified in Top_score_BED2.py to retain the two-digit format(sorted_evalues[0],'0.2e).
# But after changing to this script, 0.0 in the original file (merge.bed) will be converted to 0.00e+00 without modifying the script
# Write a script to replace 0.0 with 0.00e+00


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



## 4. Obtain the genome sequence of the species compared

```bash
bedtools getfasta -s -name -fi genome.fa -bed $i\_EVEs.top.bed -fo $i\_EVEs.top.fasta
```

## 5. Compare to the NR database

```bash
# Build the NR database index and compare the sequence to the NR database To get staxids you must add --taxonmap and corresponding index id file when building the index
diamond blastx --db /public/home/zhongxiaoyun/database/02_database/NR/NR/nr_diamond.dmnd  -e 1e-6 --threads 8 -f 6 qseqid qstart qend salltitles evalue qframe pident qcovhsp sstart send slen staxids -q  $i\_EVEs.top.fasta -o $i\_EVEs.top_nr.blastx
#In order to make the result of blastx match the content of top bed (to meet the needs of the script), the first column of blastx needs to be processed, only the Chr:a-b step is required, if this step is missing, there is no detail in the subsequent result
awk -F '::' '{print \$2}' $i\_EVEs.top_nr.blastx |sed  's/(-)//g'  |sed 's/(+)//g' >  $i\_EVEs.top_nr.blastx.out

```

## 6. Filter and compare the results to the NR library and classify them to remove false positives (not known virus sequences, eukaryotic sequences, low complexity sequences).

```bash
bash /public/home/zhongxiaoyun/software/EVE_software/Refine_EVEs_annotation/Refine_EVE_Annotation.sh \ # To manually change TOOl= 'blastx' in this script to TOOL= 'diamond'
                -pipeline_directory /public/home/zhongxiaoyun/software/EVE_software/Refine_EVEs_annotation \
                -tool diamond \ #diamond 
                -VHC_directory /public/home/zhongxiaoyun/software/EVE_software/VHost-Classifier \
                -file_blastx /public/home/zhongxiaoyun/work/04_EVE/01_blastx/$i\_EVEs.top_nr.blastx \
                -file_bed_tophit /public/home/zhongxiaoyun/work/04_EVE/01_blastx/$i\_EVEs.top.bed \
                -output_directory /public/home/zhongxiaoyun/work/04_EVE/01_blastx/Result/$i\_Output_directory \ #If the same sample is executed repeatedly, the file will not be overwritten and needs to be deleted and executed again
                -taxonkit_exe taxonkit
  
```

## 7.Filter EVE

```bash
#Sequences that are aligned to non-viruses are filtered out based on database annotations that are aligned to NR
awk '$9>0{print}' CompleteTable_classified.table > CompleteTable_classified.filter
#Filter unified source EVE based on virus taxid
```

