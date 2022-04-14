# Tn5 _M. maripaludis_ S2 library analysis 
## Tn5 Library from [Genome-scale analysis of gene function in the hydrogenotrophic methanogenic archaeon Methanococcus maripaludis, Sarmiento et. al. 2013](https://www.pnas.org/doi/10.1073/pnas.1220225110) 

## Parameters used for trimming transposon sequencing  
## Performed using cutadapt 1.18 with Python 3.7.2
```
cutadapt -j 0 -m 30 --pair-filter=any --nextseq-trim=30 -u 27 -U 27 -g GGTTGAGATGTGTATAAGAGACNG -g CTACAAGAGCGGTGAG -a ACTCACCGCTCTTGTAG -a CTGTCTCTTATA -A CTGTCTCTTATA -o S5_R1.fastq.gz -p S5_R2.fastq.gz S5_R1.fastq.gz S5_R2.fastq.gz
```
### cutadapt use all available cores
-j 0 
### discard all reads <30 bp 
-m 30 
### if either R1 or R2 is <30 bp discard both
--pair-filter=any
### remove bp from 3' end until reach base with quality score of 30 or greater
--nextseq-trim=30
### remove first 27 bases from both R1 and R2 (removes tn seq) so first base of both R1 and R2 should now be genomic sequence immediately following tn insertion
-u 27 -U 27
### there is still some transposon sequence present at 5' end of R1 reads (likely do to redundcancy in seq of library prep primers
#### remove remaining transposon seq
-g GGTTGAGATGTGTATAAGAGACNG
### both R1 has sequence contaminantion of illumina primers that are annealed following shearing of the gDNA - remove from 3' and 5' end
-g CTACAAGAGCGGTGAG -a ACTCACCGCTCTTGTAG
### remove Nextera transposase sequence from R1 and R2 reads
-a CTGTCTCTTATA -A CTGTCTCTTATA


## map remaining reads to _M. maripaludis_ S2 genome using bowtie2 v2.3.4.1
```
bowtie2 -x MM_index -1 ~/tnseq01/concat/trim_02/S5_R1.fastq.gz -2 ~/tnseq01/concat/trim_02/S5_R2.fastq.gz -S S5.sam
```
## check quality of mapping, samtools v.1.9
```
samtools stats S5.sam > S5_stats
```
## pull reads that only map to a single location 
```
awk '$5 > 1' S5.sam > S5_uniq.sam
```
# following one liner and tabulation script ammended from the incredibly helpful [Jon Badalmenti](https://github.com/jbadomics/tnseq)
## generate hit file that compiles the insertion counts 
```
/panfs/roc/groups/2/kcosta/day00094/bin/bioawk/bioawk -c sam '{ if($flag==99) print $qname,$rname,$flag,$pos,$pos+1,length($seq); if($flag==83) print $qname,$rname,$flag,$pos,$pos+length($seq)-1,length($seq) }' ~/tnseq01/concat/trim_02/bowtie2/S5_uniq.sam | cut -f2,5 | sort -g -k2 | uniq -c | sed 's/^ *//' | awk -v OFS='\t' '{ print $2,$3,$1}' > S5_uniq_bowtie2.hits.txt
```
```
cd tabulate_insertions_S5/
```
## must use python2, I used v2.7.12_anaconda4.2, biopython v.1.61
```
python tabulate_insertions8.py MMS2.gbk ~/tnseq01/concat/trim_02/bowtie2/S5_uniq_bowtie2.hits.txt 'MMP' 0.05 0.20
```
