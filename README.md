# Pac-Bio-AND-Hybrid-assembly
a collection of codes for genomes assembly and cleaning compiled and designed by Victor Jimenez (vr.jimenez.vs@gmail.com

## 1. runing CANU (including, trimming, reads correction and assembly) ##
```r

mkdir pac_bio_contigs/ ; 

for seq in *.gz
do
r1=$(basename $seq .subreads.fastq.gz)
canu -p ${r1} -d ${r1} genomeSize=1500000 corThreads=30 -pacbio-raw $seq ;
mv ${r1}/${r1}.contigs.fasta ${r1}/${r1}.contigs.gfa . ;
mv ${r1}.contigs.fasta ${r1}.contigs.gfa pac_bio_contigs/ ;
done; 

ls ;
```

## 2. check read length ##
```r
## contig length  ##
zcat 2_PB.subreads.fastq.gz | grep "^@" | cut -d"=" -f4 | cut -d " " -f1 | uniq -d | sort -n -r | head -n 10 ; 
grep ">" 6_PB.contigs.fasta | cut -d "=" -f2 | cut -d " " -f1 | uniq | sort -n -r ;
ls ;

## install circlator ##
conda create -n circlator ; 
conda activate circlator ;
 
sudo apt-get install apt-file
sudo apt-file update
apt-file search libcrypto.so.1.0.0
```

## 3. circularize contigs with circlator (check: reads prefix) ##
```r
mkdir circlator.res/ ; 

for sample in *.gz
do
r1=$(basename $sample .correctedReads.fasta.gz)
circlator all --threads 30 --verbose ${r1}.contigs.fasta ${r1}.correctedReads.fasta.gz ${r1}.circlator.out ; 
mv ${r1}.circlator.out/06.fixstart.fasta ${r1}.circlator.out/${r1}.06.fixstart.fasta ; 
cp ${r1}.circlator.out/${r1}.06.fixstart.fasta circlator.res/ ; 
done ; 

cd circlator.res/ ; 
ls -lh ; 
```

## 4. correct PACBIO-contigs (check: genome names and prefix, threads, reads prefix) ##
```r
mkdir pilon_contigs/ ; 

for s1 in *.fasta
do
r1=$(basename $s1 .06.fixstart.fasta)

bwa index $s1 ; 
bwa mem -t 30 ${r1}.06.fixstart.fasta ${r1}.raw_1.fastq.gz ${r1}.raw_2.fastq.gz | samtools sort > ${r1}.pilon.aln.bam ;
samtools index ${r1}.pilon.aln.bam ; 
samtools faidx ${r1}.06.fixstart.fasta ;
conda activate pilon ;
pilon --genome ${r1}.06.fixstart.fasta --frags ${r1}.pilon.aln.bam --output ${r1}.pilon --outdir ${r1}.pilon --fix all --mindepth 0.5 --changes --verbose --threads 30 ;
cp ${r1}.pilon/${r1}.pilon.fasta . ;
mv ${r1}.pilon.fasta pilon_contigs/ ;
rm ${r1}.pilon.aln.bam ;
done ; 

ls ;
```

## 5. correct SPADES-contigs (check: genome names and prefix, threads, reads prefix) ##
```r

mkdir pilon_contigs/ ;
mkdir pilon_intermediates/ ; 

for s1 in *.fasta
do
r1=$(basename $s1 _spades_scaffolds.fasta)

bwa index $s1 ; 
bwa mem -t 30 ${r1}_spades_scaffolds.fasta ${r1}_f_paired.fq.gz ${r1}_r_paired.fq.gz | samtools sort -@ 30 > ${r1}.pilon.aln.bam ;
samtools index -@ 30 ${r1}.pilon.aln.bam ; 
samtools faidx ${r1}_spades_scaffolds.fasta ;
conda activate pilon ;
pilon --genome ${r1}_spades_scaffolds.fasta --frags ${r1}.pilon.aln.bam --output ${r1}.pilon --outdir ${r1}.pilon --fix all --mindepth 0.5 --changes --verbose --threads 30 ;
cp ${r1}.pilon/${r1}.pilon.fasta . ;
mv ${r1}.pilon.fasta pilon_contigs/ ;
rm ${r1}.pilon.aln.bam ;
mv ${r1}.pilon/ pilon_intermediates/ ;
rm *.fasta.* ;
done ;

rm *.bam.bai ; 

ls ;
```
