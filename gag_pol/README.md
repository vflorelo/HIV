Assembly of gag/pol fragments from NGS
======================================
Amplicons were obtained by RT-PCR using primers (#SAM) followed by NGS sequencing using (#Caro & #Paco).
We followed the procedure from Shambhu *et. al.* [[ref](https://pubmed.ncbi.nlm.nih.gov/27448822/)], with slight modifications.

Workspace setup
----------------
The workspace consisted of the following structure:

`${workdir}/HIV/assembly`,`${workdir}/HIV/clusters`,`${workdir}/HIV/mappings`,`${workdir}/HIV/qc`,`${workdir}/HIV/reads` and `${workdir}/HIV/tree`. Reads were placed in individual directories inside the reads subdirectory.
Quality assesment
-----------------
NGS-quality was inspected using fastQC [[ref](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)]. Plots and data can be accessed [[here](qc.md)].
```bash
cd ${workdir}/HIV/qc
ln -s ../reads/*/*.gz .
fastqc *.gz
```
Quality and length based trimming
---------------------------------
We used [sickle](https://github.com/najoshi/sickle) to discard reads with a phred score of 28 or lower and a length of less than 190 nucleotides. In order to optimize the assembly process, we removed redundant reads using [fastuniq](http://sourceforge.net/projects/fastuniq/) [[ref](https://pubmed.ncbi.nlm.nih.gov/23284954/)].
```bash
cd ${workdir}/HIV/reads
for sample in $(ls)
do
  cd $sample;
  sample_num=$(ls | grep gz | cut -d_ -f2 | sort -V | uniq | head -n1);
  sickle pe -f ${sample}_${sample_num}_L001_R1_001.fastq.gz -r ${sample}_${sample_num}_L001_R2_001.fastq.gz -o ${sample}_R1.fastq.gz -p ${sample}_R2.fastq.gz -s ${sample}_failed.fastq.gz -q 30 -l 226 -t sanger -g
  zcat ${sample}_R1.fastq.gz > ${sample}_R1.fastq
  zcat ${sample}_R2.fastq.gz > ${sample}_R2.fastq
  echo -e "${sample}_R1.fastq\n${sample}_R2.fastq" > list
  fastuniq -i list -t q -o ${sample}_R1_nr.fastq -p ${sample}_R2_nr.fastq
  cd ..
done
```
Fragment assembly
-----------------
Major HIV-subtypes were assembled using [iva](https://github.com/sanger-pathogens/iva) [[ref](https://pubmed.ncbi.nlm.nih.gov/25725497/)]
```bash
cd ${workdir}/HIV/assembly
mkdir iva
cd iva
ln -s ../../reads/*/*nr.fastq .
ls | grep fastq | cut -d_ -f1 | sort -V | uniq > sample_list
for sample in $(cat sample_list)
do
  iva -f ${sample}_R1_nr.fastq -r ${sample}_R2_nr.fastq -v --threads 12 ${sample}
done
```
After fragment assembly was completed, sequences were merged into a single file. Fragments that were assembled in more than a contiguous sequence (contig) were merged using the [seaview](http://pbil.univ-lyon1.fr/software/seaview3) consensus tool. Complete contigs were aligned to the HIV HXB2 genomic sequence corresponding to the amplified region from nucleotides 790 to 2549. The alignment was trimmed so all sites were informative for phylogenetic reconstruction
```bash
cd ${workdir}/HIV/assembly
cd iva
for sample in $(cat sample_list)
do
  perl -pe "s/\>contig\.0000/\>$sample\_ctg0/" $sample/contigs.fasta >> contigs.fasta
done
esearch -db nucleotide -query "K03455 [accn]" | efetch -format fasta > K03455.fasta
seqret K03455.fasta -sbegin 790 -send 2549 fasta::stdout >> contigs.fasta
```
Fragment assembly
-----------------
Minor HIV-subtypes were assembled using [atram2](https://github.com/juliema/aTRAM) [[ref](https://pubmed.ncbi.nlm.nih.gov/29881251/)]. Contig assembly stats were parsed from the fasta headers
```bash
cd ${workdir}/HIV/assembly
mkdir atram2
cd atram2
ln -s ../../reads/*/*nr.fastq .
ls | grep fastq | cut -d_ -f1 | sort -V | uniq > sample_list
esearch -db nucleotide -query "K03455 [accn]" | efetch -format fasta > K03455.fasta
for sample in $(cat sample_list)
do
  atram_preprocessor.py --end-1 ${sample}_R1_nr.fastq --end-2 ${sample}_R2_nr.fastq --blast-db ${sample} --cpus 12 --fastq;
  atram.py --blast-db ${sample} --query K03455.fasta --output-prefix ${sample} --assembler trinity --cpus 12 --iterations 5
  grep \> ${sample}.${sample}_K03455.filtered_contigs.fasta | perl -pe 's/\>\>/\>$sample\t/;s/\ .*TRINITY/\t/;s/\ .*\=/\t/' | awk 'BEGIN{FS="\t"}{print $1 FS $2 FS $1$3 FS $4}' > contig_stats.tsv
done
for sample in $(cut -f1 contig_stats.tsv | sort -V | uniq)
do
  sample_datablock=$(grep -w ^$sample contig_stats.tsv)
  num_contigs=$(echo "$sample_datablock" | wc -l)
  for i in $(seq 1 $num_contigs)
  do
    contig_id=$(echo "$sample_datablock" | tail -n+$i | head -n1 | cut -f2)
    contig_name=$(echo "$sample_datablock" | tail -n+$i | head -n1 | cut -f3)
    contig_score=$(echo "$sample_datablock" | tail -n+$i | head -n1 | cut -f4)
    contig_iteration=$(echo "$contig_id" | cut -d_ -f1)
    seqret ${sample}.${sample}_K03455.filtered_contigs.fasta:$contig_id fasta::stdout | perl -pe "s/\>.*/\>$contig_name\_$contig_iteration\ $contig_score/"
  done
done > contigs.fasta
```
We then selected contigs aligning to the gag region or the pol region using BLASTn [[ref](https://pubmed.ncbi.nlm.nih.gov/10890397/)], only contigs fully matching the corresponding regions were kept. The resulting sequences were aligned using [muscle](https://drive5.com/muscle/) [[ref](https://pubmed.ncbi.nlm.nih.gov/15318951/)] via seaview
```bash
makeblastdb -in contigs.fasta -title contigs.fasta -dbtype nucl
blastn -query K03455_gag.fasta -db contigs.fasta -max_target_seqs 1000 -out gag.blastn -outfmt "6 qseqid sseqid qlen slen qstart sstart qend send score evalue length positive"
#Column 7 corresponds to the 3' end of the gag region, contigs mapping after position 25 or before position 1400 were discarded
awk 'BEGIN{FS="\t"}{if($5<=25 && $7>=1400){print $AF}}' gag.blastn > gag_filtered.blastn
for seq_id in $(cut -f2 gag_filtered.blastn | sort -V | uniq)
do
  seq_start=$(grep -w $seq_id gag.blastn | cut -f8 | sort -n | head -n1 )
  seq_end=$(grep -w $seq_id gag.blastn | cut -f8 | sort -n | tail -n1 )
  seqret contigs.fasta:$seq_id -sbegin $seq_start -send $seq_end fasta::stdout
done > gag.fasta
blastn -query K03455_pol.fasta -db contigs.fasta -max_target_seqs 1000 -out pol.blastn -outfmt "6 qseqid sseqid qlen slen qstart sstart qend send score evalue length positive"
for seq_id in $(cut -f2 pol.blastn | sort -V | uniq)
do
  seq_start=$(grep -w $seq_id pol.blastn | cut -f8 | sort -n | head -n1 )
  seq_end=$(grep -w $seq_id pol.blastn | cut -f8 | sort -n | tail -n1 )
  seqret contigs.fasta:$seq_id -sbegin $seq_start -send $seq_end fasta::stdout
done > pol.fasta
```
Clustering and reduction of contigs
-----------------------------------
In order to reduce the number of sequences to analyse we constructed clusters of contigs, these clusters consisted of sequences at least 99% identical, cd-hit was employed for such purpose
```bash
cd HIV/clusters
cp ${workdir}/HIV/assembly/atram2/gag.fasta .
cp ${workdir}/HIV/assembly/atram2/pol.fasta .
for sample in $(grep \> gag.fasta | cut -d_ -f1 | cut -d\> -f2 | sort -V | uniq)
do
  contig_list=$(grep ^${sample}_ gag.fasta | cut -d\> -f2 | sort -V | uniq)
  for contig in $contig_list
  do
    seqret gag.fasta:$contig fasta::stdout
  done > $sample.gag.tmp.fasta
  cd-hit -i $sample.gag.tmp.fasta -o $sample.gag.fasta -T 12 -c 0.99
done

for sample in $(grep \> pol.fasta | cut -d_ -f1 | cut -d\> -f2 | sort -V | uniq)
do
  contig_list=$(grep ^${sample}_ pol.fasta | cut -d\> -f2 | sort -V | uniq)
  for contig in $contig_list
  do
    seqret pol.fasta:$contig fasta::stdout
  done > $sample.pol.tmp.fasta
  cd-hit -i $sample.pol.tmp.fasta -o $sample.pol.fasta -T 12 -c 0.99
done
```
Abundance of transcripts
-----------------------------------
Once we constructed the clusters of nonredundant sequences, we proceeded to quantify such sequences by mapping the reads using [hisat2](https://github.com/DaehwanKimLab/hisat2) [[ref](https://pubmed.ncbi.nlm.nih.gov/25751142/)], alignments were processed using [samtools](https://github.com/samtools/samtools) and abundances were estimated using [htseq](https://htseq.readthedocs.io/en/release_0.11.1/index.html)
Summarized results for [gag transcripts](gag_counts.md) and [pol transcripts](pol_counts.md)
```bash
cd ${workdir}/HIV/mapping
mkdir gag
cd gag
cp ${workdir}/HIV/clusters/*.gag.fasta .
ln -s ${workdir}/HIV/reads/*/*R1.fastq.gz .
ln -s ${workdir}/HIV/reads/*/*R2.fastq.gz .
for sample in $(ls | grep fasta | cut -f1 -d\. | sort -V | uniq)
do
  infoseq ${sample}.pol.fasta | awk 'BEGIN{sep="\t"}{print $1 sep "atram2" sep "contig" sep 1 sep $2 sep "." sep "+" sep "." sep "contig_id=\""$1"\"" }' > ${sample}.gff
  hisat2-build -p 12 ${sample}.gag.fasta ${sample}
  hisat2 -p 12 -x ${sample} -1 ${sample}_R1.fastq.gz -2 ${sample}_R2.fastq.gz --very-sensitive > ${sample}.sam
  samtools view -f 3 -F 256 -q 60 -h -b -o ${sample}.tmp.bam -@ 12 ${sample}.sam
  samtools sort -l 9 -@ 12 -o ${sample}.bam ${sample}.tmp.bam
  samtools index ${sample}.bam
  htseq-count -n 12 -r pos -t contig -i contig_id ${sample}.bam ${sample}.gff > ${sample}.tsv
	contig_list=$(cut -f2 ${sample}.tsv | sort -V | uniq )
	total_reads=$(awk 'BEGIN{FS="\t"}{sum+=$3}END{print sum}' ${sample}.tsv)
	for contig in $contig_list
	do
		contig_reads=$(grep -w $contig ${sample}.tsv| cut -f3)
		echo -e "$sample\t$contig\t$contig_reads\t$total_reads" | awk 'BEGIN{FS="\t"}{print $1 FS $2 FS $3 FS $3/$4}'
	done
done >> gag_counts.tsv

cd ${workdir}/HIV/mapping
mkdir pol
cd pol
cp ${workdir}/HIV/clusters/*.pol.fasta .
ln -s ${workdir}/HIV/reads/*/*R1.fastq.gz .
ln -s ${workdir}/HIV/reads/*/*R2.fastq.gz .
for sample in $(ls | grep fasta | cut -f1 -d\. | sort -V | uniq)
do
  infoseq ${sample}.pol.fasta | awk 'BEGIN{sep="\t"}{print $1 sep "atram2" sep "contig" sep 1 sep $2 sep "." sep "+" sep "." sep "contig_id=\""$1"\"" }' > ${sample}.gff
  hisat2-build -p 12 ${sample}.pol.fasta ${sample}
  hisat2 -p 12 -x ${sample} -1 ${sample}_R1.fastq.gz -2 ${sample}_R2.fastq.gz --very-sensitive > ${sample}.sam
  samtools view -f 3 -F 256 -q 60 -h -b -o ${sample}.tmp.bam -@ 12 ${sample}.sam
  samtools sort -l 9 -@ 12 -o ${sample}.bam ${sample}.tmp.bam
  samtools index ${sample}.bam
  htseq-count -n 12 -r pos -t contig -i contig_id ${sample}.bam ${sample}.gff > ${sample}.tsv
	contig_list=$(cut -f2 ${sample}.tsv | sort -V | uniq )
	total_reads=$(awk 'BEGIN{FS="\t"}{sum+=$3}END{print sum}' ${sample}.tsv)
	for contig in $contig_list
	do
		contig_reads=$(grep -w $contig ${sample}.tsv| cut -f3)
		echo -e "$sample\t$contig\t$contig_reads\t$total_reads" | awk 'BEGIN{FS="\t"}{print $1 FS $2 FS $3 FS $3/$4}'
	done
done >> pol_counts.tsv
```
