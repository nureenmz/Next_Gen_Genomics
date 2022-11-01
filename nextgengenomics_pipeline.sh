## Week 2 - FastQC
# Create symbolic link to file ie. a shortcut
ln -s /path/to/file /path/copied/to
# fastqc on .gz file, report runtime
time fastqc *.gz
# Check for reads containing the adapter sequence to be trimmed
# for Nextera, Illumina PRep, Illumina PCR kits
zgrep CTGTCTCTTATACACATC *.gz | wc
# Install python cutadapt
python3 -m pip install --user --upgrade cutadapt 
echo PATH=~/.local/bin/:$PATH >> ~/.profile 
source ~/.profile
# Trim adapters using cutadapt
cutadapt -q 30 -m 50 -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT --trim-n -o ERR537186_1_trimmed.fastq -p ERR537186_2_trimmed.fastq ERR537186_1.fastq.gz ERR537186_2.fastq.gz
	# (-q) trim LQ base from 3' end before adapter removal
	# (-m) Discard trimmed reads shorter than 50
	# (-a) Sequence of the adapter ligated to the 3' end of first read
	# (-A) Sequence of the adapter ligated to the 3' end of SECOND read
	# (--trim-n) REmove flanking N bases from each read
	# -o first output file
	# -p second output file
# Run fastqc in parallel for trimmed files
fastqc ERR537186_1_trimmed.fastq.gz ERR537186_2_trimmed.fastq.gz



## Week 5 - Assembly
# Create NGG env (prepared in server) using conda
/localdisk/software/anaconda3/bin/conda env create -n NGG2 --file /localdisk/data/NGG/conda_envs/NGG2.yml
# Activate the env (2 options)
source /localdisk/software/anaconda3/bin/activate NGG2
source /localdisk/software/anaconda3/bin/deactivate NGG2 # to leave env
# Create symbolic link to E.coli sequence files
ln -s /localdisk/data/NGG/ecoli_data/6991_1.fastq ./6991_1.fastq
ln -s /localdisk/data/NGG/ecoli_data/6991_2.fastq ./6991_2.fastq
# fastqc on .fastq files
fastqc *.fastq
firefox *.html &
# Trim sequences using cutadapt
	# CGGTTCAGCAGGAATGCCGAGCAGTAGATCGGAAGAGCGGTTCAG  overrepped seq in seq 1
	# expr length "string" = 45 bases = length of the reads??
	# CGGCATTCCTGCTGAACCGAGCAGTAGATCGGAAGAGCGTCGTGT overrepped seq in seq 2
cutadapt -a CGGTTCAGCAGGAATGCCGAGCAGTAGATCGGAAGAGCGGTTCAG -A CGGCATTCCTGCTGAACCGAGCAGTAGATCGGAAGAGCGTCGTGT --trim-n -o 6991_1_trimmed.fastq -p 6991_2_trimmed.fastq 6991_1.fastq 6991_2.fastq
# Assemble genome using Spades 'standard' mode
spades.py
time spades.py -1 6991_1.fastq -2 6991_2.fastq -t 2 -o new_assembly
ls new_assembly
# Assess the assembly for the number of contigs, the span, and N50
# Process more than one assembly at one time
cd
wget https://www.dropbox.com/s/m50fs79vs0l9ezi/scaffold_stats.pl 
chmod +x scaffold_stats.pl
# Report data for all contigs greater than 100 bases
~/scaffold_stats.pl -f new_assembly/contigs.fasta -t 100
# Visualize DB graph using Bandage
Bandage &
# Tune Spades parameters in 'careful' mode
time spades.py --careful -1 6991_1.fastq -2 6991_2.fastq -t 2 -o careful_assembly
	# --careful reduces the number of mismatches and short indels, but increases run time
	# then run scaffolding and bandage
### Mapping to a reference
# Index the assembled genome using BWA
mkdir aligned_reads
cp XX/contigs.fasta aligned_reads/
bwa index -a is aligned_reads/contigs.fasta
# Align the reads to the genome by find the suffix array coordinates of the single-end reads
bwa aln aligned_reads/contigs.fasta 6991_1.fastq > aligned_reads/6991_1.sai
bwa aln aligned_reads/contigs.fasta 6991_2.fastq > aligned_reads/6991_2.sai
# Generate paired-end alignments in SAM format
bwa sampe aligned_reads_careassem/contigs.fasta \
aligned_reads_careassem/6991_1.sai aligned_reads_careassem/6991_2.sai \
6991_1.fastq 6991_2.fastq \
> aligned_reads_careassem/6991.sam
# Convert SAM to BAM, sort and index BAM
samtools view -bS aligned_reads/6991.sam > aligned_reads/6991.bam
samtools sort aligned_reads/6991.bam > aligned_reads/6991_sorted.bam
samtools index aligned_reads/6991_sorted.bam
igv &
# Genomes > load genome > contigs.fasta
# File > 6991_sorted.bam > select contig and explore the alignment


## Week 6 - Long-read Assembly
# Create NGG env (prepared in server) using conda
/localdisk/software/anaconda3/bin/conda env create -n NGG3 --file /localdisk/data/NGG/conda_envs/NGG3.yml
# Activate the env
source /localdisk/software/anaconda3/bin/activate NGG3
# Create symbolic link to C.elegans data files
ln -s /localdisk/data/NGG/pREADS_data-corrected.fasta.gz .
ln -s /localdisk/data/NGG/CMONO_bothruns_allreads_nanopore.fastq.gz .
#
# QC of long-read data using NanoPlot, time the run
time NanoPlot --fasta pREADS_data-corrected.fasta.gz -o PBio_nanoplot --loglength --N50
firefox PBio_nanoplot/NanoPlot-report.html &
# Make initial assembly with wtdbg2, Run C.elegans 40-fold coverage dataset
# "-x rs" specifies sequencing technology, "-g100m" specifies genome size
# This generates an FDG of the 256 base bins
mkdir Celegans_wtdbg2
time wtdbg2 -i pREADS_data-corrected.fasta.gz \
-o Celegans_wtdbg2/Celegans_PacBio \
-g 100m -t 4 -x rs
# .lay.gz is a layout file that describes how the FDG nodes might be turned into contigs
# Generate the contigs using wtpoa-cns script
time wtpoa-cns -t 8 -i Celegans_wtdbg2/Celegans_PacBio.ctg.lay.gz \
-fo Celegans_wtdbg2/Celegans_PacBio.ctg.fa
# Assess assembly using scaffold_stats.pl
~/scaffold_stats.pl -f Celegans_wtdbg2/Celegans_PacBio.ctg.fa
# Modify assembly by removing -x flag and generate contigs
mkdir Celegans_wtdbg2_2
time wtdbg2 -i pREADS_data-corrected.fasta.gz \
-o Celegans_wtdbg2_2/Celegans_none \
-g 100m -t 4
time wtpoa-cns -t 20 -i Celegans_wtdbg2_2/Celegans_none.ctg.lay.gz \
-fo Celegans_wtdbg2_2/Celegans_none.ctg.fa
# Assess both specified and unspecified assemblies
~/scaffold_stats.pl -f Celegans_wtdbg2/Celegans_PacBio.ctg.fa -f Celegans_wtdbg2_2/Celegans_none.ctg.fa
#
# Do QC of C.monodelphis minION data
pip install NanoPlot --upgrade
time NanoPlot --fastq CMONO_bothruns_allreads_nanopore.fastq.gz -o minION_nanoplot --loglength --N50 -t 8
firefox minION_nanoplot/NanoPlot-report.html &
# Make initial assembly, generate contigs, assess assembly
mkdir Cmonodelphis_wtdbg2
time wtdbg2 -i CMONO_bothruns_allreads_nanopore.fastq.gz \
-o Cmonodelphis_wtdbg2/Cmonodelphis_minION \
-g 115m -t 20 -x ont
time wtpoa-cns -t 20 -i Cmonodelphis_wtdbg2/Cmonodelphis_minION.ctg.lay.gz \
-fo Cmonodelphis_wtdbg2/Cmonodelphis_minION.ctg.fa
~/scaffold_stats.pl -f Cmonodelphis_wtdbg2/Cmonodelphis_minION.ctg.fa
# Make another assembly of minION data, but modify the assembler parameters (from ont to rs), generate contigs, and compare assemblies
mkdir Cmonodelphis_wtdbg2_2
time wtdbg2 -i CMONO_bothruns_allreads_nanopore.fastq.gz \
-o Cmonodelphis_wtdbg2_2/Cmonodelphis_PacBio \
-g 115m -t 20 -x rs \
| time wtpoa-cns -t 20 -i Cmonodelphis_wtdbg2_2/Cmonodelphis_PacBio.ctg.lay.gz \
-fo Cmonodelphis_wtdbg2/Cmonodelphis_PacBio.ctg.fa \
| ~/scaffold_stats.pl -f Cmonodelphis_wtdbg2/Cmonodelphis_minION.ctg.fa -f Cmonodelphis_wtdbg2/Cmonodelphis_PacBio.ctg.fa \
> scaffstats_cmonodelphis.txt


## Week 8 - Metagenome assembly
# Create NGG env (prepared in server) using conda
 /localdisk/software/anaconda3/bin/conda env create -n NGG4 --file /localdisk/data/NGG/conda_envs/NGG4.yml
# Activate the env
source /localdisk/software/anaconda3/bin/activate NGG4
# Create symbolic link to sequence files, taxon ID map, genomes, and blobtools
ln -s /localdisk/software/blobtools-1.1/data/nodesDB.txt ./nodesDB.txt
ln -s /localdisk/data/NGG/ERR260505_1.fastq.gz ./ERR260505_1.fastq.gz
ln -s /localdisk/data/NGG/ERR260505_2.fastq.gz ./ERR260505_2.fastq.gz
ln -s /localdisk/data/NGG/NexteraPE-PE.fa ./NexteraPE-PE.fa
ln -s /localdisk/data/NGG/taxid_map ./taxid_map
ln -s /localdisk/data/NGG/genomes.fasta ./genomes.fasta
##
# fastqc of raw data
fastqc -t 2 ERR260505_1.fastq.gz ERR260505_2.fastq.gz
firefox ERR260505_1_fastqc.html ERR260505_2_fastqc.html &
# Trim and filter raw reads using Skewer using specified Fasta file
mkdir trimmed_reads_skewer
skewer-0.2.2-linux-x86_64 -n -Q 20 -l 75 -t 2 -m any \
-x NexteraPE-PE.fa ERR260505_1.fastq.gz ERR260505_2.fastq.gz \
-o trimmed_reads_skewer/ERR260505
# Repeat FastQC to check skewer output
fastqc -t 2 trimmed_reads_skewer/*.fastq
##
# SPAdes assembly in std genome mode [-m ?]
spades.py --only-assembler -m 30 -t 2 \
-1 trimmed_reads_skewer/*pair1.fastq \
-2 trimmed_reads_skewer/*pair2.fastq \
-o spades_standard
# SPAdes assembly in metagenome mode
spades.py --meta --only-assembler -m 30 -t 4 \
-1 trimmed_reads_skewer/*pair1.fastq \
-2 trimmed_reads_skewer/*pair2.fastq \
-o spades_meta
# Examine output statistics using scaffold_stats.pl
~/scaffold_stats.pl -f spades_standard/scaffolds.fasta -t 100 -h > stats_spades.txt
~/scaffold_stats.pl -f spades_meta/scaffolds.fasta -t 100 -h >> stats_spades.txt
# Rename SPAdes output fasta files
mv spades_standard/scaffolds.fasta spades_standard/spades_std.fasta
mv spades_meta/scaffolds.fasta spades_meta/spades_meta.fasta
##
# Compare the two assemblies using MetaQUAST
metaquast -o metaquast/ --no-plots -t 4 \
spades_standard/spades_std.fasta spades_meta/spades_meta.fasta
# Genome content is identified from the SILVA rRNA database using BLASTN by default
# Outputs 30 reference genomes with the best scores
# Check metaquast output report and taxonomy chart
firefox metaquast/report.html metaquast/krona_charts/summary_taxonomy_chart.html &
##
# Make database of 120 different bacterial genomes (makeblastdb)
# Information is added to each sequence such as taxon id
makeblastdb -taxid_map taxid_map -parse_seqids -in genomes.fasta -dbtype nucl
# blastn of assembled sequences against the database
blastn -query spades_meta/spades_meta.fasta \
-db genomes.fasta -outfmt '6 qseqid staxids bitscore std' \
-num_threads 2 > blast.out
##
# Visualize taxonomic contents of metagenome assembly
# Create an index of the SPAdes assembly 
bwa index -a is spades_meta/spades_meta.fasta
# Find the suffix array coordinates of the single-end reads (.sai files)
mkdir aligned_reads_meta
bwa aln -t 4 spades_meta/spades_meta.fasta trimmed_reads_skewer/*pair1.fastq \
> aligned_reads_meta/ERR260505_1.sai
bwa aln -t 4 spades_meta/spades_meta.fasta trimmed_reads_skewer/*pair2.fastq \
> aligned_reads_meta/ERR260505_2.sai
# Generate paired-end alignments in SAM format
bwa sampe spades_meta/spades_meta.fasta \
aligned_reads_meta/ERR260505_1.sai aligned_reads_meta/ERR260505_2.sai \
trimmed_reads_skewer/ERR260505-trimmed-pair1.fastq trimmed_reads_skewer/ERR260505-trimmed-pair2.fastq \
> aligned_reads_meta/ERR260505.sam
# Convert sam to bam, then sort and index
samtools view -S -b aligned_reads_meta/ERR260505.sam > aligned_reads_meta/ERR260505.bam
samtools sort aligned_reads_meta/ERR260505.bam -o aligned_reads_meta/ERR260505.sorted.bam
samtools index aligned_reads_meta/ERR260505.sorted.bam
# Collate data in the three files to make a blobDB.json file using blobtools
blobtools create -i spades_meta/spades_meta.fasta \
-b aligned_reads_meta/ERR260505.sorted.bam \
-t blast.out \
--db nodesDB.txt
# Draw blob plots [-m? -r?] of .json file
blobtools plot -i blobDB.json -m -r species
mkdir blobplots
mv *.png blobplots/
# Check blobtools outputs in order (last: read_cov file)
xdg-open <file> &
