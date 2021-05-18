# XRR-pipeline-for-sequencing-data-aDNA
After sequencing, on phoenix (/hpcfs/groups/acad_users/NGS) create a folder named:
Date_TypeOfSeq_Owner_InfoOfTheProject_SequencingID

Example: 200212_IXT_XRR_Teopancazco_SequencingID

# Download the sequencing run

```
curl -O url_of_the_text_file_containing_downloading_links
wget -i previous_file.txt
```
Move files in a directory titled BCL and untar:
```
for i in *.tar.gz; do
       bn=`basename -s .tar.gz $i`
       mkdir -pv $bn
       tar xvzf $i -C $bn
done
```
Remove untared files (they are too heavy)

Work in your own /hpcfs/users/aXXXXXXX/

# Demultiplexing
## Make a sample sheet

The sample sheet contains information about the study, the index sequences. It should have a format like this:
```
[Header],,,,
IEMFileVersion,4,,,
Investigator Name,Xavi,,,
Experiment Name,lane5_shotgun,,,
Date,17/02/20,,,
Workflow,GenerateFASTQ,,,
,,,,
[Reads],,,,
151,,,,
151,,,,

[Settings]
CreateFastqForIndexReads,1,,,
,,,,
[Data],,,,
Sample_ID,Sample_Name,I7_Index_ID,index
A24193,LP110_1,P7_acad_8nt_0032,CCTGGTAT
```

## bcl2fastq

```
module load arch/haswell
module load bcl2fastq2/2.19.1

bcl2fastq \
--runfolder-dir /hpcfs/groups/acad_users/NGS/directory \
--output-dir /hpcfs/users/aXXXXXXX/directory \
--sample-sheet SampleSheet.csv \
--ignore-missing-positions \
--ignore-missing-controls \
--ignore-missing-filter \
--ignore-missing-bcls \
--barcode-mismatches=0,1,2
```
# Trimming 
## AdapterRemoval
Do not recommend to use. It recovers less reads than fastp. 
(Each sample is in its own folder). 
```
module load arch/haswell
module purge
module load AdapterRemoval/2.2.1-foss-2016b

for d in * ; do
	if [ -d "$d" ]; then
		(cd "$d" &&
			AdapterRemoval \
				--adapter1	AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
        			--adapter2	AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT \
				--file1 *R1_*.fastq.gz \
				--file2 *R2_*.fastq.gz \
				--basename "$d" \
				--threads 20 \
				--gzip)
	fi
done
```
## fastp
### how to use conda? 
Setup miniconda as https://github.com/ACAD-UofA/Bioinformatics-Wiki/wiki/Setup-Miniconda
To install new softwares, such as fastp or multiqc:
```
conda activate
conda create -n fastp/multiqc
conda install -c bioconda fastp / conda install -c bioconda multiqc
```
### fastp for demultiplexed data
```
conda activate fastp
for d in * ; do
	if [ -d "$d" ]; then
		(cd "$d" &&
		fastp --verbose \
			 --merge --overlap_len_require 10 \
			 --detect_adapter_for_pe  \
			 --trim_poly_x \
			 --length_required 30 \
			 --low_complexity_filter --complexity_threshold 55 \
			 --correction \
			 --overrepresentation_analysis \
			 --cut_front --cut_front_mean_quality 20 --cut_tail --cut_tail_mean_quality 20 \
			 --json ${d}_fastp.json \
			 --html ${d}_fastp.html \
			 --report_title ${d}  \
			 --thread 8 \
			 --merged_out ${d}_fastp_collapsed_filtered.fq.gz \
			 --out1 ${d}_fastp_filtered_pair1.fq.gz \
			 --out2 ${d}_fastp_filtered_pair2.fq.gz \
			 --in1 ${d}*R1_*.fastq.gz \
			 --in2 ${d}*R2_*.fastq.gz
		)
	fi
done
```
### fastp for trimmed data (adapter removal) 
```
for d in * ; do
	if [ -d "$d" ]; then
		(cd "$d" &&
		fastp --verbose \
			 --merge --overlap_len_require 10 \
			 --detect_adapter_for_pe  \
			 --trim_poly_x \
			 --length_required 30 \
			 --low_complexity_filter --complexity_threshold 55 \
			 --correction \
			 --overrepresentation_analysis \
			 --cut_front --cut_front_mean_quality 20 --cut_tail --cut_tail_mean_quality 20 \
			 --json ${d}_fastp.json \
			 --html ${d}_fastp.html \
			 --report_title ${d}  \
			 --thread 8 \
			 --merged_out ${d}_fastp_collapsed_filtered.fq.gz \
			 --out1 ${d}_fastp_filtered_pair1.fq.gz \
			 --out2 ${d}_fastp_filtered_pair2.fq.gz \
			 --in1 ${d}.pair1.truncated.gz \
			 --in2 ${d}.pair2.truncated.gz
		)
	fi
done
```
### fastp visualisation (multiQC)
Gather all the .json files of each sample and move them in the same folder and run:
```
conda activate multiqc

multiqc . 
```
# Mapping
## bwa aln
Create a folder per sample and run a job loop:
```
for d in * ; do
	if [ -d "$d" ]; then
			
		(cd "$d" && sbatch ../04_bwa_aln_aDNA_default_xrr.sh)
	fi
done
```
Job script for shotgun data and nuclear capture (human genome reference) or mitochondrial capture (MT reference):
```
#!/bin/bash
#SBATCH -p batch        	                                # partition (this is the queue your job will be added to) 
#SBATCH -N 1               	                                # number of nodes (use a single node)
#SBATCH -n 8              	                                # number of cores (sequential job => uses 1 core)
#SBATCH --time=72:00:00    	                                # time allocation, which has the format (D-HH:MM:SS), here set to 1 hour
#SBATCH --mem=16GB 

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xxxxxxxxxx@adelaide.edu.au
#SBATCH --gres=tmpfs:300G

module load arch/haswell
module load BWA/0.7.15-foss-2017a
module load SAMtools/1.8-foss-2016b 

mkdir ${TMPDIR}/ref/
cp -rv /PATH/reference.fast* ${TMPDIR}/ref/
cp -rv *.fq.gz ${TMPDIR}/
#cp -rv BWA_aln_aDNA_Def_indexs.sai ${TMPDIR}/

bwa aln  -t ${SLURM_NTASKS} -l 1024 -n 0.01 -o 2 ${TMPDIR}/ref/reference.fasta ${TMPDIR}/*.fq.gz  > ${TMPDIR}/BWA_aln_aDNA_Def_indexs.sai

cp ${TMPDIR}/BWA_aln_aDNA_Def_indexs.sai .

bwa samse ${TMPDIR}/ref/reference.fasta \
	${TMPDIR}/BWA_aln_aDNA_Def_indexs.sai \
	${TMPDIR}/*.fq.gz \
	| samtools sort -l 5 \
	-O BAM \
	-o map.bam 
```
## Mitochondrial circular mapping:
Do not recommend to use. It recovers less reads than fastp. 
Step 1:
```
module load arch/haswell
module load BWA/0.7.15-foss-2017a

bwa aln /PATH/ref.fasta *.gz  -n 0.01 -l 1024 -f BWA_aln_aDNA_circular.sai

bwa samse /PATH/ref.fasta \
	BWA_aln_aDNA_circular.sai \
	*.gz \
	-f circmapper.sam
```
Step 2:
```
module load arch/haswell
module load Java/10.0.1

java -jar /PATH/realign-1.93.5.jar -e 500 -i *.sam -r human_g1k_v37_MT.fasta
```
Step 3:
```
module load arch/haswell
module load SAMtools/1.8-foss-2016b

samtools sort -@ 4 *.bam -o output_circmapper_realigned.sorted.bam
samtools index output_circmapper_realigned.sorted.bam
```
# Number of reads 
fastq files:
```
for d in * ; do
        if [ -d "$d" ]; then
                (cd "$d" && zcat *.fq.gz | echo $((`wc -l`/4)) > "$d"_fq_reads_number.txt)
        fi
done
```
bam files:
```
for d in * ; do
        if [ -d "$d" ]; then
                (cd "$d" && samtools view -c -F 260 *.bam > "$d"_bam_reads_number.txt)
        fi
done
```
# Preseq
Use to predict and estimate the complexity of a genomic sequencing library, equivalent to predicting and estimating the number of redundant reads from a given sequencing depth and how many will be expected from additional sequencing using an initial sequencing experiment. 
The file has to be sorted. 
```
preseq c_curve -o output_c_curve_complexity_output.txt -B input.sorted.bam

preseq lc_extrap -o output_lc_extrap_future_yield.txt -B input.sorted.bam

bam2mr input.sorted.bam -o name.mr -L 10000 -R 1000000
sort -k 1,1 -k 2,2n -k 3,3n name.mr > output.mr

preseq gc_extrap -o output_gc_extrap_future_coverage.txt output.mr

preseq bound_pop -o output_bound_pop_species_richness.txt output.mr
```
# Mark Duplicates (Picard)
```
module load arch/haswell
module load Java/10.0.1

for f in *.bam; do 

java -Xmx2g -jar /PATH/picard/build/libs/picardcloud.jar MarkDuplicates I=${f} O=${f}.mdup.bam M=${f}.mdup_metrics.txt TAGGING_POLICY=All VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=queryname

done
```
# Read Sorting (Samtools sort)
```
module load arch/haswell
module load SAMtools/1.9-foss-2016b

for f in *.mdup.bam; do 

samtools sort -O BAM -o ${f}.sorted.bam -@ 2 ${f}

done
```
# Indexing sorted bam (Samtools index)
```
module load arch/haswell
module load SAMtools/1.9-foss-2016b

for f in *.sorted.bam; do 

samtools index ${f}

done
```
# Add or replace read groups (Picard)
```
module load arch/haswell
module load Java/10.0.1

for f in *.sorted.bam ; do 

java -Xmx2g -jar /PATH/picard/build/libs/picardcloud.jar AddOrReplaceReadGroups I=${f} O=${f}.rg.bam RGID=${f} RGLB=libSim RGPL=Illumina RGPU=unitSim RGSM=${f} SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT

done
```
# Indel realignment (GATK)
Step 1:
```
module load arch/haswell
Java/1.8.0_121
module load GATK/3.7-Java-1.8.0_121

for f in *.rg.bam; do 

java -Xmx2g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ref_path -o ${f}.intervals --num_threads 2 --allow_potentially_misencoded_quality_scores -I ${f} 

done
```
Step 2:
```
module load arch/haswell
Java/1.8.0_121
module load GATK/3.7-Java-1.8.0_121

for f in *.rg.bam; do 

java -Xmx2g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner --filter_bases_not_stored --allow_potentially_misencoded_quality_scores -model USE_READS -compress 4 --targetIntervals ${f}.intervals --out ${f}.indelreal.bam -R ref_path --input_file ${f} 

done
```
# Filtering/Cleanup (Samtools view)
```
module load arch/haswell
module load SAMtools/1.9-foss-2016b

for f in *indelreal.bam; do 

samtools view -b -u -q 1 \
        -F `python -c 'print(hex(0x4|0x100|0x200|0x400|0x800))'` \
        ${f} \
        | samtools sort -O bam -o ${f}.filtered.bam

done
```
It is important to index the filtered bam file. 
# Stats (Samtools idxstats)
```
for i in *.filtered.bam; do samtools idxstats $i >${i/filtered\.bam/filtered\.bam\.idxstats};done
```
If you want to extract the mitochondrial data from the shotgun data:
```
for i in *.bam; do samtools view -bh $i MT > ${I}_MT.bam; done
```
# Authenticity
## mapDamage
```
conda activate mapDamage

for f in *.bam ; do

mapDamage -d mapDamage_${f} -i ${f} -y 0.1 -r ref_path

done
```
## pmdtools
Use to extract reads with post-mortem damage.
```
module load arch/haswell
module activate pmdtools
module load SAMtools/1.9-foss-2016b

for f in *.sorted.bam; do 

samtools view -h ${f} | python /PATH/pmdtools/share/pmdtools-0.60-2/pmdtools.0.60.py --customterminus 0,-1 --header | samtools view -Sb - > ${f}.pmds3filter.bam

done

for f in *.pmds3filter.bam; do

samtools index ${f}

done

for f in *.pmds3filter.bam; do

samtools view -c -F 260 ${f} > ${f}_reads_number.txt

done
```
# Sex assign
Based on the idxstats files of the shotgun data. 
```
module load arch/haswell
module load Python/2.7.13-foss-2016b
module load numpy/1.9.2-foss-2016b-Python-2.7.13
module load scipy/0.17.1-foss-2016b-Python-2.7.13
module load matplotlib/1.5.2-foss-2016b-Python-2.7.13

python sexassign.py -w -o sexes.pdf *.idxstats > sexes.txt
```
# ANGSD
Contamination estimate. 

Step 1 (calmd):
```
module load arch/haswell
module load SAMtools/1.9-foss-2016b

for f in *.bam; do 

samtools calmd -uAr ${f} /PATH/REF.fasta > ${f}.baq.bam

done
```
Step 2:
```
module load arch/haswell
module load ANGSD/0.931-foss-2016b
module load  R/3.6.0-foss-2016b

for f in *.bam; do

angsd -i ${f} -r X:5000000-154900000 -doCounts 1 -iCounts 1 -minMapQ 30 -minQ 30 -out ${f}

Rscript /PATH/ANGSD_program/angsd/R/contamination.R mapFile="/PATH/ANGSD_program/angsd/RES/map100.chrX.gz" hapFile="/PATH/ANGSD_program/angsd/RES/hapMapCeuXlift.map.gz" countFile="${f}.icnts.gz" mc.cores=24 > ${f}.txt

done
```
