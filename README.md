# SNP_calling_Arabidopsis
Pipeline for SNP calling using sequencing data from Arabidopsis thaliana


# SHORE pipeline

This pipeline is based on SHORE ([Ossowski et al., 2008](https://genome.cshlp.org/content/18/12/2024.long)). Check documentation [here](http://shore.sourceforge.net/wiki/). The pipeline described here derived from [Durvasula et al., 2017]() with minor modifications.


## Alignment


The first step is to align the data on a reference genome. In this case, we use the TAIR10 reference for *Arabidopsis thaliana* (accession Col-0). The fasta file can be downloaded on the TAIR10 [website](ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/). One should download the files for each chromosomes and concatenate them into one fasta file.

```
# Download individual fasta files
wget ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chr1.fas
wget ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chr2.fas
wget ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chr3.fas
wget ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chr4.fas
wget ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chr5.fas
wget ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chrC.fas
wget ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chrM.fas

# Concatenate into one file
cat TAIR10_chr1.fas TAIR10_chr2.fas TAIR10_chr3.fas TAIR10_chr4.fas TAIR10_chr5.fas TAIR10_chrC.fas TAIR10_chrM.fas > TAIR10.fa

```


### Generate index files for BWA

We used BWA for aligning the reads. Check pros and cons for different aligner on [SHORE wiki](http://shore.sourceforge.net/wiki/index.php/Supported_Short_Read_Aligners).

```
# Index fasta file of the reference genome (here <reference.fa> should 
# be the TAIR10.fa file generated above
bwa index -a bwtsw <reference.fa>
``` 

### Perform the alignment

We used the ALN algorithm but the MEM algorithm can be used as well and seems to be better than ALN (check this [paper](https://arxiv.org/pdf/1303.3997.pdf) from one of BWA creator Heng Li).


For single-end read data:
```
# Align reads
bwa aln -n 0.1 -o 1 <reference.fa>  <fastq_file.fastq> -f <fastq_file.fastq.sai>

# Repetitive hits will be randomly chosen
bwa samse -r "@RG\tID:$SAMPLE_NAME\tSM:$SAMPLE_NAME" <TAIR10.fa> \
	<fastq_file.sai> <fastq_file.fastq.sai> > <alignment.sam>  

```


For paired-end read data:
```
# Align the first read of the pairs
bwa aln -n 0.1 -o 1 <reference.fa>  <fastq_file_read_1.fastq> -f <fastq_file_read_1.fastq.sai> 

# Align the second read of the pairs
bwa aln -n 0.1 -o 1 <reference.fa>  <fastq_file_read_2.fastq> -f <fastq_file_read_1.fastq.sai> 

# Repetitive read pairs will be placed randomly
bwa sampe -a 500 -r "@RG\tID:$SAMPLE_NAME\tSM:$SAMPLE_NAME" <TAIR10.fa> \
	<fastq_file_read_1.fastq.sai> <fastq_file_read_2.fastq.sai> \
	 <fastq_file_read_1.fastq>  <fastq_file_read_2.fastq> > <alignment.sam>  


```

### Convert SAM file to MapList format

SHORE cannot process the SAM file derived from the BWA alignment but uses a MapList file. Check [here](http://shore.sourceforge.net/wiki/index.php/SHORE_File_Formats) for more information on SHORE formats.


```
# Convert SAM to MapList format 
shore convert Alignment2Maplist <alignment.sam> --refseq <reference.fa> > <alignment.map.list>

# Sort the MapList file
shore sort --preset maplist --infiles <file.map.list> --inplace

```


### Convert fasta reference in SHORE format

SHORE creates from the reference FASTA file mapping indices and caculates GC content and sequence complexity. Check more [here](http://shore.sourceforge.net/wiki/index.php/Shore_preprocess)

```
shore preprocess -C --indexes BWA,SuffixArray --fastafile <reference.fa> --indexfolder <output_directory>

# Check with Andrea with --indexes parameters, not in documentation http://shore.sourceforge.net/wiki/index.php/Shore_preprocess

```

The output file which should be used in the next step is named `reference.fa.shore`.


### Call SNPs

The `consensus` function identifies polymorphisms (SNPs). See documentation [here](http://shore.sourceforge.net/wiki/index.php/Shore_consensus). The two important out files generated are `quality_variant.txt` and `quality_reference.txt`.

```
shore/./shore consensus \
    -n ${numIDactual} \
    -f <reference.fa.shore>
    -o <output_directory> \
    -i <alignment.map.list> \
    -a <scoring_matrix_hom.txt> \
    -b 0.7 -g 4 -N

```
Add details about arguments here.


### Convert SHORE output to VCF

The file `quality_variant.txt` and `quality_reference.txt` should be converted into VCF file and merged.

```
shore convert Variant2VCF <reference.fa> < <quality_variant.txt> > <quality_variant.vcf>


shore convert Variant2VCF <reference.fa> < <quality_reference.txt> > <quality_reference.vcf>


```




