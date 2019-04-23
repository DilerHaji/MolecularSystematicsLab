### Logging into the cluster 

Log into the UConn Xanadu Computing Cluster. The -Y flag allows you to use a graphical interface on your local computer to visualize outputs in the cluster. 

```
ssh -Y USERNAME@xanadu-submit-ext.cam.uchc.edu
```

Make a folder for this lab and then navigate into that folder

```
mkdir molsys
cd molsys
```

Copy fastq files into the folder you are currently in. We will copy them from Eric's folder beacause he has already changed his permissions to allow others to access his files. If you are in your home directory, the path to the folder containing the data in ../egordon/Dropbox/. The "." at the end of the following code tells "cp" to copy the file into the current directory that you are in. 

```
cp WHERE_THE_FILE_IS . 
```

### The following procedure to assembly is mostly from Eric Gordon's bioinformatics workshop (https://github.com/erg55/Simon-lab-workshop).


### Copying raw sequence files (.fastq.gz) to the current directory.

The sequence data is split across two lanes. Let's combine each set into the two paired ends reads before assembly:

```
cat *SAMPLE*L001_R1_001.fastq.gz *SAMPLE*L002_R1_001.fastq.gz > SAMPLE_R1.fastq.gz
cat *SAMPLE*L001_R2_001.fastq.gz *SAMPLE*L002_R2_001.fastq.gz > SAMPLE_R2.fastq.gz
```

### Deduplication
 Deduplication makes assembly faster by getting rid of optical or pcr duplicates which don't contribute to coverage.

```
module load bbmap
clumpify.sh in1=SAMPLE_R1.fastq.gz in2=SAMPLE_R2.fastq.gz out1=SAMPLE_dedup_R1.fastq.gz out2=SAMPLE_dedup_R2.fastq.gz dedupe
```

### Trimming
We will also need to download a file with the sequence of the Illumina adaptor we think was used into our working directory. This command tells trimmomatic to remove any sequences matching Illumina adaptors, remove low quality (< 3 quality score) trailing or leading bases, using a sliding window of 4 bases removing windows where the quality score is less than 20 on average and finally discarding any read less than 50 bp long after all trimming.

```
module load Trimmomatic
cp ../../egordon/metagenomes/files/TruSeq3-PE.fa .
trimmomatic-0.36.jar PE -phred33 SAMPLE_dedup_R1.fastq.gz SAMPLE_dedup_R2.fastq.gz SAMPLE_forward_paired.fq.gz SAMPLE_forward_unpaired.fq.gz SAMPLE_reverse_paired.fq.gz SAMPLE_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50;
```

### Merging paired reads
This can only occur if the insert size was short enough that the two 150 bp reads on either side overlapped. I have been told merging about half is a good percentage.

```
module load bbmap
bbmerge.sh in1=SAMPLE_forward_paired.fq.gz in2=SAMPLE_reverse_paired.fq.gz out=SAMPLE_merged.fq.gz outu1=SAMPLE_unmergedF.fq.gz outu2=SAMPLE_unmergedR.fq.gz
```

### Combine all single read files for assembly:
These are the reads that could not be merged in the previous step. Even through they couldn't be merged into a longer contiguous sequence, we want to use them in our assembly in order to get the best possible contig set that we can get given our data. 

```
cat SAMPLE_unmergedF.fq.gz SAMPLE_unmergedR.fq.gz SAMPLE_forward_unpaired.fq.gz SAMPLE_reverse_unpaired.fq.gz > SAMPLE_allsinglereadscombined.fq.gz
```

### Run SPAdes assembler

Let's first just make sure the syntax is right for this to work. Run the command below specifying only one thread under -t:

```
/home/CAM/egordon/spades/SPAdes-3.12.0-Linux/bin/spades.py -t 1 --merged SAMPLE_merged.fq.gz -s SAMPLE_allsinglereadscombined.fq.gz -o SAMPLE_trimmedspades.assembly/
```

If it looks like it's running correctly, stop it with control+C and lets try and submit it as a job because it may require some time and more memory than available on the head node. Generally you shouldn't run a bunch heavy tasks on this node.

Make a script like below named spades.sh (make sure to update your email). You can open an editable txt editor within your terminal window and copy the below text into a file named spades.sh by typing "nano spades.sh."

```
#!/bin/bash
#SBATCH --job-name=spades
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=100G
#SBATCH --mail-user=YOUR_UCONN_EMAIL_ADDRESS
#SBATCH -o myscript_%j.out
#SBATCH -e myscript_%j.err
/home/CAM/egordon/spades/SPAdes-3.12.0-Linux/bin/spades.py -t 16 --merged SAMPLE_merged.fq.gz -s SAMPLE_allsinglereadscombined.fq.gz -o SAMPLE_trimmedspades.assembly/
```

Submit it.

```
sbatch spades.sh 
```

### QUAST

We can use a program to get some basic assembly stats. This can be useful comparing the effectiveness of various programs or parameters.

```
module load quast
quast.py contigs.fasta
```

### Mapping reads to assembled contigs

Let's map our original reverse and forward reads back to the assembly to get an idea of the depth per contig (i.e., how many of the original raw reads map back to each base in each contig. 

```
module load bwa
bwa index contigs.fasta
bwa mem -t 2 -k 50 -B 10 -O 10 -T 90 contigs.fasta ../SAMPLE_dedup_R1.fastq.gz ../SAMPLE_dedup_R2.fastq.gz > bwafile
```

### Processsing mapping output 

Convert the resulting file into a BAM file and then use samtools to get the depth per contig.  

```
module load samtools
samtools view -b -F 4 bwafile > mapped.bam
samtools sort mapped.bam > mapped.sorted.bam
samtools depth mapped.sorted.bam > depth.txt
```

### Using R for further processing

Import the resulting file into R. 

```
module load R/3.4.1
R
depth <- read.table("depth.txt")
colnames(depth) <- c("contig_name", "depth", "position")
```

Use the aggregate function to calculate the per-contig mean depth and the contig length.

```
depth_depth <- aggregate(depth$depth, by = list(depth$contig_name), FUN = mean)
depth_length <- aggregate(depth$position, by = list(depth$contig_name), FUN = max) 
```

Plot a histogram of the per-contig mean depth and a plot of the relationship between per-contig mean depth and contig length. 

```
hist(depth_depth$x) 
hist(depth_depth$x, breaks = 1000, xlim = c(0, 1000))

plot(depth_depth$x, depth_length$x)
```

Pick a cutoff for filterig your contigs and save the list of contig names in a file within your current directory. 

```
contig_names <- depth_depth[depth_depth$x > 200, 1]
write.table(contig_names, "contig_names.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

### Extracting final contig set 

Exit R by typing quit(). Use SeqKit to extract all the contigs that you outputed to your contig_names.txt file. 

```
module load seqkit/0.10.0
seqkit grep --pattern-file contig_names.txt YOUR_CONTIGS.FASTA_FILE > new_contigs.fasta
```

### BLAST 
The next part assumes that you have downloaded and copied the BUSCO database into your current directory. https://busco.ezlab.org/

Make a blast database of your assembled contigs 

```
makeblastdb -in contigs.fasta -dbtype 'nucl' 
```

Use tblastn to blast each BUSCO against the assembled contigs database. The format of the output is controlled by the "-outfmt" flag. Use this reference for deciding what you want from the output – http://www.metagenomics.wiki/tools/blast/blastn-output-format-6

```
module load blast
tblastn -query ancestral -db contigs.fasta -outfmt "6 qseqid sseqid length evalue bitscore sframe sseq sstart send" -evalue 0.01 -word_size 3 -out blastout.txt
```

Enter R while in the directory where of blastout.txt. 

```
module load R/3.4.1
R
```

Load in the blast output as a table and turn it into a data.frame (This is a table where each column can be a different data class). Save that into an object called "blast". 

```
blast <- data.frame(read.table("blastout.txt"))
```

















