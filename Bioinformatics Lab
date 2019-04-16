Log into the UConn Xanadu Computing Cluster 

```
ssh USERNAME@xanadu-submit-ext.cam.uchc.edu

```

Make a folder for this lab and then navigate into that folder

```
mkdir molsys
cd molsys

```

Copy fastq files into the folder you are currently in. We will copy them from Eric's folder beacause he has already changed his poermissions to allow others to access his files.

```
cp WHERE_THE_FILE_IS . 
```

The following procedure to assemby is from Eric Gordon's bioinformatics workshop (https://github.com/erg55/Simon-lab-workshop).


Deduplication

Deduplication makes assembly faster by getting rid of optical or pcr duplicates which don't contribute to coverage.

```
module load bbmap
clumpify.sh in1=S125_R1.fastq.gz in2=S125_R2.fastq.gz out1=S125_dedup_R1.fastq.gz out2=S125_dedup_R2.fastq.gz dedupe
```

Trimming
We will also need to download a file with the sequence of the Illumina adaptor we think was used into our working directory. This command tells trimmomatic to remove any sequences matching Illumina adaptors, remove low quality (< 3 quality score) trailing or leading bases, using a sliding window of 4 bases removing windows where the quality score is less than 20 on average and finally discarding any read less than 50 bp long after all trimming.

```
module load Trimmomatic
cp ../../egordon/metagenomes/files/TruSeq3-PE.fa .
trimmomatic-0.36.jar PE -phred33 S125_dedup_R1.fastq.gz S125_dedup_R2.fastq.gz S125_forward_paired.fq.gz S125_forward_unpaired.fq.gz S125_reverse_paired.fq.gz S125_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50;
```

Merging paired reads
This can only occur if the insert size was short enough that the two 150 bp reads on either side overlapped. I have been told merging about half is a good percentage.

```
module load bbmap
bbmerge.sh in1=S125_forward_paired.fq.gz in2=S125_reverse_paired.fq.gz out=S125_merged.fq.gz outu1=S125_unmergedF.fq.gz outu2=S125_unmergedR.fq.gz
```

Combine all single read files for assembly:

```
cat S125_unmergedF.fq.gz S125_unmergedR.fq.gz S125_forward_unpaired.fq.gz S125_reverse_unpaired.fq.gz > S125_allsinglereadscombined.fq.gz
```

Run SPAdes assembler

Let's first just make sure the syntax is right for this to work. Run the command below specifying only one thread under -t:

```
/home/CAM/egordon/spades/SPAdes-3.12.0-Linux/bin/spades.py -t 1 --merged S125_merged.fq.gz -s S125_allsinglereadscombined.fq.gz -o S125trimmedspades.assembly/
```

If it looks like it's running correctly, stop it with control+C and lets try and submit it as a job because it may require some time and more memory than available on the head node. Generally you shouldn't run a bunch heavy tasks on this node.

Make a script like below named spades.sh (make sure to update your email):

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
/home/CAM/egordon/spades/SPAdes-3.12.0-Linux/bin/spades.py -t 16 --merged S125_merged.fq.gz -s S125_allsinglereadscombined.fq.gz -o S125trimmedspades.assembly/
```

Submit it:

```
sbatch spades.sh 
```

QUAST
We can use a program to get some basic assembly stats. This can be useful comparing the effectiveness of various programs or parameters.

```
module load quast
quast.py contigs.fasta
```










