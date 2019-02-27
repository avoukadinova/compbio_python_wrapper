# Import libraries
import os
from Bio import SeqIO

# QUESTION 1

# We are working with 4 E. Coli strains: HM27, HM46, HM65, HM69 
# For each, we use wget to retreive the sequence from ftp NCBI, we unzip using gunzip, and we rename the *.fna file to the strain name

# HM 27
# Get the sequence using wget and ftp link
os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_genomic.fna.gz")
# Unzip the sequence using gunzip
os.system("gunzip GCF_000387825.2_ASM38782v2_genomic.fna.gz")
# Rename the *.fna file to the strain name
os.system("mv GCF_000387825.2_ASM38782v2_genomic.fna HM27.fna")

# HM 46
os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/845/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_genomic.fna.gz")
os.system("gunzip GCF_000387845.2_ASM38784v2_genomic.fna.gz")
os.system("mv GCF_000387845.2_ASM38784v2_genomic.fna HM46.fna")

# HM 65
os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/785/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_genomic.fna.gz")
os.system("gunzip GCF_000387785.2_ASM38778v2_genomic.fna.gz")
os.system("mv GCF_000387785.2_ASM38778v2_genomic.fna HM65.fna") 

# HM 69
os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/865/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_genomic.fna.gz")
os.system("gunzip GCF_000387865.2_ASM38786v2_genomic.fna.gz")
os.system("mv GCF_000387865.2_ASM38786v2_genomic.fna HM69.fna")

# QUESTION 2 AND 3

# We are going to write our output to a UPEC.log file
out_file = open("UPEC.log", "w")

# Use gene_names to keep track of the strains
gene_names = ["HM27", "HM46","HM65", "HM69"]
name_index = 0

# For each strain...
for name in gene_names:
    # Initialize the contig count and total base pairs to 0
    total_bp = 0
    contig_count = 0
    
    # For each record in the fna file
    for record in SeqIO.parse(str(name) + ".fna", "fasta"):
        # Increment the contig count and base pair count
        contig_count += 1
        total_bp = total_bp + len(record.seq)
   
    # Write the output to the log file
    out_file.write("The " + str(name) + " gene has " + str(contig_count) + " contigs.\n")
    out_file.write("The " + str(name) + " gene has " + str(total_bp) + " base pairs in its assembly.\n\n")

# QUESTION 4, 5, AND 6
# We are going to use Prokka to annotate the four strains

cwd = os.getcwd()

# Initialize the CDS and tRNA counts from NCBI
CDS_ncbi = [4984, 4783, 5018, 5299]
tRNA_ncbi = [74, 78, 77, 73]

# Initialize the CDS and tRNA counts found using Prokka
CDS_prokka = []
tRNA_prokka = []

index = 0
# Go through each strain
for name in gene_names:
    # Run Prokka
    os.system("prokka --prefix " + str(name) + " --usegenus --genus Escherichia " + str(name) + ".fna")
    # Write the command to the log file
    out_file.write("Prokka command used for "+str(name)+" is: prokka --prefix "+str(name)+" --usegenus --genus Escherichia "+str(name)+".fna\n")
    # Go into the Prokka directory
    os.chdir(str(name))
    
    # Read in the data from the Prokka output
    raw_data = open(str(name)+".txt").read()
    prokka_data = list(raw_data.split("\n"))
    
    # For each line in the Prokka output
    for line in prokka_data:
        # If the line refers to the CDS
        if line[0:4] == "CDS:":
            new_line = line.split(" ")
            #append the CDS
            CDS_prokka.append(int(new_line[1]))
        # If the line refers to the tRNA
        elif line[0:5] == "tRNA:":
            new_line = line.split(" ")
            #append the tRNA
            tRNA_prokka.append(int(new_line[1]))
        else:
            continue

    #return to the main directory
    os.chdir(cwd)

    # Write the output to the log file
    out_file.write("Prokka annotation of " + str(name)+ " is: \n")
    out_file.write(raw_data)
    out_file.write("\n")
    
    out_file.write("Prokka found "+str(CDS_ncbi[index]-CDS_prokka[index])+" less CDS and "+str(tRNA_prokka[index]-tRNA_ncbi[index])+" more tRNA than the RefSeq in assembly "+str(name)+"\n\n\n")

    index += 1

# QUESTION 7
# We will now use tophat and cufflinks to map the reads

# Get the *.sra files using the ftp NCBI links
# Use fast-dump to split the *.sra files

# HM 27
os.system("wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR127/SRR1278956/SRR1278956.sra")
os.system("fastq-dump -I --split-files SRR1278956.sra")

# HM 46
os.system("wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR127/SRR1278960/SRR1278960.sra")
os.system("fastq-dump -I --split-files SRR1278960.sra")

# HM 65
os.system("wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR128/SRR1283106/SRR1283106.sra")
os.system("fastq-dump -I --split-files SRR1283106.sra")

# HM 69
os.system("wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR127/SRR1278963/SRR1278963.sra")
os.system("fastq-dump -I --split-files SRR1278963.sra")

# Move the *.gff files into the main directory
os.system("cp HM27/HM27.gff .")
os.system("cp HM46/HM46.gff .")
os.system("cp HM65/HM65.gff .")
os.system("cp HM69/HM69.gff .")

# For each strain, rename the two *.fastq files to the strain name
# Then, build a bowtie dir with the strain name + bow
# Then, run tophat

# HM 27
os.system("mv SRR1278956_1.fastq HM27_1.fastq")
os.system("mv SRR1278956_2.fastq HM27_2.fastq")
os.system("bowtie2-build HM27.fna HM27bow")
os.system("tophat -G HM27.gff --transcriptome-index=HM27bow HM27bow")
os.system("tophat2 -p 4 -o HM27out HM27bow HM27_1.fastq HM27_2.fastq")

# HM 46
os.system("mv SRR1278960_1.fastq HM46_1.fastq")
os.system("mv SRR1278960_2.fastq HM46_2.fastq")
os.system("bowtie2-build HM46.fna HM46bow")
os.system("tophat -G HM46.gff --transcriptome-index=HM46bow HM46bow")
os.system("tophat2 -p 4 -o HM46out HM46bow HM46_1.fastq HM46_2.fastq")

# HM 65
os.system("mv SRR1283106_1.fastq HM65_1.fastq")
os.system("mv SRR1283106_2.fastq HM65_2.fastq")
os.system("bowtie2-build HM65.fna HM65bow")
os.system("tophat -G HM65.gff --transcriptome-index=HM65bow HM65bow")
os.system("tophat2 -p 4 -o HM65out HM65bow HM65_1.fastq HM65_2.fastq")

# HM 69
os.system("mv SRR1278963_1.fastq HM69_1.fastq")
os.system("mv SRR1278963_2.fastq HM69_2.fastq")
os.system("bowtie2-build HM69.fna HM69bow")
os.system("tophat -G HM69.gff --transcriptome-index=HM69bow HM69bow")
os.system("tophat2 -p 4 -o HM69out HM69bow HM69_1.fastq HM69_2.fastq")

# Run cufflinks
os.system("cufflinks -p 4 -G HM27.gff -o HM27 HM27out/accepted_hits.bam")
os.system("cufflinks -p 4 -G HM46.gff -o HM46 HM46out/accepted_hits.bam")
os.system("cufflinks -p 4 -G HM65.gff -o HM65 HM65out/accepted_hits.bam")
os.system("cufflinks -p 4 -G HM69.gff -o HM69 HM69out/accepted_hits.bam")

# QUESTION 8
# I could not complete question 8 due to the server crashing
# If I could've run my code, these are the commands I would've run

#First, copy the bam file into the current directory
#Then, rename the bam file to include the strain name
#Finally, use samtools to sort the bam files

#Commands I would've run:

#os.system("cp HM27out/accepted_hits.bam .")
#os.system("mv accepted_hits.bam accepted_hits_27.bam")
#os.system("samtools sort accepted_hits_27.bam -o accepted_hits_27.sorted.bam")

#os.system("cp HM46out/accepted_hits.bam .")
#os.system("mv accepted_hits.bam accepted_hits_46.bam")
#os.system("samtools -o accepted_hits_46.sam accepted_hits_46.bam")
#os.system("samtools sort accepted_hits_46.sam > accepted_hits_46_sorted.bam")

#os.system("cp HM65out/accepted_hits.bam .")
#os.system("mv accepted_hits.bam accepted_hits_65.bam")
#os.system("samtools -o accepted_hits_65.sam accepted_hits_65.bam")
#os.system("samtools sort accepted_hits_65.sam > accepted_hits_65_sorted.bam")

#os.system("cp HM69out/accepted_hits.bam .")
#os.system("mv accepted_hits.bam accepted_hits_69.bam")
#os.system("samtools -o accepted_hits_69.sam accepted_hits_69.bam")
#os.system("samtools sort accepted_hits_69.sam > accepted_hits_69_sorted.bam")

out_file.close()
