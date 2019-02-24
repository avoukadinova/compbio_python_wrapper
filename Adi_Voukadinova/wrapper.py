import os
from Bio import SeqIO

os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_genomic.fna.gz")
os.system("gunzip GCF_000387825.2_ASM38782v2_genomic.fna.gz")
os.system("mv GCF_000387825.2_ASM38782v2_genomic.fna HM27.fna")

os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/845/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_genomic.fna.gz")
os.system("gunzip GCF_000387845.2_ASM38784v2_genomic.fna.gz")
os.system("mv GCF_000387845.2_ASM38784v2_genomic.fna HM46.fna")

os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/785/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_genomic.fna.gz")
os.system("gunzip GCF_000387785.2_ASM38778v2_genomic.fna.gz")
os.system("mv GCF_000387785.2_ASM38778v2_genomic.fna HM65.fna") 

os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/865/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_genomic.fna.gz")
os.system("gunzip GCF_000387865.2_ASM38786v2_genomic.fna.gz")
os.system("mv GCF_000387865.2_ASM38786v2_genomic.fna HM69.fna")

out_file = open("UPEC.log", "w")
gene_names = ["HM27", "HM46","HM65", "HM69"]
name_index = 0

for name in gene_names:
    total_bp = 0
    contig_count = 0
    
    for record in SeqIO.parse(str(name) + ".fna", "fasta"):
        contig_count += 1
        total_bp = total_bp + len(record.seq)
    out_file.write("The " + str(name) + " gene has " + str(contig_count) + " contigs.\n")
    out_file.write("The " + str(name) + " gene has " + str(total_bp) + " base pairs in its assembly.\n\n")

cwd = os.getcwd()
CDS_ncbi = [4984, 4783, 5018, 5299]
tRNA_ncbi = [74, 78, 77, 73]
CDS_prokka = []
tRNA_prokka = []
index = 0

for name in gene_names:
    os.system("prokka --prefix " + str(name) + " --usegenus --genus Escherichia " + str(name) + ".fna")
    out_file.write("Prokka command used for "+str(name)+" is: prokka --prefix "+str(name)+" --usegenus --genus Escherichia "+str(name)+".fna\n")
    os.chdir(str(name))

    raw_data = open(str(name)+".txt").read()
    prokka_data = list(raw_data.split("\n"))

    for line in prokka_data:
        if line[0:4] == "CDS:":
            new_line = line.split(" ")
            CDS_prokka.append(int(new_line[1]))
        elif line[0:5] == "tRNA:":
            new_line = line.split(" ")
            tRNA_prokka.append(int(new_line[1]))
        else:
            continue

    os.chdir(cwd)

    out_file.write("Prokka annotation of " + str(name)+ " is: \n")
    out_file.write(raw_data)
    out_file.write("\n")
    
    out_file.write("Prokka found "+str(CDS_ncbi[index]-CDS_prokka[index])+" less CDS and "+str(tRNA_prokka[index]-tRNA_ncbi[index])+" more tRNA than the RefSeq in assembly "+str(name)+"\n\n\n")

    index += 1


#os.system("wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR127/SRR1278956/SRR1278956.sra")
#os.system("fastq-dump -I --split-files SRR1278956.sra")

#os.system("wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR127/SRR1278960/SRR1278960.sra")
#os.system("fastq-dump -I --split-files SRR1278960.sra")

#os.system("wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR128/SRR1283106/SRR1283106.sra")
#os.system("fastq-dump -I --split-files SRR1283106.sra")

#os.system("wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR127/SRR1278963/SRR1278963.sra")
#os.system("fastq-dump -I --split-files SRR1278963.sra")

#os.system("cp HM27/HM27.gff .")
#os.system("cp HM46/HM46.gff .")
#os.system("cp HM65/HM65.gff .")
#os.system("cp HM69/HM69.gff .")

#os.system("mv SRR1278956_1.fastq HM27_1.fastq")
#os.system("mv SRR1278956_2.fastq HM27_2.fastq")
#os.system("bowtie2-build HM27.fna HM27bow")
#os.system("tophat -G HM27.gff --transcriptome-index=HM27bow HM27bow")
#os.system("tophat2 -p 4 -o HM27out HM27bow HM27_1.fastq HM27_2.fastq")

#os.system("mv SRR1278960_1.fastq HM46_1.fastq")
#os.system("mv SRR1278960_2.fastq HM46_2.fastq")
#os.system("bowtie2-build HM46.fna HM46bow")
#os.system("tophat -G HM46.gff --transcriptome-index=HM46bow HM46bow")
#os.system("tophat2 -p 4 -o HM46out HM46bow HM46_1.fastq HM46_2.fastq")

#os.system("mv SRR1283106_1.fastq HM65_1.fastq")
#os.system("mv SRR1283106_2.fastq HM65_2.fastq")
#os.system("bowtie2-build HM65.fna HM65bow")
#os.system("tophat -G HM65.gff --transcriptome-index=HM65bow HM65bow")
#os.system("tophat2 -p 4 -o HM65out HM65bow HM65_1.fastq HM65_2.fastq")

#os.system("mv SRR1278963_1.fastq HM69_1.fastq")
#os.system("mv SRR1278963_2.fastq HM69_2.fastq")
#os.system("bowtie2-build HM69.fna HM69bow")
#os.system("tophat -G HM69.gff --transcriptome-index=HM69bow HM69bow")
#os.system("tophat2 -p 4 -o HM69out HM69bow HM69_1.fastq HM69_2.fastq")

#os.system("cufflinks -p 4 -G HM27.gff -o HM27 HM27out/accepted_hits.bam")
#os.system("cufflinks -p 4 -G HM46.gff -o HM46 HM46out/accepted_hits.bam")
#os.system("cufflinks -p 4 -G HM65.gff -o HM65 HM65out/accepted_hits.bam")
#os.system("cufflinks -p 4 -G HM69.gff -o HM69 HM69out/accepted_hits.bam")

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
