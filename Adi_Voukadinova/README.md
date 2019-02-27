# compbio_python_wrapper

PLEASE READ PRIOR TO RUNNING THE PYTHON WRAPPER!

** Instructions **
1. Download the *ENTIRE* Adi_Voukadinova directory. This becomes the main directory in the wrapper.
2. Run the wrapper from the command line using python3
3. The output is a UPEC.log file for the genome annotation using Prokka 

** Reqired software **
1. Prokka
2. Cufflinks
3. Tophat/Tophat2
4. Biopython

** Things to note **
1. The os.system commands use 4 threads. This is to prevent crashing.
2. I did not complete #8 because the server kept crashing. I got as far as sorting my bam files. 
