
# Nucleotide Sequence Assembly from DNA Microarray Data

## Purpose

The purpose of this project is to assemble and determine DNA sequences from 
DNA-microarray data.
The assembly of these sequences from DNA-microarray data is achieved via
Sequencing by Hybridization (SBH)

## Project Description

Short segments of nucleotides called _reads_ are provided in a fasta file.
Reads are simply varied-length segments of nucleotides which are produced
by breaking up strands of DNA. In the case of this project however, all
reads are of length 30. 

The final goal is to be able to combine and assemble these reads through 
SBH assembly to form _contigs_. A contig is a contiguous region of DNA that is formed
by combining overlapping DNA segments.
What this means is that the reads provided are going to be combined where there 
is significant overlap, until we are left with contigs.

## Instructions

- This program is implemented using `java`. As a result, it will require a
Java Runtime Environment (JRE) in order to run
- Create the fasta file that contains the reads you want to assemble and
save it in the source folder
- In your JRE, compile the program using `javac Sequencer.java`
- Run the program using `java Sequencer`
- Enter the name of the fasta file you provided to be sequenced

## Test Run

Fasta file called test.fasta is provided as follows:
```
    > read1
    GCATTTATTATCGGGCGACAATCCAATAGT
    > read2
    TCCTCTTTTAACTACAAAGCGTGTTTTTTG
    > read3
    TTAATCTGCCGTTTGTTATGGTTCTGAGAA
    > read4
    GCATTTATTATCGGGAAAGCGTGTTTTTTG
    > read5
    AAAAAAAAAAAAAAAAAAGCGTGTTTTTTG
    > read6
    AAAAAAAAAAAAAAAAAAGCGTGTTTTTTG
    > read7
    AAAGCGTGTTTTTTGAAAAAAAAAAAAAAA
    > read8
    AAAAAACTGTTGATACAGAAAAACTTTTCG
    > read9
    AAAAAACTGTTGATACAGAAAAACTTTTCG
    > read10
    GGGGAAAATTTTCCCCGGGGTTTTAAAACC
    > read11
    CGGGGTTTTAAAACCCCCCCGGGGGGGGGG
    > read12
    TTTTTTTTTTGGGGAAAATTTTCCCCGGGG
```

Run the following lines of code:

 `javac Sequencer.java`
 `java Sequencer`
 
```
    Enter File Name:
    test.fasta
```


The output is a file called resultsTest.fasta and it is located in the source folder.
Its contents are the following:
```
    > contig1
    GCATTTATTATCGGGAAAGCGTGTTTTTTGAAAAAAAAAAAAAAAAAAGCGTGTTTTTTG
    > contig2
    GCATTTATTATCGGGCGACAATCCAATAGT
    > contig3
    TCCTCTTTTAACTACAAAGCGTGTTTTTTG
    > contig4
    TTAATCTGCCGTTTGTTATGGTTCTGAGAA
    > contig5
    TTTTTTTTTTGGGGAAAATTTTCCCCGGGGTTTTAAAACCCCCCCGGGGGGGGGG
```

## Author

### Bonor Ayambem

Github - @bonor-ayambem
