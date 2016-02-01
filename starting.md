---
layout: page
title: "Getting Started"
description: ""
group: navigation
---
{% include JB/setup %}

The short tutorial below explains how to run __Shannon__ using a small example distributed with the program. 

The program includes some short read files to test if the program is running properly.


~~~
We will assemble Samples/SE_read.fasta first:

python shannon.py  -o ../ShannonOutput --single Samples/SE_read.fasta

Similarly we can assemble a fastq file now:

python shannon.py -o ../ShannonOutput --single Samples/SE_reads.fastq 

We can also assemble paired-end reads as follows:

python shannon.py -o ../ShannonOutput --left Samples/PE_reads_1.fastq --right Samples/PE_reads_2.fastq

~~~

