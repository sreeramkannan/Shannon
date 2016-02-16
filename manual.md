---
layout: page
title: "Manual"
description: ""
group: navigation
---
{% include JB/setup %}

~~~
Shannon Version 0.0.0

~~~
### Usage

For single-ended reads,

Usage: python shannon.py -o running_directory \-\-single read_file  [options]

For paired-end reads,

Usage: python shannon.py -o running_directory \-\-left read_pair1 \-\-right read_pair2  [options]

The running_directory mentions the name of a directory where Shannon can run. This directory should be empty or non-existent while starting the run. 

The reads should be in fasta or fastq format. 

~~~

### Options

The string [options] can be either empty or one or more of the following: 

-p nJobs

This option is used in order to specify the number of parallel jobs. Needs GNU parallel installed. 

-K kmerSize

This option is used to set the Kmer size.

--partition partitionSize
This option is used to set the maximum size of each partition.

--compare reference_file_name

This option is used to compare the produced output to the reference and create a log. To run this option blat has to be installed. 

~~~

### Output 

The main output is in running_directory/shannon.fasta which contains the list of reconstructed transcripts in fasta format.

There is a log file in running_directory/log.txt 

The output of the --compare option is in running_directory/compare_log.txt

The directory running_directory/TEMP contains intermediate running files and can be deleted after the run.

~~~
### Requirements

Memory: Please reserve atleast 1GB / 1 million single-end reads.

Cores: The program is partially multi-threaded and needs GNU parallel to run.

Disk Space: The program will use upto 5 times the amount of space required for storing the reads (in FASTA format). Please ensure you have this amount of space before running.
