---
layout: page
title: "FAQ"
group: navigation
---

{% include JB/setup %}

- I'm having trouble with __Shannon__. Can I get help?
  - Yes. If you think you have discovered a bug that needs to be fixed please
    file a report on the GitHub page. If you have a question about installing
    or running the program please ask on the [shannon-users Google user
    group](https://groups.google.com/forum/#!forum/shannon-users).


- Why does Shannon need a pre-error correction step?
  - The main reason for running a pre-error correction software is *not* accuracy but computational efficiency. Any pre-error correction software that can quickly weed out most erroneous Kmers without getting rid of many true Kmers can be used along with Shannon. 

- Shannon seems to run pre-error correction sometimes but not other times. Why? 
  - If the file extension is FASTQ or FQ, then Shannon runs Quorum for error correction. If not, Shannon directly proceeds without any pre-error correction. 

- Can I use Shannon with another error correction software? 
  - Yes. While Quorum is the default eror correction software, if you have preference for another error correction software, you can first run that software and use the output FASTA file to run Shannon. In that case, Shannon will not run Quorum again. 

- Are there automatic checkpoints to restart Shannon with a prior run? 
  - At the moment, no. But Shannon is under continuous development, and we will see many features come up soon.

 

   
