---
layout: page
title: "Source"
description: ""
group: navigation
---
{% include JB/setup %}

The __Shannon__ GitHub repository is [here](https://github.com/sreeramkannan/shannon). Source code can also be downloaded from the [download page](download.html). Currently, __Shannon__ can be built on Linux and Mac. 

If building on Mac, we suggest using a package manager such as [Homebrew](http://brew.sh) to download dependencies. Homebrew is easily installed by copying and pasting the command below at a terminal prompt:

`ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"`


Other dependencies are either included, or can be installed using package managers on the system.

#### Requirements: 

- A 64-bit operating system
- Python 2.7 or greater
- __Metis__ 
    - See [here](http://glaros.dtc.umn.edu/gkhome/metis/metis/download) for instructions on installing metis.
- __CVXOPT__ 
    - See [here](http://cvxopt.org/install/index.html) for installing. The following command will work:
     `pip install cvxopt` 
- __Jellyfish__ v2.0 or higher. 
    - See [here](http://www.genome.umd.edu/jellyfish.html) for installing.
- __GNU-Parallel__
	- See [here](https://www.gnu.org/software/parallel/) for installing.
- __Quorum__ 
    - See [here](http://www.genome.umd.edu/quorum.html) for installing.




#### Download

__Shannon__ is hosted on GitHub. The source code can be obtained by cloning the repository as follows:

`git clone https://github.com/sreeramkannan/Shannon.git`


#### Test



The __Shannon__ package comes with a small test transcriptome and read files that can be used to test that the package was compiled and installed correctly. The [Getting started](starting.html) page provides a complete description of how to run __Shannon__ on these files.

