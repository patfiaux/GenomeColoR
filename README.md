# GenomeColoR
A package for coloring genome tracks by segment types

The main feature of GenomeColoR is to highlighting different regions of the same track with different colors as specified by a set of labels, doing so for multiple tracks while also displaying genes and genome coordinates.

# Installation
GenomeColoR requires R version 3.51 or higher

# Obtain source code
Clone source code to your desired location with the following command: `git clone https://github.com/patfiaux/GenomeColoR.git`. Alternatively, download the repository.

# Install requirements
ou will need the following packages to run RELICS. If you don't have them, install them using the following commands. Installations will take about 5 minutes on a standard laptop.

ggplot2

```install.packages('ggplot2')```

# Bioconductor packages

# Input data format

# Quickstart with example data
In this example we will be using data from a CRISPR inhibition (CRISPRi) screen by Fulco et al. (2016). The study densly tiled single guide RNAs across the GATA1 locus to detect functional sequences that had an inhibitory effect on GATA1. The data provided here has already been processed and subsequently analyzed by [RELICS](https://github.com/patfiaux/RELICS), a tool developped specifically for analyzing tiling CRISPR screens for the detection of functional sequences. 

We recommend that you navigate to the GenomeColoR directory. In an interactive R session

```setwd('path/to/GenomeColoR')```

## 1. Source the script
This will load in all the libraries and necessary functions

```source('path/to/Code/GenomeColoR.R')

# if you are already in the GenomeColoR directory:
# source('Code/GenomeColoR.R')```

## 2. Load in the data

We load in the RELICS scores as well as the raw data given as log2 fold change.
```gata1.scores <- read.csv('Data/GATA1_FS_track.csv', stringsAsFactors = F)
gata1.l2fc <- read.csv('Data/GATA1_l2fc_track.csv', stringsAsFactors = F)
```

Finally, we also load in H3K27ac data, an epigentic marker for open chromatin
```gata1.h3k27ac <- read.table('Data/GATA1_H3K27ac_track.txt', header = F, sep = '\t', stringsAsFactors = F)```

# Advanced flags
