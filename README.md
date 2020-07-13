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
In this example we will recreate the following figure:

`![alt text](https://github.com/patfiaux/GenomeColoR/Figures/GATA1_GenomeColoR_plot.png?raw=true)`

The figure shows the results of analyzing data from a CRISPR inhibition (CRISPRi) screen by [Fulco et al. (2016)](https://science.sciencemag.org/content/354/6313/769.abstract) with [RELICS](https://github.com/patfiaux/RELICS), a tool developped specifically for analyzing tiling CRISPR screens for the detection of functional sequences (FS). RELICS reports all functional sequences detected and their location. To keep apart the different FS we will give each of them a distinct color.

We will add in the log2 fold change of the raw data as a separate track as well as an epigentic track that is commonly used as a heuristic for detecting functional sequences.

We recommend that you navigate to the GenomeColoR directory. In an interactive R session

```
setwd('path/to/GenomeColoR')
```

## 1. Source the script
This will load in all the libraries and necessary functions

```
source('path/to/Code/GenomeColoR.R')

# if you are already in the GenomeColoR directory:
# source('Code/GenomeColoR.R')
```

## 2. Load in the data
We load in the RELICS scores as well as the raw data given as log2 fold change.
```
gata1.scores <- read.csv('Data/GATA1_FS_track.csv', stringsAsFactors = F)
gata1.l2fc <- read.csv('Data/GATA1_l2fc_track.csv', stringsAsFactors = F)
```

Finally, we also load in H3K27ac data, an epigentic marker for open chromatin.
```
gata1.h3k27ac <- read.table('Data/GATA1_H3K27ac_track.txt', header = T, sep = '\t', stringsAsFactors = F)
```

## 3. Set up the input
The main function `plot_tracks` requires at least 3 inputs. 

### 3.1 Track colors 
Specify the colors for each of the tracks from 3.1. Keep in mind that the colors are assigned to the track in the order they are given, not by the name in the list. The latter just helps keeping track of what tracks already have been colored.

The colors for each track are given as `data.frame` where each row corresponds to a label present in the track to plot. 
```
relics.colors <- data.frame(label = c("FS 0", "FS 1", "FS 2", "FS 3", "FS 4", "FS 5", "non_Functional"),
                            color = c("#FF0000", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#000000"))
                            
l2fc.colors <- data.frame(label = c("l2fc"),
                          color = c('grey30'), 
                          stringsAsFactors = F)
                          
h3k27.colors <- data.frame(label = c("H3K27ac"), 
                           color = c('darkred'), 
                           stringsAsFactors = F)
                           
track.colors <- list(RELICS = relics.colors,
                     L2FC = l2fc.colors,
                     H3K27ac = h3k27.colors)

```

Note that because the labels are treated as factors, it's important to fix their ordering according to the order of the color data frame, else color will not match its label. This is only necessary if there are multiple colors / labels within a track.

```
gata1.scores$label <- factor(gata1.scores$label, levels = relics.colors$label)
```

### 3.2 Track score type
For each track, provide the type of track you want to plot.
In this case, both the RELICS scores and the H3K27ac track are plotted as continous tracks coloring the height from 0 to the given score. However, the log2 fold change track is plotted with individual points.

```
track.types <- c('genomeScore', 'pointCloud', 'genomeScore')
```

### 3.3 Tracks to plot
The tracks to plot are provided as a list. Make sure you place the tracks into the list after dealing with the factor issue mentioned in 3.1.

```
track.scores <- list(RELICS = gata1.scores,
                    L2FC = gata1.l2fc,
                    H3K27ac = gata1.h3k27ac)
```

# 4 Track height
Differnt tracks might have varying importance in your plot. Incrrease the size of on track relative to another to draw attention to the key elements of your figure.

```
# the track heights are relative, so the numbers here are selected somewhat arbitrarily.
track.heights <- c(4,4,0.5)
```

## 5 Plotting
If nothing is specified, GenomeColoR will use hg37 as reference genome. 

```
# This can take a few seconds to load
plot_tracks(track.scores,
            track.types,
            track.colors,
            track.height = track.heights)
```

# Advanced flags
