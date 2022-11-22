# Adrenal Gland ENCODE single-cell analysis
Contact information: jschaepe@stanford.edu

Integrated analysis for the following ENCODE 10x multiome adrenal gland experiments:
- ENCSR194KHA (scATAC)
- ENCSR420EWQ (scATAC)
- ENCSR693GAD (scATAC)
- ENCSR726IPC (snRNA)
- ENCSR362YDM (snRNA)
- ENCSR724KET (snRNA)

The samples are referred to in the analysis as:
- W62 (ENCSR693GAD, ENCSR724KET)
- W64 (ENCSR420EWQ, ENCSR362YDM)
- W73 (ENCSR194KHA, ENCSR726IPC)

Analysis steps in the code:
1. Loading scATAC data for each sample and creating arrow files with ArchR.
2. Filtering scATAC data with a minTSS of 5 and minFrags of 3,162 and outputting cell barcodes that pass.
3. Loading snRNA data with Seurat. 
4. Filtering snRNA data with number RNA reads between 1,000 and 30,000 and percent of mitochondrial reads less than 10%.
5. Additional filtering of snRNA for only cells that pass scATAC QC.
6. Integrating snRNA datasets with scTransform and performing dimensionality reduction and clustering.
7. Visualizing snRNA data.

Required packages:
- library(ArchR), version 1.0.1
- library(Seurat), version 4.0.6
- library(dplyr)
- library(patchwork)
- library(ggplot2)
- library(sctransform)

Installation instructions:
1. Install ArchR following instructions here: https://www.archrproject.com/
2. Install Seurat following instructions here: https://satijalab.org/seurat/articles/install.html
