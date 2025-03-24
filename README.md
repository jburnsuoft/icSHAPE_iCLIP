# icSHAPE_iCLIP

This repository contains tools to analyze the overlap of experimental RNA-Binding Protein (RBP) binding coordinates (iCLIP/eCLIP) with in vivo and in vitro icSHAPE RNA structural profiles. These scripts facilitate the integration of RNA structural data with RBP binding data to uncover insights into post-transcriptional regulatory processes.

## Publication Reference

This repository supports the analysis presented in the following publication:

*C2H2-zinc-finger transcription factors bind RNA and function in diverse post-transcriptional regulatory processes*,  
**Molecular Cell**, Volume 84, Issue 19, 2024, Pages 3810-3825.e10,  
ISSN 1097-2765, [https://doi.org/10.1016/j.molcel.2024.08.037](https://doi.org/10.1016/j.molcel.2024.08.037).

Program conception, final development, and troubleshooting were performed by **James D. Burns** (Greenblatt and Zhang Labs, UofT). Initial development was performed in collaboration with student mentee **Ryan Denniston** (Zhang Lab, UofT).

---

## Overview of Scripts

### 1. `CLIPSHAPE.R`
**Purpose:**  
This script identifies overlapping regions between iCLIP/eCLIP peak coordinates and icSHAPE RNA structural data. It extracts sequences and icSHAPE scores for regions of interest.

**Inputs:**  
- `clipfile`: iCLIP/eCLIP peak coordinate file in `.bed` format.  
- `plusShapeName` and `minusShapeName`: `.bed` files containing icSHAPE scores for the positive and negative strands, respectively.  
- `shapeRange`: Range around iCLIP peaks to search for overlaps.  
- `clusterDistance`: Maximum distance between peaks to consider them part of the same cluster.  
- Reference genome file (`GRCh37.p13.genome.fa`).  

**Outputs:**  
- Fasta file containing sequences around iCLIP peaks that overlap with icSHAPE data.  
- `.SHAPE` files containing icSHAPE scores for overlapping regions.  

---

### 2. `MAPSHAPE.R`
**Purpose:**  
This script averages and visualizes the overlap between icSHAPE and iCLIP coordinates. It generates median icSHAPE profiles for in vivo, in vitro, and shuffled datasets.

**Inputs:**  
- `vitroShapeFolder`: Folder containing in vitro icSHAPE `.SHAPE` files.  
- `vivoShapeFolder`: Folder containing in vivo icSHAPE `.SHAPE` files.  
- `plFolder`: Folder containing RNAplFold unpaired probability `.lunp` files.  
- `plshuffle`: Folder containing shuffled RNAplFold `.lunp` files.  

**Outputs:**  
- Median icSHAPE profile plots saved as a PDF (`<protein_name>_median_plot.pdf`).  
![image](https://github.com/user-attachments/assets/56203b12-4eb4-4e80-a6a5-b19da562ded7)


---

### 3. `icSHAPE_HEATMAP.R`
**Purpose:**  
This script generates heatmaps and performs clustering analysis on the differences between in vivo and in vitro icSHAPE profiles. It identifies optimal clustering parameters using various methods.

**Inputs:**  
- `vitroShapeFolder`: Folder containing in vitro icSHAPE `.SHAPE` files.  
- `vivoShapeFolder`: Folder containing in vivo icSHAPE `.SHAPE` files.  

**Outputs:**  
- Heatmaps of icSHAPE differences between in vivo and in vitro datasets.  
- Clustering analysis results using methods such as silhouette analysis, gap statistics, hierarchical clustering, and Bayesian maximization.  

---

## How to Use

1. **Prepare Input Files:**  
   - Ensure that iCLIP/eCLIP peak coordinates are in `.bed` format.  
   - Prepare icSHAPE data in `.bed` or `.SHAPE` formats as required.  
   - Place input files in appropriate folders.

2. **Run Scripts:**  
   - Use `CLIPSHAPE.R` to extract overlapping regions and generate `.SHAPE` files.  
   - Use `MAPSHAPE.R` to visualize median icSHAPE profiles.  
   - Use `icSHAPE_HEATMAP.R` to perform clustering analysis and generate heatmaps.

3. **Output Files:**  
   - Outputs will be saved in the specified directories or as PDF files.

---

## Dependencies

The following R libraries are required to run the scripts:
- `seqinr`
- `Biostrings`
- `pheatmap`
- `factoextra`
- `NbClust`
- `fpc`
- `cluster`
- `vegan`
- `mclust`
- `apcluster`

Install missing libraries using `install.packages("<package_name>")` or `BiocManager::install("<package_name>")` for Bioconductor packages.

---

## Example Workflow

1. Run `CLIPSHAPE.R` to extract overlapping regions:
   ```r
   CLIPSHAPE("example_clip.bed", 25, 0, "example_plus.bed", "example_minus.bed", 
             "output_sequences.fa", "output_seq_folder/", 
             "output_shape_folder/", "output_ct_folder/", "output_db_folder/")
   
