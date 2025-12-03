# FastProtFlex
![alt text](https://github.com/jure-praznikar/FastProtFlex/blob/main/proteinRMSF_small.png)  <br>
Lightweight, fast, and purely coordinate-based GDV linear model predicts protein flexibility. <br>
It can be applied in near real time  (on the order of 10 seconds) even for <br>
large proteins with 20,000 atoms on a standard desktop or laptop. <br>
Zenodo https://doi.org/10.5281/zenodo.17771418
bioRxiv(preprint) https://www.biorxiv.org/content/10.64898/2025.11.30.691417v1

### **Prerequisites for running R scripts**  

**install packages**  
* install.packages("igraph")  
* install.packages("bio3d")
* install.packages("pdist")
* install.packages("remotes")
* remotes::install_github("alan-turing-institute/network-comparison")

**load libraries**  
* library(igraph)  
* library(bio3d)
* library(pdist)
* library(netdist)

**Usage**
<br>
Save the files: FUNCTION_GDV.r, predict.r, and protein.pdb into the same folder. <br> 
Run the R script "predict.r". <br>
The file "protein.pdb" is used as input. <br> 
The script creates an additional file, "proteinRMSF.pdb". <br>
Normalized RMSF values are written in the B-factor column of the standard PDB file. <br>
<br>
<br>
<br>
**Licence**
This work is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0). <br>
You may copy, redistribute, remix, and build upon the material in any medium, provided that you do not use it for <br>
commercial purposes. Proper attribution must be given to the original creator. <br>
For full license details, please visit: <br>
https://creativecommons.org/licenses/by-nc/4.0/ <br>
