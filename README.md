# FastProtFlex
![alt text](https://github.com/jure-praznikar/FastProtFlex/blob/main/proteinRMSF_small.png)  <br>
Lightweight, fast, and purely coordinate-based GDV linear model predicts protein flexibility. <br>
It can be applied in near real time  (on the order of 10 seconds) even for <br>
large proteins with 20,000 atoms on a standard desktop or laptop. <br>
Zenodo https://doi.org/10.5281/zenodo.17771418

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
Save the files FUNCTION_GDV.r, predict.r, and protein.pdb into the folder, <br> 
then run the R script "predict.r". <br>
The file "protein.pdb" is used as input. <br> 
The script creates an additional file, "proteinRMSF.pdb". <br>
Normalized RMSF values are written in the B-factor column of the standard PDB file. <br>
