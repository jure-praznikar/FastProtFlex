# This work is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0). 
# You may copy, redistribute, remix, and build upon the material in any medium, provided that you do not use it for 
# commercial purposes. Proper attribution must be given to the original creator.
# For full license details, please visit:
# https://creativecommons.org/licenses/by-nc/4.0/

rm(list=ls())
#load libraries
library(igraph)
library(bio3d)
library(pdist)
library(netdist)
# INPUT file
pdbfile<-"protein.pdb"
# GDV Graphlet Degree Vector function
source("FUNCTION_GDV.r")
#
# LINEAR model, Regression coeff.
beta0 <- 0
beta <- c(-1.9063749, -0.9075780,  2.9751331,  2.9623554, -0.9232815,
          -0.1610720,  0.1171952, -0.4793347,  0.2411250,  0.3110215,
           0.4226776, -0.7755721,  0.5147505, -1.7265063, -1.4770230)
# read pdb file
pdb<-read.pdb(pdbfile)
# exclude Hydrogen and non-protein atoms 
inds<-atom.select(pdb,'noh')
pdb<-trim.pdb(pdb,inds)
inds<-atom.select(pdb,'protein')
pdb<-trim.pdb(pdb,inds)
# use function, by-parts if more than 10000 atoms, else all-at-once
if (length(pdb$atom$b)<=10000) {    
      cat('ALL at once! \n')
      goSM<-GDVallatonce(pdb)
    } else {
      cat('By-parts! \n')
      goSM<-GDVperpartes_new2(pdb) # split to chunks of  "equal" sizes (100 AAs) 
}
# goSM is Nx15 matrix, N is number of atoms, and 15 number of orbits
goSM<-scale(log(goSM+1)) # scale and log of GDV
RMSFPR <- beta0 + goSM %*% beta # predicted RMSF-values (linear model)     
#
cat('Atoms: ',length(RMSFPR),'\n')
# write predicted values in "B-factor" column of pdb file
pdb$atom$b<-RMSFPR
write.pdb(pdb,file='proteinRMSF.pdb')
cat('Normalized (predicted) RMSF data \n')
cat('is written in B-val column in the file: proteinRMSF.pdb \n')
cat('Done! \n')
