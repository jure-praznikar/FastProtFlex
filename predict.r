rm(list=ls())
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
cat('Normalized (predicted) RMSF profile \n')
cat('is written in B-val column in file:proteinRMSF.pdb \n')
write.pdb(pdb,file='proteinRMSF.pdb')
