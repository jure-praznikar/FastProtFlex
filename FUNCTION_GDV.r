# This work is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0). 
# You may copy, redistribute, remix, and build upon the material in any medium, provided that you do not use it for 
# commercial purposes. Proper attribution must be given to the original creator.
# For full license details, please visit:
# https://creativecommons.org/licenses/by-nc/4.0/

# ALL AT ONCE ####################
GDVallatonce <- function(pdb) { 
  xyz<-t(matrix(pdb$xyz,nrow=3)) #vector to Nx3 matrix  
  DDsave<-dist(xyz)
  DD<-DDsave*0
  threshold<-7 #treshold=7.0Å to create graph
  DD[DDsave<threshold]<-1
  adj<-as.matrix(DD)
  #g<-graph.adjacency(adj,mode="undirected",weighted=NULL)
  g<-graph_from_adjacency_matrix(adj,mode="undirected",weighted=NULL)
  g<-delete_vertices(g,which(degree(g)==0)) #delete isolated vertices
  rm(adj,DD)
  go<-count_orbits_per_node(g,max_graphlet_size=4) # orbit
  #smooth
  DD<-DDsave*0
  DD[DDsave<=2.1]<-1 # covalent bond cutoff (2.1A)for smoothing
  adj<-as.matrix(DD)
  goSM<-go*0 # for smoothed
  for (c in 1:dim(adj)[1]){
       n<-which(adj[c,]==1)
       if (length(n)>=2) {goSM[c,]<-(go[c,] + colMeans(go[n,]))*0.5}
            else if (length(n)==1) {goSM[c,]<-(go[c,] + go[n,])*0.5}
            else if (length(n)==0) {goSM[c,]<-go[c,]}
  }
  # end smooth
  return(goSM)
}

###############################################################
###############################################################
###############################################################

# BY-PARTS #################### split to chunks of "equal" sizes
# select residues around the chain, cutoff is 12-15 A
GDVperpartes_new2 <- function(pdb) {
  #pdb$atom$chain<-"A" 
  pdb<-suppressWarnings(clean.pdb(pdb, consecutive = TRUE, force.renumber = TRUE))
  pdb$atom$chain<-"A" 
  L<-sum(pdb$calpha) # number of residues
  # split to "target size" chunks
  target_size<-100 # the size of chunck
  cat("Size of chunk: ",target_size," residues. \n")
  k <- ceiling(L / target_size)  # number of parts
  indices <- ceiling(seq_along(1:L)/L*k)
  SPL<-split(1:L, indices)  
  aL<-which(pdb$calpha==TRUE) # all CA atoms
  cat('Total number of chunks: ',length(SPL),'\n')
  goSM<-numeric(0)
  xyz<-t(matrix(pdb$xyz,nrow=3)) #vector to Nx3 matrix
  for (r in 1:length(SPL)) { 
          mySEL<-SPL[[r]]
          CAinds<-which(pdb$calpha==TRUE)[mySEL]
          
      #cat('Chunk: ',r,' RESIDUE index (from-to): ',pdb$atom$resno[CAinds[1]],
      #                                pdb$atom$resno[CAinds[length(CAinds)]],'\n')
      
      inds1<-atom.select(pdb,resno=mySEL) # select core (chunk)
      diffRESI<-setdiff(aL,CAinds) # not in CAinds/mySEL but it is in aL<-which(pdb$calpha==TRUE)
      # distance between CAinds(core) and the rest of protein CA-atoms diffRESI
      MAT <- as.matrix(pdist( xyz[CAinds,], xyz[diffRESI,] ))
      temp1<-apply(MAT,2,min) # find min value pair-wise
      k<-which(temp1<=15) # if less than treshold=15, then select
      inds2<-atom.select(pdb,resno=pdb$atom$resno[diffRESI[k]])
      #
      my.atoms<-rbind(xyz[inds1$atom,],xyz[inds2$atom,])
      #cat('Core: ',length(inds1$atom),
      #           'Adjacent: ',length(inds2$atom),
      #           'Total: ',length(inds1$atom)+length(inds2$atom),'\n')
      DDsave<-dist(my.atoms)
      DD<-DDsave*0
      threshold<-7.0 #treshold=7.0Å to create graph
      DD[DDsave<threshold]<-1
      adj<-as.matrix(DD)
      g<-graph_from_adjacency_matrix(adj,mode="undirected",weighted=NULL)
      g<-delete_vertices(g,which(degree(g)==0)) #delete isolated vertices
      rm(adj,DD)
      goRES<-count_orbits_per_node(g,max_graphlet_size=4) # GDV
      #smooth
      DD<-DDsave*0
      DD[DDsave<=2.1]<-1 # covalent bond cutoff (2.1A)for smoothing
      adj<-as.matrix(DD)
      goRES_SM<-goRES*0 # for smoothed
      for (c in 1:length(inds1$atom)){
          n<-which(adj[c,]==1)
          if (length(n)>=2) {goRES_SM[c,]<-(goRES[c,] + colMeans(goRES[n,]))*0.5}
          else if (length(n)==1) {goRES_SM[c,]<-(goRES[c,] + goRES[n,])*0.5}
          else if (length(n)==0) {goRES_SM[c,]<-goRES[c,]}
      }
      #end smooth
      goRES_SM<-goRES_SM[1:length(inds1$atom),] # only selected atoms "by-parts"
      goSM<-rbind(goSM,goRES_SM)
      #
      p1<-r/length(SPL)*100
      cat("\r",'Progress: ',round(p1,2),' % ')
  }
  cat('\n')
  return(goSM)
}
###########################
