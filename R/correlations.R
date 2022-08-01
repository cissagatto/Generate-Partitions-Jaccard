##############################################################################
# Generate Jaccard Partitions                                                #
# Copyright (C) 2022                                                         #
#                                                                            #
# This program is free software: you can redistribute it and/or modify it    #
# under the terms of the GNU General Public License as published by the      #
# Free Software Foundation, either version 3 of the License, or (at your     #
# option) any later version. This program is distributed in the hope that    #
# it will be useful, but WITHOUT ANY WARRANTY; without even the implied      #
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the   #
# GNU General Public License for more details.                               #     
#                                                                            #
# Elaine Cecilia Gatto | Prof. Dr. Ricardo Cerri | Prof. Dr. Mauri Ferrandin #
# Federal University of Sao Carlos (UFSCar: https://www2.ufscar.br/) Campus  #
# Sao Carlos Computer Department (DC: https://site.dc.ufscar.br/)            # 
# Program of Post Graduation in Computer Science                             #
# (PPG-CC: http://ppgcc.dc.ufscar.br/) Bioinformatics and Machine Learning   #
# Group (BIOMAL: http://www.biomal.ufscar.br/)                               #
#                                                                            #
##############################################################################


########################################################################
# WORSKSPACE
########################################################################
FolderRoot = "~/Generate-Partitions-Jaccard"
FolderScripts = "~/Generate-Partitions-Jaccard/R"


############################################################################
# FUNCTION COMPUTE JACCARD                                                 #
#   Objective:                                                             #
#      Modeling correlations with index Jaccard                            #
#   Parameters                                                             #
#       ds: specific dataset information                                   #
#       resLS: specific dataset label space                                #
#       dataset_name: dataset name. It is used to save files.              #
#       namesLabels: label names                                           #
#       FolderHC: hclust and cutree folder path                            #
#       number_folds: number of folds created                              #
#   Return                                                                 #
#       correlation matrix                                                 #
############################################################################
computeJaccard <- function(ds,
                           dataset_name,
                           number_dataset, 
                           number_cores, 
                           number_folds, 
                           folderResults,
                           resLS,
                           namesLabels){
  
  diretorios = directories(dataset_name, folderResults)
  
    s = 1
    jaccardParalel <- foreach(s = 1:number_folds) %dopar%{
      
      cat("\n\nFold: ", s)
      
      #####################################################################
      FolderRoot = "~/Generate-Partitions-Jaccard"
      FolderScripts = "~/Generate-Partitions-Jaccard/R"
      
      #####################################################################
      # LOAD LIBRARIES
      setwd(FolderScripts)
      source("libraries.R")
      
      setwd(FolderScripts)
      source("utils.R")
      
      #####################################################################
      merge_matrix <- function(matrix_correlation){
        #cat("\n\tMerge matrix")
        matrix_correlation <- round(matrix_correlation,4)
        melt_mat_cor <- melt(matrix_correlation)
        return (melt_mat_cor)
        cat("\n")
        gc()
      }
      
      #####################################################################
      get_lower_tri<-function(matrix_correlation){
        #cat("\n\tGet Lower Tri")
        matrix_correlation_1[upper.tri(matrix_correlation)] <- NA
        return(matrix_correlation_1)
        cat("\n")
        gc()
      }
      
      #####################################################################
      get_upper_tri <- function(matrix_correlation){
        #cat("\n\tGet Upper Tri")
        matrix_correlation[lower.tri(matrix_correlation)]<- NA
        return(matrix_correlation)
        cat("\n")
        gc()
      }
      
      #####################################################################
      cut_matrix <- function(matrix_correlation, measure){
        #cat("\n\tCut Matrix")
        upper_tri <- get_upper_tri(matrix_correlation)
        melt_mat_cor <- melt(upper_tri, na.rm = TRUE) 
        return(melt_mat_cor)
        cat("\n")
        gc()
      }
      
      #####################################################################
      reorder_mat_cor <- function(matrix_correlation){
        #cat("\n\tReorder Matrix")
        dd <- as.dist((1-matrix_correlation)/2)
        hc <- hclust(dd)
        print(hc)
        matrix_correlation <- matrix_correlation[hc$order, hc$order]
        return(matrix_correlation)
        cat("\n")
        gc()
      }
      
      #####################################################################
      # created folder for the split
      FolderHCES = paste(diretorios$folderPartitions, "/Split-", s, sep="")
      if(dir.exists(FolderHCES)==TRUE){
        cat("\n")
      } else{
        dir.create(FolderHCES)  
      }
      setwd(FolderHCES)
      
      #cat("\nGET LABELS SPACE\n")
      classes = resLS$Classes[s]
      classes = data.frame(classes)
      classes = t(classes)
      
      #cat("\nCOMPUTES JACCARD\n")
      matrix_correlation = distance(classes, method = "jaccard", use.row.names = TRUE)
      write.csv(matrix_correlation, "matrix_correlation.csv", row.names = FALSE)
      nomesRotulos = colnames(matrix_correlation)
      
      #cat("\nCHECK NA\n")
      matrix_correlation_na = is.na(matrix_correlation)
      matrix_correlation_na_2 = replace(x = matrix_correlation, list = is.na(matrix_correlation), values = 0)
      
      #cat("\nGET COL NAMES\n")
      rownames(matrix_correlation) <- namesLabels
      matrix_correlation_2 = as.matrix(matrix_correlation)
      
      #cat("\nREORGANIZE\n")
      matrix_correlation_order <- reorder_mat_cor(matrix_correlation_2)
      upper_tri <- get_upper_tri(matrix_correlation_order)
      
      #cat("\nMELT MATRIX\n")
      melt_mat_cor <- melt(upper_tri, na.rm = TRUE)
      write.csv(melt_mat_cor, "melt_mat_cor.csv", row.names = FALSE)
    
      gc()
    }
  
  
  gc()
  cat("\n##################################################################")
  cat("\n# END OF COMPUTE JACCARD FUNCTION                                #")
  cat("\n##################################################################")
  cat("\n\n\n\n")
}


############################################################################
# FUNCTION CUTREE HCLUST                                                   #
#   Objective                                                              #
#       Partitions the correlation matrix using a hierarchical             #
#     clustering algorithm                                                 #
#   Parameters                                                             #
#       ds: specific dataset information                                   #
#       resLS: specific dataset label space                                #
#       dataset_name: dataset name. It is used to save files.              #
#       namesLabels: label names                                           #
#       FolderHClust: hclust and cutree folder path                        #
#       number_folds: number of folds created                              #
#   Return                                                                 #
#       partitions and graphics                                            #
############################################################################
CutreeHClust <- function(ds,
                         dataset_name,
                         number_dataset, 
                         number_cores, 
                         number_folds, 
                         folderResults,
                         namesLabels,
                         resLS){
  
  diretorios = directories(dataset_name, folderResults)
  
  s = 1
  cutreeParalel <- foreach(s = 1:number_folds) %dopar%{
    
    cat("\n\nFold: ", s)
    
    #####################################################################
    FolderRoot = "~/Generate-Partitions-Jaccard"
    FolderScripts = "~/Generate-Partitions-Jaccard/R"
    
    #####################################################################
    # LOAD LIBRARIES
    setwd(FolderScripts)
    source("libraries.R")
    
    setwd(FolderScripts)
    source("utils.R")   
      
    
    #####################################################################
    num.fold = c(0)
    num.part = c(0)
    num.group = c(0)
    names.labels = c(0)
    AllPartitions = data.frame(num.fold, num.part, num.group, names.labels)
    
    #####################################################################
    group = c(0)
    label = c(0)
    allPartitions2 = data.frame(label, group)
    
    #####################################################################
    # cat("\nData frame")
    fold = c(0)
    partition = c(0)
    num.groups = c(0)
    resumePartitions = data.frame(fold, partition, num.groups)
    
      #cat("\nCreate folder\n")
      FolderHCES = paste(diretorios$folderPartitions, "/Split-", s, sep="")
      
      FolderOS = paste(diretorios$folderOutputDataset, "/Split-", s, sep="")
      if(dir.exists(FolderOS)==FALSE){
        dir.create(FolderOS)
      }
      
      #cat("\nOpen matrix correlation\n")
      setwd(FolderHCES)
      matrix_correlation = data.frame(read.csv("matrix_correlation.csv"))
      rownames(matrix_correlation) <- namesLabels
      matrix_correlation_2 = as.matrix(matrix_correlation)    
        
      #cat("\nCreates the folder to save graphics\n")            
      FolderGraphics = paste(FolderHCES, "/Graphics", sep="")
      if(dir.exists(FolderGraphics)==TRUE){
        cat("\n")
      } else {
        dir.create(FolderGraphics)
      }
        
      #cat("\nCreates the folder to save clusters\n")                  
      FolderClusters = paste(FolderHCES, "/Clusters", sep="")
      if(dir.exists(FolderClusters)==TRUE){
        cat("\n")
      } else {
        dir.create(FolderClusters)
      }
        
        #cat("\nDEND\n")                  
      Dend <- matrix_correlation_2 %>% as.dist %>% hclust(method = "single") %>% as.dendrogram
      
      #cat("\nOTTER DENDRO\n")                  
      OtterDendro = as.dendrogram(hclust(d = as.dist(matrix_correlation_2), 
                                         method="single"))
      
      #cat("\nAsDist = as.dist(matrix_correlation)\n")                  
      AsDist = as.dist(matrix_correlation)
      
      #cat("\nAsDistMatrix = as.matrix(AsDist)\n")                  
      AsDistMatrix = as.matrix(AsDist)
      
      #cat("\nHC = hclust(AsDist, method=metodos[i])\n")                  
      HC = hclust(AsDist, method="single")
      
      #cat("\nDendro = as.dendrogram(HC)\n")                  
      Dendro = as.dendrogram(HC)
      
      #cat("\nCreates the folder to save clusters\n")                  
      DendData <- dendro_data(Dendro, type = "rectangle")
        
      #cat("\nGRAPHIC: RADIAL\n")
      setwd(FolderGraphics)
      pdf("radial.pdf", width = 10, height = 8)
      print(plot(as.phylo(HC), type = "radial", cex = 0.6, no.margin = TRUE))
      dev.off()
      cat("\n")
        
      #cat("\nGRAPHIC: FAN\n")
      pdf("fan.pdf", width = 10, height = 8)
      print(plot(as.phylo(HC), type = "fan", cex = 0.6, no.margin = TRUE))
      dev.off()
      cat("\n")
      
      #cat("\nGRAPHIC: UNROOT\n")      
      pdf("unroot.pdf", width = 10, height = 8)
      print(plot(as.phylo(HC), type = "unrooted", cex = 0.6, no.margin = TRUE))
      dev.off()
      cat("\n")
      
      #cat("\nGRAPHIC: CLADOGRAM\n")
      pdf("cladogram.pdf", width = 10, height = 8)
      print(plot(as.phylo(HC), type = "cladogram", cex = 0.6, no.margin = TRUE))
      dev.off()
      cat("\n")
        
      #cat("\nGRAPHIC: DENDRO\n")
      pdf("hc_plot.pdf", width = 10, height = 8)
      print(plot(Dendro))
      print(with(pvclust:::hc2axes(as.hclust(Dendro)), 
                 text(x.axis, y.axis, round(y.axis, 2),col = "red", adj = c(0.5, 1.5), cex = 0.5)))
      dev.off()
      cat("\n")     
        
      #####################################################################
      #cat("\nClustering: from the first to the last label\n")                  	  
      k = 1
      for(k in 1:ds$Labels){
        cat("\ncluster: ", k)   
        
        group = c(0)
        label = c("")
        clusters3 = cbind(group, label)
        
        setwd(FolderClusters)
        
        #cat("\nCUTREE\n")        
        cutLabels = cutree(HC, k)
        clusters = data.frame(cutree(HC, k))
        names(clusters) = "grupo"
        label = c(rownames(clusters))
        
        group = c(clusters$grupo)
        label = label
        clusters3 = data.frame(group, label)
        
        #cat("\nSAVE CUTREE\n")
        setwd(FolderClusters)
        write.csv(clusters3, paste("partition-", k, ".csv", sep=""), row.names = FALSE)		
        
        FolderOSP = paste(FolderOS, "/Partition-", k, sep="")
        if(dir.exists(FolderOSP)==FALSE){
          dir.create(FolderOSP)
        }
        
        setwd(FolderOSP)
        write.csv(clusters3, paste("partition-", k, ".csv", sep=""), row.names = FALSE)		
        
        cat("\nFrequencia")
        frequencia1 = count(clusters3, clusters3$group)
        names(frequencia1) = c("group", "totalLabels")
        setwd(FolderOSP)
        write.csv(frequencia1, paste("fold-", s, "-labels-per-group-partition-", k, ".csv", sep=""), row.names = FALSE)
        
        ############################################################################################################
        # cat("\nData frame")
        fold = s
        partition = k
        num.groups = k
        resumePartitions = rbind(resumePartitions, data.frame(fold, partition, num.groups))
        
        ############################################################################################################
        num.fold = s
        num.part = k
        num.group = clusters3$group
        names.labels = clusters3$label
        AllPartitions = rbind(AllPartitions, data.frame(num.fold, num.part, num.group, names.labels))
        
        ############################################################################################################
        nomesDosRotulos = clusters3$rotulos
        group = clusters3$group
        allPartitions2 = cbind(allPartitions2, group)
        b = k + 2
        names(allPartitions2)[b] = paste("partition-", k, sep="")
        
        k = k + 1 
        gc()
        
      } # fim do cluster
      
      allPartitions2 = allPartitions2[,-2]
      allPartitions2$label = nomesDosRotulos
      setwd(FolderOS)
      write.csv(allPartitions2, paste("fold-", s, "-all-partitions.csv", sep=""), row.names = FALSE)      
      
      setwd(FolderOS)
      resumePartitions2 = resumePartitions[c(-1,-2),]
      w = nrow(resumePartitions2)
      resumePartitions3 = resumePartitions2[-w,]
      write.csv(resumePartitions3, paste("fold-", s, "-groups-per-partition.csv", sep=""), row.names = FALSE)
      
      print(system(paste("rm -r ", FolderOS, "/Partition-1", sep="")))
      print(system(paste("rm -r ", FolderOS, "/Partition-", as.numeric(ds$Labels), sep="")))
        
      gc()
    } # fim do fold
  
  
  
  gc()
  cat("\n###################################################################")
  cat("\n# END OF THE FUNCTION HCLUST AND CUTREE                           #")
  cat("\n###################################################################")
  cat("\n\n\n\n")
}

#############################################################################
# Please, any errors, contact us: elainececiliagatto@gmail.com              #
# Thank you very much!                                                      #
#############################################################################