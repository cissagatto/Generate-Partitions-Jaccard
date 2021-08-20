cat("\n\n################################################################################################")
cat("\n# START JACCARD PARTITIONS                                                                       #")
cat("\n##################################################################################################\n\n") 

##################################################################################################
# Generate Jaccard Partitions                                                                    #
# Copyright (C) 2021                                                                             #
#                                                                                                #
# This program is free software: you can redistribute it and/or modify it under the terms of the #
# GNU General Public License as published by the Free Software Foundation, either version 3 of   #  
# the License, or (at your option) any later version. This program is distributed in the hope    #
# that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of         #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    #
# more details.                                                                                  #     
#                                                                                                #
# Elaine Cecilia Gatto | Prof. Dr. Ricardo Cerri | Prof. Dr. Mauri Ferrandin                     #
# Federal University of Sao Carlos (UFSCar: https://www2.ufscar.br/) Campus Sao Carlos           #
# Computer Department (DC: https://site.dc.ufscar.br/)                                           #
# Program of Post Graduation in Computer Science (PPG-CC: http://ppgcc.dc.ufscar.br/)            #
# Bioinformatics and Machine Learning Group (BIOMAL: http://www.biomal.ufscar.br/)               #
#                                                                                                #
##################################################################################################

##################################################################################################
# Script
##################################################################################################

#################################################################################################
# Configures the workspace according to the operating system                                     #
##################################################################################################
sistema = c(Sys.info())
FolderRoot = ""
if (sistema[1] == "Linux"){
  FolderRoot = paste("/home/", sistema[7], "/Generate-Partitions-Jaccard", sep="")
} else {
  FolderRoot = paste("C:/Users/", sistema[7], "/Generate-Partitions-Jaccard", sep="")
}
FolderScripts = paste(FolderRoot, "/scripts", sep="")


##################################################################################################
# Options Configuration                                                                          #
##################################################################################################
options(java.parameters = "-Xmx32g")
options(show.error.messages = TRUE)
options(scipen=30)



##################################################################################################
# Read the dataset file with the information for each dataset                                    #
##################################################################################################
setwd(FolderRoot)
datasets <- data.frame(read.csv("datasets.csv"))



##################################################################################################
# ARGS COMMAND LINE                                                                              #
##################################################################################################
cat("\nGet Args")
args <- commandArgs(TRUE)



##################################################################################################
# Get dataset information                                                                        #
##################################################################################################
ds <- datasets[args[1],]


##################################################################################################
# Get dataset information                                                                        #
##################################################################################################
number_dataset <- as.numeric(args[1])
cat("\nGPJ \t number_dataset: ", number_dataset)



##################################################################################################
# Get the number of cores                                                                        #
##################################################################################################
number_cores <- as.numeric(args[2])
cat("\nGPJ \t cores: ", number_cores)



##################################################################################################
# Get the number of folds                                                                        #
##################################################################################################
number_folds <- as.numeric(args[3])
cat("\nGPJ \t folds: ", number_folds)



##################################################################################################
# Get the number of folds                                                                        #
##################################################################################################
folderResults <- toString(args[4])
cat("\nGPJ \t  folder: ", folderResults)



##################################################################################################
# Get dataset name                                                                               #
##################################################################################################
dataset_name <- toString(ds$Name) 
cat("\nGPJ \t nome: ", dataset_name)



##################################################################################################
# DON'T RUN -- it's only for test the code
# ds <- datasets[29,]
# dataset_name = ds$Name
# number_dataset = ds$Id
# number_cores = 10
# number_folds = 10
# folderResults = "/dev/shm/res"
##################################################################################################



##################################################################################################
# cat("\n\nCopy FROM google drive \n")
# destino = paste(FolderRoot, "/datasets/", dataset_name, sep="")
# origem = paste("cloud:elaine/Datasets/CrossValidation_WithValidation/", dataset_name, sep="")
# comando = paste("rclone -v copy ", origem, " ", destino, sep="")
# cat("\n", comando, "\n") 
# a = print(system(comando))
# a = as.numeric(a)
# if(a != 0) {
# stop("Erro RCLONE")
# quit("yes")
# }



##################################################################################################
cat("\nCreate Folder")
if(dir.exists(folderResults)==FALSE){
  dir.create(folderResults)
}


##################################################################################################
# LOAD RUN.R                                                                                     #
##################################################################################################
setwd(FolderScripts)
source("run.R") 


##################################################################################################
# GET THE DIRECTORIES                                                                            #
##################################################################################################
cat("\nGet directories\n")
diretorios = directories(dataset_name, folderResults)


##################################################################################################
# execute the code and get the total execution time                                              #
# n_dataset, number_cores, number_folds                                                          #
##################################################################################################
cat("\nExecute Jaccard Partitions \n")
timeFinal <- system.time(results <- executeGPJ(number_dataset, number_cores, number_folds, folderResults))
print(timeFinal)


Folder = paste(diretorios$folderResults, "/", dataset_name, sep="")

##################################################################################################
# save the total time in rds format in the dataset folder                                        #
##################################################################################################
cat("\nSave Rds\n")
str0 <- paste(Folder, "/", dataset_name, "-jaccard-results.rds", sep="")
save(results, file = str0)


##################################################################################################
# save results in RDATA form in the dataset folder                                               #
##################################################################################################
cat("\nSave Rdata \n")
str1 <- paste(Folder, "/", dataset_name, "-jaccard-results.RData", sep="")
save(results, file = str1)


##################################################################################################
cat("\nCompress results")
setwd(Folder)
str3 = paste("tar -zcvf ", dataset_name, "-results-gpj.tar.gz ", Folder, sep="")
print(system(str3))


########################################################################################################################
# cat("\n Copy to OUTPUT google drive")
# origem = diretorios$folderOutputDataset
# destino = paste("cloud:elaine/[2021]ResultadosExperimentos/Generate-Partitions-Jaccard/", dataset_name, sep="")
# comando = paste("rclone -v copy ", origem, " ", destino, sep="")
# print(system(comando))
# cat("\n", comando, "\n") 
# a = print(system(comando))
# a = as.numeric(a)
# if(a != 0) {
# stop("Erro RCLONE")
# quit("yes")
# }


########################################################################################################################
# cat("\n Copy to RESULTS google drive")
# origem = paste(Folder, "/", dataset_name, "-results-gpj.tar.gz ", sep="")
# destino = paste("cloud:elaine/[2021]ResultadosExperimentos/Generate-Partitions-Jaccard/", dataset_name, sep="")
# comando = paste("rclone -v copy ", origem, " ", destino, sep="")
# cat("\n", comando, "\n") 
# a = print(system(comando))
# a = as.numeric(a)
# if(a != 0) {
# stop("Erro RCLONE")
# quit("yes")
# }


##################################################################################################
#cat("\nDelete folder output dataset \n")
#str6 = paste("rm -r ", diretorios$folderOutputDataset, sep="")
#print(system(str6))


##################################################################################################
#cat("\nDelete folder specific dataset \n")
#str7 = paste("rm -r ", diretorios$folderSpecificDataset, sep="")
#print(system(str7))


rm(list = ls())

gc()

cat("\n##################################################################################################")
cat("\n# END OF RANDOM PARTITIONS VERSION 1. Thanks God!                                                #") 
cat("\n##################################################################################################")
cat("\n\n\n\n") 

##################################################################################################
# Please, any errors, contact us: elainececiliagatto@gmail.com                                   #
# Thank you very much!                                                                           #
##################################################################################################
