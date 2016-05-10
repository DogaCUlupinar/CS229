library(stringr)
library(MASS)

set.seed(12345)
setwd("/Users/dulupinar/Documents/UCLA/Classes/Spring16/CS229/project")

#Constants
MISSING_VAL        = -1         #the value we use to denote that it is missing
MISSING_PERCENTAGE = .2  #the percentage of missing information we want
FEATURE_COUNT      = 50        #the columns of features we would like
NUM_ITERATIONS     = 10      #changes this to something more clever
LAMBDA             = 1

if ("data_init" %in% ls() == FALSE){
  print("Reintializing Missing Data")
  #read in data
  haploid = read.csv("./data/chr-22.geno.reduced.csv",header = FALSE)
  individuals = read.table("./data/chr-22.ind")
  if("snps" %in% ls() == FALSE) snps = read.table("./data/chr-22.snp")
  
  #process data
  num_individuals = dim(haploid)[2]
  num_snps = dim(haploid)[1]
  snp_id = paste("snp",str_split_fixed(snps$V1, ":", 2)[c(1:num_snps),2],sep = "")
  diploid = haploid[,seq(1,num_individuals,2)] + haploid[,seq(2,num_individuals,2)]
  colnames(diploid) = individuals[seq(1,num_individuals,2),1]
  rownames(diploid) = snp_id

  num_individuals = num_individuals/2
  
  diploid_incomp = data.matrix(diploid)
  data_init = 1
}

#create missing values
if (abs(length(which(diploid_incomp == -1))/length(diploid_incomp) - MISSING_PERCENTAGE) > .05){
  print("Recrating missing values")
  diploid_incomp = data.matrix(diploid)
  for (i in c(1:num_snps)){
    for (j in c(1:num_individuals)){
      if(runif(1) < MISSING_PERCENTAGE){
        diploid_incomp[i,j] = MISSING_VAL
      }
    }
  }
}

##now lets do the regression
learned_individual_features = matrix(runif(num_individuals*FEATURE_COUNT),nrow = num_individuals, ncol = FEATURE_COUNT)
learned_snp_features        = matrix(runif(num_snps*FEATURE_COUNT),nrow= num_snps, ncol = FEATURE_COUNT)
#ridge regression to find movie feature vectors

if (TRUE){
  for (i in c(1:10)){
    #ridge regression to find the snp feature matrix
    for (i in c(1:num_snps)){
      y = diploid_incomp[i,diploid_incomp[i,] != MISSING_VAL]
      x = learned_individual_features[diploid_incomp[i,] != MISSING_VAL,]
      ridge = lm.ridge(y ~. + 0, data = cbind.data.frame(y,x), lambda=LAMBDA)
      learned_snp_features[i,] = coef(ridge)
    }
    
    #ridge regression to find individual feature matrix
    for (i in c(1:num_individuals)){
      y = diploid_incomp[diploid_incomp[,i] != -1,i]
      x = learned_snp_features[diploid_incomp[,i] != -1,]
      ridge = lm.ridge( y~. + 0, data = cbind.data.frame(y,x), lambda=LAMBDA)
      learned_individual_features[i,] = coef(ridge)
    }
    print(sum(abs(diploid - learned_snp_features %*% t(learned_individual_features))))
  }
}

#accuracy
my_diploid = round(learned_snp_features %*% t(learned_individual_features))
diploid_matrix = data.matrix(diploid)
real_values = diploid_matrix[which(diploid_incomp == -1)]
imputed_values = my_diploid[which(diploid_incomp == -1)]
print(length(which(real_values == imputed_values))/length(real_values))