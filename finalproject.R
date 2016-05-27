library(stringr)
library(MASS)
library(mixtools)

#functions
F1score = function(imputed_values,real_values){
  diff_real = real_values[which(imputed_values != real_values)]
  diff_imp  = imputed_values[which(imputed_values != real_values)]
  
  tp = length(which(imputed_values == real_values))
  fp = length(which((diff_imp - diff_real) > 0))
  fn = length(which((diff_imp - diff_real) < 0))
  precision = tp/(tp +fp)
  recall = tp/(tp+fn)
  #print(fp/length(imputed_values))
  #print(fn)
  return(2 *(precision*recall)/(precision + recall))
}

#given a k for the kfold it gives that many divisions in data and makes the nth the test data
#indexing starts at 1
kfold = function(data,k = 10,n = 1){
  data_size = length(data)
  test_size = data_size/k
  start = max((n-1)*test_size,1)
  end = min(start + test_size,data_size)
  return(data[start:end])
}

set.seed(12345)
setwd("/Users/dulupinar/Documents/UCLA/Classes/Spring16/CS229/project")

#Constants
MISSING_VAL        = -1         #the value we use to denote that it is missing
MISSING_PERCENTAGE = .1  #the percentage of missing information we want
FEATURE_COUNT      = 5        #the columns of features we would like
LAMBDA             = 3

if ("data_init" %in% ls() == FALSE){
  print("Reintializing Missing Data")
  #read in data
  haploid = read.csv("./data/chr-22.geno.reduced.csv",header = FALSE)
  diploid_benchmark = read.table("./data/testytest.txt",sep = " ")
  diploid_benchmark_ref = read.table("./data/reftest.txt",sep = " ")
  individuals = read.table("./data/chr-22.ind")
  if("snps" %in% ls() == FALSE) snps = read.table("./data/chr-22.snp")
  
  #process data
  num_individuals = dim(haploid)[2]
  num_snps = dim(haploid)[1]
  snp_id = paste("snp",str_split_fixed(snps$V1, ":", 2)[c(1:num_snps),2],sep = "")
  diploid = haploid[,seq(1,num_individuals,2)] + haploid[,seq(2,num_individuals,2)]
  colnames(diploid) = individuals[seq(1,num_individuals,2),1]
  rownames(diploid) = snp_id

  global_populations = str_split_fixed(individuals[seq(1, num_individuals, 2), 3], ":", 2)[,2]
  
  diploid_incomp = data.matrix(diploid)
  
  data_init = 1
}
diploid_matrix = data.matrix(diploid)

#create missing values
maskValues = function(diploid,MISSING_PERCENTAGE){
  #print("Recerating missing values")
  diploid_incomp = data.matrix(diploid)
  num_snps = dim(diploid)[1]
  num_individuals = dim(diploid)[2]
  for (i in c(1:num_snps)){
    for (j in c(1:num_individuals)){
      if(runif(1) < MISSING_PERCENTAGE){
        diploid_incomp[i,j] = MISSING_VAL
      }
    }
  }

  ref_missing = which(diploid_incomp == -1)
  not_missing = which(diploid_incomp != -1)

  return(list(diploid_incomp = diploid_incomp,ref_missing = ref_missing, not_missing = not_missing))
}

helper = function(x){
  if(x < 0){
    x = 0
  }
  if(x>2){
    x = 2
  }
  return(x)
}

maskValues3 = function(diploid,missing_rows,ref_end){
  #missing_rows is a vector of snps that you want to lose
  diploid_incomp = data.matrix(diploid)
  num_individuals = dim(diploid)[2]
  for (i in missing_rows){
    diploid_incomp[i,(ref_end+1):num_individuals] = -1
  }
  
  ref_missing = which(diploid_incomp == -1)
  not_missing = which(diploid_incomp != -1)
  
  return(list(diploid_incomp = diploid_incomp,ref_missing = ref_missing, not_missing = not_missing))
}

maskValues2 = function(diploid,MISSING_PERCENTAGE){
  #missing percentage of rows
  num_snps = dim(diploid)[1]
  cols_to_mask_num = max(1,MISSING_PERCENTAGE*num_snps) # number of cols to maske
  cols_to_mask_list = sample(num_snps,size = cols_to_mask_num,replace = FALSE) #the indicies of the cols that we will be masking
  diploid_incomp = data.matrix(diploid)
  for (i in cols_to_mask_list){
    diploid_incomp[i,] = 0
  }
  all = c(1:num_snps)
  return(list(diploid_incomp = diploid_incomp,ref_missing = sort(cols_to_mask_list), not_missing = all[!(all %in% cols_to_mask_list)]))
}

#if (abs(length(which(diploid_incomp == -1))/length(diploid_incomp) - MISSING_PERCENTAGE) > .05) diploid_incomp = maskValues(diploid)

ridgeImpute = function(diploid_incomp,missing,not_missing,FEATURE_COUNT,LAMBDA1,LAMBDA2,verbosity = 0){
  ##now lets do the regression
  num_individuals = dim(diploid_incomp)[2]
  num_snps = dim(diploid_incomp)[1]

  #intialize the vectors
  learned_individual_features = matrix(runif(num_individuals*FEATURE_COUNT),nrow = num_individuals, ncol = FEATURE_COUNT)
  learned_snp_features        = matrix(runif(num_snps*FEATURE_COUNT),nrow= num_snps, ncol = FEATURE_COUNT)
  
  #ridge regression to find movie feature vectors
  converge = Inf
  old_converge = -Inf
  iter = 0 
  while(converge != 0 && iter < 10){
    iter = iter + 1
    old_converge = converge
    #ridge regression to find the snp feature matrix
    for (i in c(1:num_snps)){
      y = diploid_incomp[i,diploid_incomp[i,] != MISSING_VAL]
      x = learned_individual_features[diploid_incomp[i,] != MISSING_VAL,]
      ridge = lm.ridge(y ~. + 0, data = cbind.data.frame(y,x), lambda=LAMBDA1) #should be 281
      learned_snp_features[i,] = coef(ridge)
    }
    
    #ridge regression to find individual feature matrix
    for (i in c(1:num_individuals)){
      y = diploid_incomp[diploid_incomp[,i] != -1,i]
      x = learned_snp_features[diploid_incomp[,i] != -1,]
      ridge = lm.ridge( y~. + 0, data = cbind.data.frame(y,x), lambda=LAMBDA2) #should be .24
      learned_individual_features[i,] = coef(ridge)
    }
    
    #determinig convergence
    my_diploid = round(learned_snp_features %*% t(learned_individual_features))
    real_values = diploid_incomp[not_missing]
    imputed_values = my_diploid[not_missing]
    a = (real_values - imputed_values)
    converge = sum(a^2)
    if(verbosity > 1 ) print(converge)
  }
  
  if(verbosity > 0 ) print(converge)
  return(learned_snp_features %*% t(learned_individual_features))
}

##plot stuff
#plot(as.numeric(names(ridge$GCV)),ridge$GCV,ylab = "GCV", xlab = expression(paste("Shrinkage Parameter (",lambda,")",sep = '')))

##mask whole rows
diploid_European = diploid[,which(global_populations == "EUR")] #select population
diploid_benchmark_test = cbind(diploid_benchmark[,15],diploid_benchmark)
diploid_benchmark_test = cbind(diploid_benchmark_ref,diploid_benchmark)
diploid_matrix = data.matrix(diploid_European)
total_corr = 0

incomp = maskValues(diploid_matrix,.25)
diploid_incomp = incomp$diploid_incomp

ref_missing = incomp$ref_missing
not_missing = incomp$not_missing

system.time(my_diploid_unrounded <- ridgeImpute(diploid_incomp,ref_missing,not_missing,25,2,.5,verbosity = 2))
my_diploid = round(my_diploid_unrounded)
my_diploid = apply(my_diploid,helper)
real_values = diploid_matrix[ref_missing]
imputed_values = my_diploid[ref_missing]

acc = length(which(imputed_values == real_values))/length(imputed_values)
print(sprintf("With feature count %d and LAMDA %d F! score is %f and accuracy %f",FEATURE_COUNT,LAMBDA,F1score(imputed_values,real_values),acc))
correlation = cor(imputed_values,real_values)^2
print(i)
print(correlation)
if(is.na(correlation)){
  correlation = 0
}
print("bartu")
total_corr = total_corr + correlation


##mask randomly
total_correlation = 0.861056
##gmm garbage
a = normalmixEM(my_diploid_unrounded[ref_missing],k = 3,mu = c(0,1),maxrestarts = 10)
table(apply(a$posterior,1,which.max))
table(t(diploid_matrix[ref_missing]))
table(t(my_diploid[ref_missing]))

##pretty gm
m_0_sd_1 = rnorm(300,mean = 0, sd = .2)
m_1_sd_1 = rnorm(2000,mean = 1, sd = .2)
m_2_sd_1 = rnorm(900,mean = 2, sd = .2)
total = c(m_0_sd_1,m_1_sd_1,m_2_sd_1)
pretty_gmm = normalmixEM(total,k = 3)
plot(pretty_gmm,which = 2)
table(apply(pretty_gmm$posterior,1,which.max))

##plot stuff for determining rank
convergence = c(2776,1971,1419,1135,857,710,553,471,394,335,307,264,254,237,181,134,118,83,70,62,62,46,44,36,25,23,22,15,13,11,9,6,4,3,2,2,1,0)
runtime     = c(4.352,6.155,6.372,11.027,11.532,10.948,14.342,13.894,14.233,14.532,15.095,14.233,10.254,13.533,14.042,16.143,17.221,17.982,18.853,17.254,18.241,18.143,18.644,19.241,20.991,21.43,24.241,22.4213,25.2143,24.42,25.214,22.1341,21.341,17.2542,19.1252,18.982,19.143,18.152)
plot(c(1:length(convergence)),convergence,ylab = "MSE",xlab = "rank", xlim = c(1,40))

par(new = T)
plot(c(1:length(runtime)),runtime,axes = F, ylab = NA, xlab = NA, cex = 1.2, pch = 16)
axis(side = 4)
mtext(side = 4,line = 3, "Number Selected")

legend("topleft",legend=c("MSE","Runtime (Seconds)"),pch = c(1,16))
#Running the algy
if (FALSE){
  LAMBDA = 0
  for( j in seq(5,50,5)){
    FEATURE_COUNT = j
    for (i in c(1:10)){
      LAMBDA = i
      diploid_incomp =maskValues(diploid,MISSING_PERCENTAGE)
      my_diploid = ridgeImpute(diploid_incomp,FEATURE_COUNT,LAMBDA)
      real_values = diploid_matrix[which(diploid_incomp == -1)]
      imputed_values = my_diploid[which(diploid_incomp == -1)]
      acc = length(which(imputed_values == real_values))/length(imputed_values)
      print(sprintf("With feature count %d and LAMDA %d F! score is %f and accuracy %f",FEATURE_COUNT,LAMBDA,F1score(imputed_values,real_values),acc))
    }
  }
}


