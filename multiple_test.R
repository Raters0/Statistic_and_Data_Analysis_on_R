#create positive test, return a matrix
create_posi <- function(test_num=20,sample_size=3){
  test_mean <- sample(seq(1,3,0.01),test_num,replace = T) #change sensitivity here!
  test_sd <- sample(seq(0.2,2,0.01),test_num,replace = T)
  mean_vector <- rep(test_mean,sample_size)
  std_vector <- rep(test_sd,sample_size)
  posi_matrix <- matrix(rnorm(test_num*sample_size,mean=mean_vector,sd=std_vector),ncol=sample_size)  # add test data
  posi_matrix <- cbind(posi_matrix,rowMeans(posi_matrix[,1:sample_size]))  # add mean columns
  posi_matrix <- cbind(posi_matrix,apply(posi_matrix[,1:sample_size],1,sd))  # add std columns
  rownames(posi_matrix) <- paste('gene',as.character(c(1:test_num)),sep='') 
  colnames(posi_matrix) <- c(paste('test',as.character(c(1:sample_size)),sep=''),'mean','sd')
  return(posi_matrix)
}

#create negative test, return a matrix
create_nega <- function(test_num=80,sample_size=3){
  nega_matrix <- matrix(rnorm(test_num*sample_size,mean=0,sd=1),ncol=sample_size)  # add test data
  nega_matrix <- cbind(nega_matrix,rowMeans(nega_matrix[,1:sample_size]))  # add mean columns
  nega_matrix <- cbind(nega_matrix,apply(nega_matrix[,1:sample_size],1,sd))  # add std columns
  rownames(nega_matrix) <- paste('gene',as.character(c(1:test_num)),sep='')
  colnames(nega_matrix) <- c(paste('test',as.character(c(1:sample_size)),sep=''),'mean','sd')
  return(nega_matrix)
}

# create a multiple test, and return a matrix
create_treatment <- function(total_num=100,posi_num=20,sample_size=3){
  posi_matrix <- create_posi(posi_num,sample_size)
  nega_matrix <- create_nega(total_num-posi_num,sample_size)
  result_matrix <- rbind(posi_matrix,nega_matrix)  # combine two matrix vertically
  rownames(result_matrix) <- paste('gene',as.character(c(1:total_num)),sep='')
  return(result_matrix)
}

# create a two-sample t test, return a list containing p-values, treatment-matrix, control matrix
t_simulation <- function(pn,nn,tss=3,css=3){
  # pn = positive number, need to >= 2
  # nn = negative number, need to >= 2
  # tss = treatment_sample_size, need to >= 2
  # css = control_sample_size, need to >=2
  treatment_matrix <- create_treatment(pn+nn,pn,tss)
  control_matrix <- create_nega(pn+nn,css)
  s_diff <- (((treatment_matrix[,'sd']^2)*(tss-1))+((control_matrix[,'sd']^2)*(css-1)))/(tss+css-2)  # two-sample t test
  t <- (treatment_matrix[,'mean'] - control_matrix[,'mean'])/sqrt(s_diff/tss+s_diff/css)
  p <- 2*pt(-abs(t),df=tss+css-2)  # two-sided t test
  return(list(p=p,treatment_matrix=treatment_matrix,control_matrix=control_matrix))
}

# Bonferroni method correction
Bon_method <- function(result_list,total_posi_names,fwer){
  fwers <- p.adjust(result_list$p, method="bonferroni")
  test_posi <- rownames(result_list$treatment_matrix)[fwers<fwer]
  cat('Bonferroni method:\n')
  cat('1) hit numbers:\n')
  print(length(test_posi))
  cat('2) true positive:\n')
  print(sum(test_posi %in% total_posi_names))
  cat('3) total positive:\n')
  print(length(total_posi_names))
  cat('4) Power: \n')
  print(length(test_posi[test_posi %in% total_posi_names])/length(total_posi_names))
  cat('5) FDR: \n')
  print(1-sum(test_posi %in% total_posi_names)/length(test_posi))
  cat('------------------------\n')
}

# Benjamini-Hochberg method correction
BH_method <- function(result_list,total_posi_names,fdr){
  fdrs <- p.adjust(result_list$p, method="BH")
  test_posi <- rownames(result_list$treatment_matrix)[fdrs<fdr]
  cat('BH method:\n')
  cat('1) hit numbers:\n')
  print(length(test_posi))
  cat('2) true positive:\n')
  print(sum(test_posi %in% total_posi_names))
  cat('3) total positive:\n')
  print(length(total_posi_names))
  cat('4) Power: \n')
  print(length(test_posi[test_posi %in% total_posi_names])/length(total_posi_names))
  cat('5) FDR: \n')
  print(1-sum(test_posi %in% total_posi_names)/length(test_posi))
  cat('------------------------\n')
}

# Storey's method correction
Storey_method <- function(result_list,total_posi_names,fdr){
  library(qvalue)
  qobj <- qvalue(result_list$p)
  fdrs <- qobj$qvalues
  test_posi <- rownames(result_list$treatment_matrix)[fdrs<fdr]
  cat('Storey method:\n')
  cat('1) hit numbers:\n')
  print(length(test_posi))
  cat('2) true positive:\n')
  print(sum(test_posi %in% total_posi_names))
  cat('3) total positive:\n')
  print(length(total_posi_names))
  cat('4) Power: \n')
  print(length(test_posi[test_posi %in% total_posi_names])/length(total_posi_names))
  cat('5) FDR: \n')
  print(1-sum(test_posi %in% total_posi_names)/length(test_posi))
  cat('------------------------\n')
}

#main
result <- t_simulation(2000,8000,3,3)  # change sample size here
p <- result$p  # positive/negative number >= 2
total_posi <- rownames(result$treatment_matrix[1:2000,])  # return the names of true-positive tests 
hist(p,breaks=seq(0,1,0.01))

# Bonferroni method
Bon_method(result,total_posi,0.3)

# BH method
BH_method(result,total_posi,0.1)

#Storey method
Storey_method(result,total_posi,0.1)