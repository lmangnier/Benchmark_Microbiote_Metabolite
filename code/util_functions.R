#Util functions
#Univariate methods
compute.pairwise.spearman = function(X,Y){
  
  ncol.X = ncol(X)
  ncol.Y = ncol(Y)
  
  lapply(seq_len(ncol.X), function(x){
    lapply(seq_len(ncol.Y), function(y){
      cor.test(X[,x], Y[,y], method="spearman")
    })
  })
  
}

#Assess global association by microbiote, while adjusting for multiplicity (ACAT)
ACAT.by.element = function(cor.test.by.element, X, Y){
  
  ncol.X = ncol(X)
  ncol.Y = ncol(Y)
  
  sapply(seq_len(ncol.X), function(x){
    ACAT::ACAT(sapply(seq_len(ncol.Y), function(y){
      cor.test.by.element[[x]][[y]]$p.value
    }))
  })
}

#Compute Jaccard index
compute.jaccard = function(p.values, index.true.associated, threshold.signi=0.05){
  
  index.signi.p.values = which(p.values<=threshold.signi)
  
  sum(index.signi.p.values%in%index.true.associated)/
    (length(index.signi.p.values)+length(index.true.associated)-sum(index.signi.p.values%in%index.true.associated))
}

find.true.associations = function(cor.test.by.element,threshold.signi=0.05, index.X, index.Y){
  titi = lapply(cor.test.by.element, function(x) {
    which(sapply(x, function(y){
      y$p.value<=threshold.signi
    }))
  })
  
  pairs = list()
  
  for(i in 1:length(index.X)){
    a = c(index.X[i],
          index.Y[i])
    pairs[[i]] = a
  }
  
  index.non.null = which(sapply(titi, function(x) length(x))>0)
  
  l=c()
  
  for(i in index.non.null){
    
    u = c()
    for(j in titi[[i]]){
      c = c(i,j)
      n = c()
      for(p in 1:length(index.X)){
        n = c(n,all(c  == pairs[[p]]))
      }
      
      u = c(u,sum(n))
      
    }
    
    l = c(l,sum(n))
  }
  
  sum(l)/(length(index.X)+length(l)-sum(l))
}

confusion.matrix = function(pvalues, true.associated, true.non.associated ,signi.threshold=0.05){
  
  significant.associations = which(pvalues<=signi.threshold)
  non.significant.associations = which(pvalues>signi.threshold)
  
  TP = sum(significant.associations%in%true.associated)
  TN = sum(non.significant.associations%in%true.non.associated)
  FP = sum(significant.associations%in%true.non.associated)
  FN = sum(non.significant.associations%in%true.associated)
  
  return(matrix(c(TP, FN,FP, TN), nrow=2, ncol=2, dimnames = list(c("Pred_Pos","Pred_Neg"), c("True_Pos","True_Neg"))))
}

sensivity.specificity.from.confusion.matrix = function(confusion.matrix){
  
  sensitivity = confusion.matrix[1,1]/(confusion.matrix[1,1]+confusion.matrix[2,1])
  specificity = confusion.matrix[2,2]/(confusion.matrix[2,2]+confusion.matrix[1,2])
  
  return(c(sensitivity,specificity))
}
       
f1.score = function(confusion.matrix){
  m = sensivity.specificity.from.confusion.matrix(confusion.matrix)
  
  return(2/((1/m[1])+(1/m[2])))
}       

build.ROC = function(pvalues, true.associated, true.non.associated, signi.levels = seq(0,1, 0.0001)){
  
  confusion.matrix.for.all.signi = lapply(signi.levels, function(rho) confusion.matrix(pvalues,true.associated,true.non.associated,rho))
  
  specificity.sensitivity.for.all.confusion.matrix = sapply(confusion.matrix.for.all.signi, function(x){
    sensivity.specificity.from.confusion.matrix(x)
  })
  
  plot(1-specificity.sensitivity.for.all.confusion.matrix[2,],specificity.sensitivity.for.all.confusion.matrix[1,], xlab="1-Specificity", ylab="Sensitivity")
  lines(1-specificity.sensitivity.for.all.confusion.matrix[2,],specificity.sensitivity.for.all.confusion.matrix[1,], lwd=2)
  abline(0,1, lty=2, lwd=2)
}

compute.auc = function(pvalues,true.associated){
  ert = 1:length(pvalues)
  
  ert = as.numeric(ert%in%true.associated)
  auc(roc(ert,pvalues))
}

#Lasso-ZINB 


#CCA

#This function returns different types of redundancy from a CCA object
compute.redundancy.cca = function(X,Y,cca.object){
  
  #Proportion of variance by each variate
  variance.explained.X = colSums(cca.object$scores$corr.X.xscores^2)/ncol(X)
  variance.explained.Y = colSums(cca.object$scores$corr.Y.yscores^2)/ncol(Y)
  
  #Proportion of variance in the first pair explained by other member of the pair
  proportion.variance.canonical.cor = cca.object$cor^2
  
  #Proportion of X explained by the canonical variates 
  cond.redundancy.X = variance.explained.X*proportion.variance.canonical.cor
  cond.redundancy.Y = variance.explained.Y*proportion.variance.canonical.cor
  
  
  return(list("Redundancy" = list("X"=variance.explained.X, "Y"=variance.explained.Y),"VarianceCanonicalCor" = proportion.variance.canonical.cor,"ConditionalRedundancy" = list("X" = cond.redundancy.X, "Y"=cond.redundancy.Y)))
}


#PLS
compute.redundancy.pls = function(X,Y,pls.object){
  
  
  X.scaled = apply(X, 2, function(col) (col-mean(col))/sd(col))
  Y.scaled = apply(Y, 2, function(col) (col-mean(col))/sd(col))
  
  #Proportion of variance by each variate
  variance.explained.X = sapply(1:ncol(X), function(y) sum(sapply(1:ncol(X), function(x) cor(X.scaled[,x], pls.object$variates$X[,y])^2))/ncol(X))
  variance.explained.Y = sapply(1:ncol(Y), function(y) sum(sapply(1:ncol(Y), function(x) cor(Y.scaled[,x], pls.object$variates$Y[,y])^2))/ncol(Y))
  
  #Proportion of variance in the first pair explained by other member of the pair
  proportion.variance.canonical.cor = sapply(1:ncol(X), function(x) cor(pls.object$variates$X[,x], pls.object$variates$Y[,x])^2)
  
  #Proportion of variance explained by the canonical variates
  cond.redundancy.X = variance.explained.X * proportion.variance.canonical.cor
  cond.redundancy.Y = variance.explained.Y * proportion.variance.canonical.cor
  
  return(list("Redundancy" = list("X"=variance.explained.X, "Y"=variance.explained.Y),"VarianceCanonicalCor" = proportion.variance.canonical.cor,"ConditionalRedundancy" = list("X" = cond.redundancy.X, "Y"=cond.redundancy.Y)))
}
