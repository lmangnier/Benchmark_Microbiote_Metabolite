#Code to assess Correlation and Regression-based approaches for scenarios 25-25 

#Correlation approaches:
#Spearman correlation

list.scenarios.dense.Microbiotes = readRDS("C:\\Users\\loicm\\Desktop\\Projets\\Mimint\\data\\data_scenarios_denses_fixed_Microbiotes.RDS")
list.scenarios.dense.Metabolites = readRDS("C:\\Users\\loicm\\Desktop\\Projets\\Mimint\\data\\data_scenarios_denses_fixed_Metabolites.RDS")
list.scenarios.dense.index = readRDS("C:\\Users\\loicm\\Desktop\\Projets\\Mimint\\data\\data_scenarios_denses_fixed_index.RDS")

list.scenarios.dense.random.Microbiotes = readRDS("C:\\Users\\loicm\\Desktop\\Projets\\Mimint\\data\\data_scenarios_denses_random_Microbiotes.RDS")
list.scenarios.dense.random.Metabolites = readRDS("C:\\Users\\loicm\\Desktop\\Projets\\Mimint\\data\\data_scenarios_denses_random_Metabolites.RDS")
list.scenarios.dense.random.index = readRDS("C:\\Users\\loicm\\Desktop\\Projets\\Mimint\\data\\data_scenarios_denses_random_index.RDS")

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


#For a same number of associated elements (3 microbiotes associated with 3 metabolites)
#Beta=-0.5
spearman.cor.Scenario_3_beta_neg0.5 = lapply(1:1000, function(rep) {
  compute.pairwise.spearman(list.scenarios.dense.Microbiotes$`Scenario_3_beta_-0.5`[[rep]], list.scenarios.dense.Metabolites$`Scenario_3_beta_-0.5`[[rep]])})
#Beta=-0.25
spearman.cor.Scenario_3_beta_neg0.25 = lapply(1:1000, function(rep) {
  compute.pairwise.spearman(list.scenarios.dense.Microbiotes$`Scenario_3_beta_-0.25`[[rep]], list.scenarios.dense.Metabolites$`Scenario_3_beta_-0.25`[[rep]])})
#Beta=0.25
spearman.cor.Scenario_3_beta_0.25 = lapply(1:1000, function(rep) {
  compute.pairwise.spearman(list.scenarios.dense.Microbiotes$Scenario_3_beta_0.25[[rep]], list.scenarios.dense.Metabolites$Scenario_3_beta_0.25[[rep]])})
#Beta=0.5
spearman.cor.Scenario_3_beta_0.5 = lapply(1:1000, function(rep) {
  compute.pairwise.spearman(list.scenarios.dense.Microbiotes$Scenario_3_beta_0.5[[rep]], list.scenarios.dense.Metabolites$Scenario_3_beta_0.5[[rep]])})
#Beta~N(0,0.2)
spearman.cor.Scenario_3_beta_mixed = lapply(1:1000, function(rep) {
  compute.pairwise.spearman(list.scenarios.dense.random.Microbiotes$Scenario_3_beta_random[[rep]], list.scenarios.dense.random.Metabolites$Scenario_3_beta_random[[rep]])})


#Jaccard Index
jaccard.Scenario3_beta_neg0.5 = sapply(1:1000, function(rep) {
  compute.jaccard(ACAT.by.element(spearman.cor.Scenario_3_beta_neg0.5[[rep]], list.scenarios.dense.Microbiotes$`Scenario_3_beta_-0.5`[[rep]],list.scenarios.dense.Metabolites$`Scenario_3_beta_-0.5`[[rep]]),
                  list.scenarios.dense.index$`Scenario_3_beta_-0.5`[[rep]]$Microbiotes)
})

jaccard.Scenario3_beta_neg0.25 = sapply(1:1000, function(rep) {
  compute.jaccard(ACAT.by.element(spearman.cor.Scenario_3_beta_neg0.25[[rep]], list.scenarios.dense.Microbiotes$`Scenario_3_beta_-0.25`[[rep]],list.scenarios.dense.Metabolites$`Scenario_3_beta_-0.25`[[rep]]),
                  list.scenarios.dense.index$`Scenario_3_beta_-0.25`[[rep]]$Microbiotes)
})

jaccard.Scenario3_beta_0.25 = sapply(1:1000, function(rep) {
  compute.jaccard(ACAT.by.element(spearman.cor.Scenario_3_beta_0.25[[rep]], list.scenarios.dense.Microbiotes$Scenario_3_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_3_beta_0.25[[rep]]),
                  list.scenarios.dense.index$Scenario_3_beta_0.25[[rep]]$Microbiotes)
})

jaccard.Scenario3_beta_0.5 = sapply(1:1000, function(rep) {
  compute.jaccard(ACAT.by.element(spearman.cor.Scenario_3_beta_0.5[[rep]], list.scenarios.dense.Microbiotes$Scenario_3_beta_0.5[[rep]],list.scenarios.dense.Metabolites$Scenario_3_beta_0.5[[rep]]),
                list.scenarios.dense.index$Scenario_3_beta_0.5[[rep]]$Microbiotes)
  })


jaccard.Scenario3_beta_mixed = sapply(1:1000, function(rep) {
  compute.jaccard(ACAT.by.element(spearman.cor.Scenario_3_beta_mixed[[rep]], list.scenarios.dense.random.Microbiotes$Scenario_3_beta_random[[rep]],list.scenarios.dense.random.Metabolites$Scenario_3_beta_random[[rep]]),
                  list.scenarios.dense.random.index$Scenario_3_beta_random[[rep]]$Microbiotes)
})

boxplot(jaccard.Scenario3_beta_neg0.5,jaccard.Scenario3_beta_neg0.25,jaccard.Scenario3_beta_0.25,jaccard.Scenario3_beta_0.5,jaccard.Scenario3_beta_mixed,
        names=c("Neg 0.5", "Neg 0.25", "0.25", "0.5", "Mixed"), horizontal=F)


#We can do the same now comparing different number of associated elements for the same level of association
#Let's take our mixed scenario as main results
spearman.cor.Scenario_1_beta_mixed = lapply(1:1000, function(rep) {
  compute.pairwise.spearman(list.scenarios.dense.random.Microbiotes$Scenario_1_beta_random[[rep]], list.scenarios.dense.random.Metabolites$Scenario_1_beta_random[[rep]])})
spearman.cor.Scenario_10_beta_mixed = lapply(1:1000, function(rep) {
  compute.pairwise.spearman(list.scenarios.dense.random.Microbiotes$Scenario_10_beta_random[[rep]], list.scenarios.dense.random.Metabolites$Scenario_10_beta_random[[rep]])})

jaccard.Scenario1_beta_mixed = sapply(1:1000, function(rep) {
  compute.jaccard(ACAT.by.element(spearman.cor.Scenario_1_beta_mixed[[rep]], list.scenarios.dense.random.Microbiotes$Scenario_1_beta_random[[rep]],list.scenarios.dense.random.Metabolites$Scenario_1_beta_random[[rep]]),
                  list.scenarios.dense.random.index$Scenario_1_beta_random[[rep]]$Microbiotes)
})

jaccard.Scenario10_beta_mixed = sapply(1:1000, function(rep) {
  compute.jaccard(ACAT.by.element(spearman.cor.Scenario_10_beta_mixed[[rep]], list.scenarios.dense.random.Microbiotes$Scenario_10_beta_random[[rep]],list.scenarios.dense.random.Metabolites$Scenario_10_beta_random[[rep]]),
                  list.scenarios.dense.random.index$Scenario_10_beta_random[[rep]]$Microbiotes)
})


boxplot(jaccard.Scenario1_beta_mixed,jaccard.Scenario3_beta_mixed,jaccard.Scenario10_beta_mixed, names=c(1,3,10))

jaccard.associations.Scenario1_beta_mixed = sapply(1:1000, function(rep) find.true.associations(spearman.cor.Scenario_1_beta_mixed[[rep]], 0.05/25, list.scenarios.dense.random.index$Scenario_1_beta_random[[rep]]$Microbiotes,list.scenarios.dense.random.index$Scenario_1_beta_random[[rep]]$Metabolites))
jaccard.associations.Scenario3_beta_mixed = sapply(1:1000, function(rep) find.true.associations(spearman.cor.Scenario_3_beta_mixed[[rep]], 0.05/25, list.scenarios.dense.random.index$Scenario_3_beta_random[[rep]]$Microbiotes,list.scenarios.dense.random.index$Scenario_3_beta_random[[rep]]$Metabolites))
jaccard.associations.Scenario10_beta_mixed = sapply(1:1000, function(rep) find.true.associations(spearman.cor.Scenario_10_beta_mixed[[rep]], 0.05/25, list.scenarios.dense.random.index$Scenario_10_beta_random[[rep]]$Microbiotes,list.scenarios.dense.random.index$Scenario_10_beta_random[[rep]]$Metabolites))

boxplot(jaccard.associations.Scenario1_beta_mixed,jaccard.associations.Scenario3_beta_mixed,jaccard.associations.Scenario10_beta_mixed, names=c(1,3,10))

jaccard.associations.Scenario3_beta_neg0.5 = sapply(1:1000, function(rep) find.true.associations(spearman.cor.Scenario_3_beta_neg0.5[[rep]], 0.05/25, list.scenarios.dense.index$`Scenario_3_beta_-0.5`[[rep]]$Microbiotes,list.scenarios.dense.index$`Scenario_3_beta_-0.5`[[rep]]$Metabolites))
jaccard.associations.Scenario3_beta_neg0.25 = sapply(1:1000, function(rep) find.true.associations(spearman.cor.Scenario_3_beta_neg0.25[[rep]], 0.05/25, list.scenarios.dense.index$`Scenario_3_beta_-0.25`[[rep]]$Microbiotes,list.scenarios.dense.index$`Scenario_3_beta_-0.25`[[rep]]$Metabolites))
jaccard.associations.Scenario3_beta_0.25 = sapply(1:1000, function(rep) find.true.associations(spearman.cor.Scenario_3_beta_0.25[[rep]], 0.05/25, list.scenarios.dense.index$Scenario_3_beta_0.25[[rep]]$Microbiotes,list.scenarios.dense.index$Scenario_3_beta_0.25[[rep]]$Metabolites))
jaccard.associations.Scenario3_beta_0.5 = sapply(1:1000, function(rep) find.true.associations(spearman.cor.Scenario_3_beta_0.5[[rep]], 0.05/25, list.scenarios.dense.index$Scenario_3_beta_0.5[[rep]]$Microbiotes,list.scenarios.dense.index$Scenario_3_beta_0.5[[rep]]$Metabolites))

boxplot(jaccard.associations.Scenario3_beta_neg0.5,jaccard.associations.Scenario3_beta_neg0.25,
        jaccard.associations.Scenario3_beta_0.25,jaccard.associations.Scenario3_beta_0.5,jaccard.associations.Scenario3_beta_mixed)

summary(jaccard.associations.Scenario3_beta_0.25)

#Regression approaches:
library(pscl)

l.r = list()

for(r in 1:100){
  print(r)
  p.all = c()
  for(j in 1:25){
    
    p = c()
    for(i in 1:25){
      
      m = summary(zeroinfl(list.scenarios.dense.Metabolites$`Scenario_3_beta_-0.25`[[r]][,j]~list.scenarios.dense.Microbiotes$`Scenario_3_beta_-0.25`[[r]][,i], dist="negbin"))
      p = c(p,m$coefficients$count[2,4])
    }
    
    p.all = c(p.all,ACAT::ACAT(p[!is.na(p)]))
  }
  l.r[[r]] = p.all
}
