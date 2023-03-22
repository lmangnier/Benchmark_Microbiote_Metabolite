#Code to assess Correlation and Regression-based approaches for scenarios 25-25 
library(mpath)
par(mar=rep(4,4))

load("util_functions.R")

#Correlation approaches:
#Spearman correlation
list.scenarios.dense.Microbiotes = readRDS("data/data_scenarios_denses_fixed_Microbiotes.RDS")
list.scenarios.dense.Metabolites = readRDS("data/data_scenarios_denses_fixed_Metabolites.RDS")
list.scenarios.dense.index = readRDS("data/data_scenarios_denses_fixed_index.RDS")

list.scenarios.dense.random.Microbiotes = readRDS("data/data_scenarios_denses_random_Microbiotes.RDS")
list.scenarios.dense.random.Metabolites = readRDS("data/data_scenarios_denses_random_Metabolites.RDS")
list.scenarios.dense.random.index = readRDS("data/data_scenarios_denses_random_index.RDS")

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
                  list.scenarios.dense.index$`Scenario_3_beta_-0.5`[[rep]]$Microbiotes,0.05/25)
})

jaccard.Scenario3_beta_neg0.25 = sapply(1:1000, function(rep) {
  compute.jaccard(ACAT.by.element(spearman.cor.Scenario_3_beta_neg0.25[[rep]], list.scenarios.dense.Microbiotes$`Scenario_3_beta_-0.25`[[rep]],list.scenarios.dense.Metabolites$`Scenario_3_beta_-0.25`[[rep]]),
                  list.scenarios.dense.index$`Scenario_3_beta_-0.25`[[rep]]$Microbiotes,0.05/25)
})

jaccard.Scenario3_beta_0.25 = sapply(1:1000, function(rep) {
  compute.jaccard(ACAT.by.element(spearman.cor.Scenario_3_beta_0.25[[rep]], list.scenarios.dense.Microbiotes$Scenario_3_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_3_beta_0.25[[rep]]),
                  list.scenarios.dense.index$Scenario_3_beta_0.25[[rep]]$Microbiotes,0.05/25)
})

jaccard.Scenario3_beta_0.5 = sapply(1:1000, function(rep) {
  compute.jaccard(ACAT.by.element(spearman.cor.Scenario_3_beta_0.5[[rep]], list.scenarios.dense.Microbiotes$Scenario_3_beta_0.5[[rep]],list.scenarios.dense.Metabolites$Scenario_3_beta_0.5[[rep]]),
                list.scenarios.dense.index$Scenario_3_beta_0.5[[rep]]$Microbiotes,0.05/25)
  })


jaccard.Scenario3_beta_mixed = sapply(1:1000, function(rep) {
  compute.jaccard(ACAT.by.element(spearman.cor.Scenario_3_beta_mixed[[rep]], list.scenarios.dense.random.Microbiotes$Scenario_3_beta_random[[rep]],list.scenarios.dense.random.Metabolites$Scenario_3_beta_random[[rep]]),
                  list.scenarios.dense.random.index$Scenario_3_beta_random[[rep]]$Microbiotes, 0.05/25)
})

boxplot(jaccard.Scenario3_beta_neg0.5,jaccard.Scenario3_beta_neg0.25,jaccard.Scenario3_beta_0.25,jaccard.Scenario3_beta_0.5,jaccard.Scenario3_beta_mixed,
        names=c("Neg 0.5", "Neg 0.25", "0.25", "0.5", "Mixed"), horizontal=F, main="Jaccard Index for Scenarios with 3 associated elements at different association strenghts (alpha=0.05/25)", xlab="Association Strengths", ylab="Jaccard Index")


#We can do the same now comparing different number of associated elements for the same level of association
#Let's take our mixed scenario as main results
spearman.cor.Scenario_1_beta_mixed = lapply(1:1000, function(rep) {
  compute.pairwise.spearman(list.scenarios.dense.random.Microbiotes$Scenario_1_beta_random[[rep]], list.scenarios.dense.random.Metabolites$Scenario_1_beta_random[[rep]])})
spearman.cor.Scenario_10_beta_mixed = lapply(1:1000, function(rep) {
  compute.pairwise.spearman(list.scenarios.dense.random.Microbiotes$Scenario_10_beta_random[[rep]], list.scenarios.dense.random.Metabolites$Scenario_10_beta_random[[rep]])})

jaccard.Scenario1_beta_mixed = sapply(1:1000, function(rep) {
  compute.jaccard(ACAT.by.element(spearman.cor.Scenario_1_beta_mixed[[rep]], list.scenarios.dense.random.Microbiotes$Scenario_1_beta_random[[rep]],list.scenarios.dense.random.Metabolites$Scenario_1_beta_random[[rep]]),
                  list.scenarios.dense.random.index$Scenario_1_beta_random[[rep]]$Microbiotes, 0.05/25)
})

jaccard.Scenario10_beta_mixed = sapply(1:1000, function(rep) {
  compute.jaccard(ACAT.by.element(spearman.cor.Scenario_10_beta_mixed[[rep]], list.scenarios.dense.random.Microbiotes$Scenario_10_beta_random[[rep]],list.scenarios.dense.random.Metabolites$Scenario_10_beta_random[[rep]]),
                  list.scenarios.dense.random.index$Scenario_10_beta_random[[rep]]$Microbiotes, 0.05/25)
})


boxplot(jaccard.Scenario1_beta_mixed,jaccard.Scenario3_beta_mixed,jaccard.Scenario10_beta_mixed, names=c(1,3,10), main="Jaccard Index for Scenarios with mixed association strenghts (alpha=0.05/25)", xlab="#Associations", ylab="Jaccard Index")



#Confusion Matrix, Sensitivity, Specificity, F1-Score
pvalues.ACAT.Scenario10_beta_mixed = lapply(1:1000, function(rep){
  ACAT.by.element(spearman.cor.Scenario_10_beta_mixed[[rep]], list.scenarios.dense.random.Microbiotes$Scenario_10_beta_random[[rep]],list.scenarios.dense.random.Metabolites$Scenario_10_beta_random[[rep]])
})

pvalues.ACAT.Scenario3_beta_mixed = lapply(1:1000, function(rep){
  ACAT.by.element(spearman.cor.Scenario_3_beta_mixed[[rep]], list.scenarios.dense.random.Microbiotes$Scenario_3_beta_random[[rep]],list.scenarios.dense.random.Metabolites$Scenario_3_beta_random[[rep]])
})

pvalues.ACAT.Scenario1_beta_mixed = lapply(1:1000, function(rep){
  ACAT.by.element(spearman.cor.Scenario_1_beta_mixed[[rep]], list.scenarios.dense.random.Microbiotes$Scenario_1_beta_random[[rep]],list.scenarios.dense.random.Metabolites$Scenario_1_beta_random[[rep]])
})


boxplot(sapply(1:1000, function(rep) f1.score(confusion.matrix(pvalues.ACAT.Scenario1_beta_mixed[[rep]], list.scenarios.dense.random.index$Scenario_1_beta_random[[rep]]$Microbiotes, c(1:25)[-list.scenarios.dense.random.index$Scenario_1_beta_random[[rep]]$Microbiotes], 0.002))),
        sapply(1:1000, function(rep) f1.score(confusion.matrix(pvalues.ACAT.Scenario3_beta_mixed[[rep]], list.scenarios.dense.random.index$Scenario_3_beta_random[[rep]]$Microbiotes, c(1:25)[-list.scenarios.dense.random.index$Scenario_3_beta_random[[rep]]$Microbiotes], 0.002))),
        sapply(1:1000, function(rep) f1.score(confusion.matrix(pvalues.ACAT.Scenario10_beta_mixed[[rep]], list.scenarios.dense.random.index$Scenario_10_beta_random[[rep]]$Microbiotes, c(1:25)[-list.scenarios.dense.random.index$Scenario_10_beta_random[[rep]]$Microbiotes], 0.002))), ylim=c(0,1), names=c(1,3,10), ylab="F1 Score", xlab="#Associations", main="Jaccard Index for Scenarios with mixed association strenghts (alpha=0.05/25)")

#If we are interested in constructing ROC curve.
#Worst scenario
par(mfrow=c(1,2))
build.ROC(pvalues.ACAT.Scenario10_beta_mixed[[577]], list.scenarios.dense.random.index$Scenario_10_beta_random[[577]]$Microbiotes, c(1:25)[-list.scenarios.dense.random.index$Scenario_10_beta_random[[577]]$Microbiotes])
build.ROC(pvalues.ACAT.Scenario10_beta_mixed[[198]], list.scenarios.dense.random.index$Scenario_10_beta_random[[198]]$Microbiotes, c(1:25)[-list.scenarios.dense.random.index$Scenario_10_beta_random[[198]]$Microbiotes])
par(mfrow=c(1,1))

#Compute AUC
auc.Scenario10_beta_mixed = sapply(1:1000, function(x) compute.auc(pvalues.ACAT.Scenario10_beta_mixed[[x]],list.scenarios.dense.random.index$Scenario_10_beta_random[[x]]$Microbiotes))
auc.Scenario3_beta_mixed = sapply(1:1000, function(x) compute.auc(pvalues.ACAT.Scenario3_beta_mixed[[x]],list.scenarios.dense.random.index$Scenario_3_beta_random[[x]]$Microbiotes))
auc.Scenario1_beta_mixed = sapply(1:1000, function(x) compute.auc(pvalues.ACAT.Scenario1_beta_mixed[[x]],list.scenarios.dense.random.index$Scenario_1_beta_random[[x]]$Microbiotes))


boxplot(auc.Scenario1_beta_mixed,auc.Scenario3_beta_mixed,auc.Scenario10_beta_mixed, ylab="AUC", xlab="#Associations",names=c(1,3,10), main="AUC for Scenarios with mixed association strenghts")

#Comparisons with other association levels:
pvalues.ACAT.Scenario3_beta_neg0.5 = lapply(1:1000, function(rep){
  ACAT.by.element(spearman.cor.Scenario_3_beta_neg0.5[[rep]], list.scenarios.dense.Microbiotes$`Scenario_3_beta_-0.5`[[rep]],list.scenarios.dense.Metabolites$`Scenario_3_beta_-0.5`[[rep]])
})

pvalues.ACAT.Scenario3_beta_neg0.25 = lapply(1:1000, function(rep){
  ACAT.by.element(spearman.cor.Scenario_3_beta_neg0.25[[rep]], list.scenarios.dense.Microbiotes$`Scenario_3_beta_-0.25`[[rep]],list.scenarios.dense.Metabolites$`Scenario_3_beta_-0.25`[[rep]])
})

pvalues.ACAT.Scenario3_beta_0.25 = lapply(1:1000, function(rep){
  ACAT.by.element(spearman.cor.Scenario_3_beta_0.25[[rep]], list.scenarios.dense.Microbiotes$Scenario_3_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_3_beta_0.25[[rep]])
})

pvalues.ACAT.Scenario3_beta_0.5 = lapply(1:1000, function(rep){
  ACAT.by.element(spearman.cor.Scenario_3_beta_0.5[[rep]], list.scenarios.dense.Microbiotes$Scenario_3_beta_0.5[[rep]],list.scenarios.dense.Metabolites$Scenario_3_beta_0.5[[rep]])
})

auc.Scenario3_beta_neg0.5 = sapply(1:1000, function(x) compute.auc(pvalues.ACAT.Scenario3_beta_neg0.5[[x]],list.scenarios.dense.index$`Scenario_3_beta_-0.5`[[x]]$Microbiotes))
auc.Scenario3_beta_neg0.25 = sapply(1:1000, function(x) compute.auc(pvalues.ACAT.Scenario3_beta_neg0.25[[x]],list.scenarios.dense.index$`Scenario_3_beta_-0.25`[[x]]$Microbiotes))
auc.Scenario3_beta_0.25 = sapply(1:1000, function(x) compute.auc(pvalues.ACAT.Scenario3_beta_0.25[[x]],list.scenarios.dense.index$Scenario_3_beta_0.25[[x]]$Microbiotes))
auc.Scenario3_beta_0.5 = sapply(1:1000, function(x) compute.auc(pvalues.ACAT.Scenario3_beta_0.5[[x]],list.scenarios.dense.index$Scenario_3_beta_0.5[[x]]$Microbiotes))

boxplot(auc.Scenario3_beta_neg0.5,
        auc.Scenario3_beta_neg0.25,
        auc.Scenario3_beta_0.25,
        auc.Scenario3_beta_0.5,
        auc.Scenario3_beta_mixed, ylim=c(0,1), ylab="AUC", xlab="Association Strenghts",names=c("Neg 0.5", "Neg 0.25", "0.25", "0.5", "Mixed"), main="AUC for Scenarios with 3 associated elements at different association strenghts")

jaccard.associations.Scenario1_beta_mixed = sapply(1:1000, function(rep) find.true.associations(spearman.cor.Scenario_1_beta_mixed[[rep]], 0.05/25, list.scenarios.dense.random.index$Scenario_1_beta_random[[rep]]$Microbiotes,list.scenarios.dense.random.index$Scenario_1_beta_random[[rep]]$Metabolites))
jaccard.associations.Scenario3_beta_mixed = sapply(1:1000, function(rep) find.true.associations(spearman.cor.Scenario_3_beta_mixed[[rep]], 0.05/25, list.scenarios.dense.random.index$Scenario_3_beta_random[[rep]]$Microbiotes,list.scenarios.dense.random.index$Scenario_3_beta_random[[rep]]$Metabolites))
jaccard.associations.Scenario10_beta_mixed = sapply(1:1000, function(rep) find.true.associations(spearman.cor.Scenario_10_beta_mixed[[rep]], 0.05/25, list.scenarios.dense.random.index$Scenario_10_beta_random[[rep]]$Microbiotes,list.scenarios.dense.random.index$Scenario_10_beta_random[[rep]]$Metabolites))

boxplot(jaccard.associations.Scenario1_beta_mixed,jaccard.associations.Scenario3_beta_mixed,jaccard.associations.Scenario10_beta_mixed, names=c(1,3,10), ylim=c(0,1), main="Jaccard Index at pair levels for Scenarios with mixed association strenghts (alpha=0.05/25)",
        xlab="#Associations", ylab="Jaccard Index")

jaccard.associations.Scenario3_beta_neg0.5 = sapply(1:1000, function(rep) find.true.associations(spearman.cor.Scenario_3_beta_neg0.5[[rep]], 0.05/25, list.scenarios.dense.index$`Scenario_3_beta_-0.5`[[rep]]$Microbiotes,list.scenarios.dense.index$`Scenario_3_beta_-0.5`[[rep]]$Metabolites))
jaccard.associations.Scenario3_beta_neg0.25 = sapply(1:1000, function(rep) find.true.associations(spearman.cor.Scenario_3_beta_neg0.25[[rep]], 0.05/25, list.scenarios.dense.index$`Scenario_3_beta_-0.25`[[rep]]$Microbiotes,list.scenarios.dense.index$`Scenario_3_beta_-0.25`[[rep]]$Metabolites))
jaccard.associations.Scenario3_beta_0.25 = sapply(1:1000, function(rep) find.true.associations(spearman.cor.Scenario_3_beta_0.25[[rep]], 0.05/25, list.scenarios.dense.index$Scenario_3_beta_0.25[[rep]]$Microbiotes,list.scenarios.dense.index$Scenario_3_beta_0.25[[rep]]$Metabolites))
jaccard.associations.Scenario3_beta_0.5 = sapply(1:1000, function(rep) find.true.associations(spearman.cor.Scenario_3_beta_0.5[[rep]], 0.05/25, list.scenarios.dense.index$Scenario_3_beta_0.5[[rep]]$Microbiotes,list.scenarios.dense.index$Scenario_3_beta_0.5[[rep]]$Metabolites))

boxplot(jaccard.associations.Scenario3_beta_neg0.5,jaccard.associations.Scenario3_beta_neg0.25,
        jaccard.associations.Scenario3_beta_0.25,jaccard.associations.Scenario3_beta_0.5,jaccard.associations.Scenario3_beta_mixed, ylim=c(0,1), 
        main="Jaccard Index at pair levels for Scenarios with 3 associated at different association strenghts (alpha=0.05/25)",
        xlab="Association Strenghts", ylab="Jaccard Index", names=c("Neg 0.5", "Neg 0.25", "0.25", "0.5", "Mixed"))

summary(jaccard.associations.Scenario3_beta_0.25)

#Regression approaches:
#Here we use ZINB approach (100 replicates):
library(pscl)

zinb.Scenario3_beta_neg0.5 = list()
zinb.Scenario3_beta_neg0.25 = list()
zinb.Scenario3_beta_0.25 = list()
zinb.Scenario3_beta_0.5 = list()
zinb.Scenario3_beta_mixed = list()

for(r in 1:100){
  print(r)
  l = list()
  for(j in 1:25){
    
    p = c()
    for(i in 1:25){
      
      m = summary(zeroinfl(list.scenarios.dense.Metabolites$`Scenario_3_beta_-0.5`[[r]][,j]~list.scenarios.dense.Microbiotes$`Scenario_3_beta_-0.5`[[r]][,i], dist="negbin"))
      p = c(p,m$coefficients$count[2,4])
    }
    l[[j]] = p
    zinb.Scenario3_beta_neg0.5[[r]] = l
    
  }
  
}

zinb.Scenario3_beta_neg0.5[[10]]
zinb.Scenario3_beta_neg0.25[[10]]
zinb.Scenario3_beta_0.25[[10]]
zinb.Scenario3_beta_0.5[[10]]
zinb.Scenario3_beta_mixed[[10]]

#LASSO-ZINB (100 replicates)
list.results.lasso.Scenario_1_beta_0.25 = list()
list.results.lasso.Scenario_3_beta_0.25 = list()
list.results.lasso.Scenario_10_beta_0.25 = list()

for(rep in 25:26){
  print(rep)
  X = list.scenarios.dense.Microbiotes$Scenario_3_beta_0.25[[rep]]
  
  l = list()
  for(y in 1:25){
    Xy = data.frame(cbind(list.scenarios.dense.Metabolites$Scenario_3_beta_0.25[[rep]][,y], X))
    fit.lasso.NB = mpath::glmregNB(X1~.,data=Xy, penalty="enet", rescale=FALSE, parallel = T)
    
    
    # fit.lasso = mpath::zipath(X1~.|.,data=Xy,family = "negbin", nlambda=100,
    #                           lambda.zero.min.ratio=0.001, maxit.em=300, maxit.theta=25,
    #                           theta.fixed=FALSE, trace=FALSE, penalty="enet", rescale=FALSE)
    # 
    
    #Each model are selected based on BIC
    minBic = which.min(BIC(fit.lasso.NB))
    
    l[[y]] = coef(fit.lasso.NB, minBic)
  }
  
  list.results.lasso.Scenario_3_beta_0.25[[rep]] = l
}

list.results.lasso.Scenario_3_beta_0.25[[26]]
find.true.associations.lasso(list.results.lasso.Scenario_3_beta_0.25[[3]], list.scenarios.dense.index$Scenario_3_beta_0.25[[3]]$Microbiotes,list.scenarios.dense.index$Scenario_3_beta_0.25[[3]]$Metabolites )




