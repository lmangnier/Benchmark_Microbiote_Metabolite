#Code to assess Correlation and Regression-based approaches for scenarios 25-25 

load("util_functions.R")

#Correlation approaches:
#Spearman correlation
list.scenarios.dense.Microbiotes = readRDS("../data_scenarios_denses_fixed_Microbiotes.RDS")
list.scenarios.dense.Metabolites = readRDS("../data_scenarios_denses_fixed_Metabolites.RDS")
list.scenarios.dense.index = readRDS("../data_scenarios_denses_fixed_index.RDS")

list.scenarios.dense.random.Microbiotes = readRDS("../data_scenarios_denses_random_Microbiotes.RDS")
list.scenarios.dense.random.Metabolites = readRDS("../data_scenarios_denses_random_Metabolites.RDS")
list.scenarios.dense.random.index = readRDS("../data_scenarios_denses_random_index.RDS")

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


#Confusion Matrix:


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
        sapply(1:1000, function(rep) f1.score(confusion.matrix(pvalues.ACAT.Scenario10_beta_mixed[[rep]], list.scenarios.dense.random.index$Scenario_10_beta_random[[rep]]$Microbiotes, c(1:25)[-list.scenarios.dense.random.index$Scenario_10_beta_random[[rep]]$Microbiotes], 0.002))))

#If we are interested in constructing ROC curve.
build.ROC(pvalues.ACAT.Scenario10_beta_mixed[[577]], list.scenarios.dense.random.index$Scenario_10_beta_random[[577]]$Microbiotes, c(1:25)[-list.scenarios.dense.random.index$Scenario_10_beta_random[[577]]$Microbiotes])

#Compute AUC
auc.Scenario10_beta_mixed = sapply(1:1000, function(x) compute.auc(pvalues.ACAT.Scenario10_beta_mixed[[x]],list.scenarios.dense.random.index$Scenario_10_beta_random[[x]]$Microbiotes))
auc.Scenario3_beta_mixed = sapply(1:1000, function(x) compute.auc(pvalues.ACAT.Scenario3_beta_mixed[[x]],list.scenarios.dense.random.index$Scenario_3_beta_random[[x]]$Microbiotes))
auc.Scenario1_beta_mixed = sapply(1:1000, function(x) compute.auc(pvalues.ACAT.Scenario1_beta_mixed[[x]],list.scenarios.dense.random.index$Scenario_1_beta_random[[x]]$Microbiotes))


boxplot(auc.Scenario1_beta_mixed,auc.Scenario3_beta_mixed,auc.Scenario10_beta_mixed)

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
        auc.Scenario3_beta_mixed)

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
#Here we use ZINB approach:
library(pscl)

zinb.Scenario3_beta_neg0.5 = list()
zinb.Scenario3_beta_neg0.25 = list()
zinb.Scenario3_beta_0.25 = list()
zinb.Scenario3_beta_0.5 = list()
zinb.Scenario3_beta_mixed = list()

for(r in 1:1000){
  print(r)
  l = list()
  for(j in 1:25){
    
    p = c()
    for(i in 1:25){
      
      m = summary(zeroinfl(list.scenarios.dense.Metabolites$Scenario_3_beta_0.25[[r]][,j]~list.scenarios.dense.Microbiotes$Scenario_3_beta_0.25[[r]][,i], dist="negbin"))
      p = c(p,m$coefficients$count[2,4])
    }
    l[[j]] = p
    zinb.Scenario3_beta_0.25[[r]] = l
    
  }
  
}


library(mpath)
list.scenarios.dense.index$Scenario_3_beta_0.5[[1]]
Xy = data.frame(cbind(list.scenarios.dense.Metabolites$Scenario_3_beta_0.5[[1]][,1], list.scenarios.dense.Microbiotes$Scenario_3_beta_0.5[[1]]))


fit.lasso = mpath::zipath(X1~.|.,data=Xy,family = "negbin", nlambda=100,
              lambda.zero.min.ratio=0.001, maxit.em=300, maxit.theta=25,
              theta.fixed=FALSE, trace=FALSE, penalty="enet", rescale=FALSE)

minBic = which.min(BIC(fit.lasso))
coef(fit.lasso, minBic)
