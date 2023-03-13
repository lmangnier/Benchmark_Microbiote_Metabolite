library(CCA)
library(vegan)
library(mixOmics)

list.scenarios.dense.Microbiotes = readRDS("../data_scenarios_denses_fixed_Microbiotes.RDS")
list.scenarios.dense.Metabolites = readRDS("../data_scenarios_denses_fixed_Metabolites.RDS")
list.scenarios.dense.index = readRDS("../data_scenarios_denses_fixed_index.RDS")

list.scenarios.dense.random.Microbiotes = readRDS("../data_scenarios_denses_random_Microbiotes.RDS")
list.scenarios.dense.random.Metabolites = readRDS("../data_scenarios_denses_random_Metabolites.RDS")
list.scenarios.dense.random.index = readRDS("../data_scenarios_denses_random_index.RDS")


#CCA
cca.Scenario_3_beta_neg0.5 = lapply(1:1000, function(rep){
  cc(list.scenarios.dense.Microbiotes$`Scenario_3_beta_-0.5`[[rep]],list.scenarios.dense.Metabolites$`Scenario_3_beta_-0.5`[[rep]])
})
cca.Scenario_3_beta_neg0.25 = lapply(1:1000, function(rep){
  cc(list.scenarios.dense.Microbiotes$`Scenario_3_beta_-0.25`[[rep]],list.scenarios.dense.Metabolites$`Scenario_3_beta_-0.25`[[rep]])
})
cca.Scenario_3_beta_0.25 = lapply(1:1000, function(rep){
  cc(list.scenarios.dense.Microbiotes$Scenario_3_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_3_beta_0.25[[rep]])
})
cca.Scenario_3_beta_0.5 = lapply(1:1000, function(rep){
  cc(list.scenarios.dense.Microbiotes$Scenario_3_beta_0.5[[rep]],list.scenarios.dense.Metabolites$Scenario_3_beta_0.5[[rep]])
})
cca.Scenario_3_beta_mixed = lapply(1:1000, function(rep){
  cc(list.scenarios.dense.random.Microbiotes$Scenario_3_beta_random[[rep]],list.scenarios.dense.random.Metabolites$Scenario_3_beta_random[[rep]])
})

cca.Scenario_1_beta_0.25 = lapply(1:1000, function(rep){
  cc(list.scenarios.dense.Microbiotes$Scenario_1_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_1_beta_0.25[[rep]])
})

cca.Scenario_10_beta_0.25 = lapply(1:1000, function(rep){
  cc(list.scenarios.dense.Microbiotes$Scenario_10_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_10_beta_0.25[[rep]])
})

redundancy.Scenario_3_beta_neg0.5 = lapply(1:1000, function(rep){
  compute.redundancy.cca(list.scenarios.dense.Microbiotes$`Scenario_3_beta_-0.5`[[rep]],list.scenarios.dense.Metabolites$`Scenario_3_beta_-0.5`[[rep]],cca.Scenario_3_beta_neg0.5[[rep]])})

redundancy.Scenario_3_beta_neg0.25 = lapply(1:1000, function(rep){
  compute.redundancy.cca(list.scenarios.dense.Microbiotes$`Scenario_3_beta_-0.25`[[rep]],list.scenarios.dense.Metabolites$`Scenario_3_beta_-0.25`[[rep]],cca.Scenario_3_beta_neg0.25[[rep]])})

redundancy.Scenario_3_beta_0.25 = lapply(1:1000, function(rep){
  compute.redundancy.cca(list.scenarios.dense.Microbiotes$Scenario_3_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_3_beta_0.25[[rep]],cca.Scenario_3_beta_0.25[[rep]])})

redundancy.Scenario_3_beta_0.5 = lapply(1:1000, function(rep){
  compute.redundancy.cca(list.scenarios.dense.Microbiotes$Scenario_3_beta_0.5[[rep]],list.scenarios.dense.Metabolites$Scenario_3_beta_0.5[[rep]],cca.Scenario_3_beta_0.5[[rep]])})

redundancy.Scenario_3_beta_mixed = lapply(1:1000, function(rep){
  compute.redundancy.cca(list.scenarios.dense.random.Microbiotes$Scenario_3_beta_random[[rep]],list.scenarios.dense.random.Metabolites$Scenario_3_beta_random[[rep]],cca.Scenario_3_beta_mixed[[rep]])})


redundancy.Scenario_1_beta_0.25 = lapply(1:1000, function(rep){
  compute.redundancy.cca(list.scenarios.dense.Microbiotes$Scenario_1_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_1_beta_0.25[[rep]],cca.Scenario_1_beta_0.25[[rep]])})

redundancy.Scenario_10_beta_0.25 = lapply(1:1000, function(rep){
  compute.redundancy.cca(list.scenarios.dense.Microbiotes$Scenario_10_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_10_beta_0.25[[rep]],cca.Scenario_10_beta_0.25[[rep]])})

#Conditional Redundancy CCA X
conditional.redundancy.X.Scenario_3_beta_neg0.5 = sapply(redundancy.Scenario_3_beta_neg0.5, function(x) sum(x$ConditionalRedundancy$X))
conditional.redundancy.X.Scenario_3_beta_neg0.25 = sapply(redundancy.Scenario_3_beta_neg0.25, function(x) sum(x$ConditionalRedundancy$X))
conditional.redundancy.X.Scenario_3_beta_0.25 = sapply(redundancy.Scenario_3_beta_0.25, function(x) sum(x$ConditionalRedundancy$X))
conditional.redundancy.X.Scenario_3_beta_0.5 = sapply(redundancy.Scenario_3_beta_0.5, function(x) sum(x$ConditionalRedundancy$X))
conditional.redundancy.X.Scenario_3_beta_mixed = sapply(redundancy.Scenario_3_beta_mixed, function(x) sum(x$ConditionalRedundancy$X))

#Conditional Redundancy CCA Y
conditional.redundancy.Y.Scenario_3_beta_neg0.5 = sapply(redundancy.Scenario_3_beta_neg0.5, function(x) sum(x$ConditionalRedundancy$Y))
conditional.redundancy.Y.Scenario_3_beta_neg0.25 = sapply(redundancy.Scenario_3_beta_neg0.25, function(x) sum(x$ConditionalRedundancy$Y))
conditional.redundancy.Y.Scenario_3_beta_0.25 = sapply(redundancy.Scenario_3_beta_0.25, function(x) sum(x$ConditionalRedundancy$Y))
conditional.redundancy.Y.Scenario_3_beta_0.5 = sapply(redundancy.Scenario_3_beta_0.5, function(x) sum(x$ConditionalRedundancy$Y))
conditional.redundancy.Y.Scenario_3_beta_mixed = sapply(redundancy.Scenario_3_beta_mixed, function(x) sum(x$ConditionalRedundancy$Y))


conditional.redundancy.X.Scenario_1_beta_0.25 = sapply(redundancy.Scenario_1_beta_0.25, function(x) sum(x$ConditionalRedundancy$X))
conditional.redundancy.X.Scenario_10_beta_0.25 = sapply(redundancy.Scenario_10_beta_0.25, function(x) sum(x$ConditionalRedundancy$X))

conditional.redundancy.Y.Scenario_1_beta_0.25 = sapply(redundancy.Scenario_1_beta_0.25, function(x) sum(x$ConditionalRedundancy$Y))
conditional.redundancy.Y.Scenario_10_beta_0.25 = sapply(redundancy.Scenario_10_beta_0.25, function(x) sum(x$ConditionalRedundancy$Y))

par(mfrow=c(1,2))
par(mar=c(1,1,1,1))
boxplot(conditional.redundancy.X.Scenario_3_beta_neg0.5,
        conditional.redundancy.X.Scenario_3_beta_neg0.25,
        conditional.redundancy.X.Scenario_3_beta_0.25,
        conditional.redundancy.X.Scenario_3_beta_0.5,
        conditional.redundancy.X.Scenario_3_beta_mixed, names=c(-0.5, -0.25, 0.25, 0.5, "Mixed"), ylim=c(0,0.5))

boxplot(conditional.redundancy.Y.Scenario_3_beta_neg0.5,
        conditional.redundancy.Y.Scenario_3_beta_neg0.25,
        conditional.redundancy.Y.Scenario_3_beta_0.25,
        conditional.redundancy.Y.Scenario_3_beta_0.5,
        conditional.redundancy.Y.Scenario_3_beta_mixed, names=c(-0.5, -0.25, 0.25, 0.5, "Mixed"), ylim=c(0,0.5))


boxplot(conditional.redundancy.X.Scenario_1_beta_0.25,conditional.redundancy.X.Scenario_3_beta_0.25,
        conditional.redundancy.X.Scenario_10_beta_0.25)

boxplot(conditional.redundancy.Y.Scenario_1_beta_0.25,conditional.redundancy.Y.Scenario_3_beta_0.25,
        conditional.redundancy.Y.Scenario_10_beta_0.25, ylim=c(0,1))

#PLS: regression

pls.regression.Scenario_3_beta_0.25 = lapply(1:1000, function(rep){
  mixOmics::pls(list.scenarios.dense.Microbiotes$Scenario_3_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_3_beta_0.25[[rep]], ncomp=25, mode="regression")
})

pls.regression.Scenario_1_beta_0.25 = lapply(1:1000, function(rep){
  mixOmics::pls(list.scenarios.dense.Microbiotes$Scenario_1_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_1_beta_0.25[[rep]], ncomp=25, mode="regression")
})

pls.regression.Scenario_10_beta_0.25 = lapply(1:1000, function(rep){
  mixOmics::pls(list.scenarios.dense.Microbiotes$Scenario_10_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_10_beta_0.25[[rep]], ncomp=25, mode="regression")
})

#PLS: canonical

pls.canonical.Scenario_3_beta_0.25 = lapply(1:1000, function(rep){
  mixOmics::pls(list.scenarios.dense.Microbiotes$Scenario_3_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_3_beta_0.25[[rep]], ncomp=25, mode="canonical")
})

pls.canonical.Scenario_1_beta_0.25 = lapply(1:1000, function(rep){
  mixOmics::pls(list.scenarios.dense.Microbiotes$Scenario_1_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_1_beta_0.25[[rep]], ncomp=25, mode="canonical")
})

pls.canonical.Scenario_10_beta_0.25 = lapply(1:1000, function(rep){
  mixOmics::pls(list.scenarios.dense.Microbiotes$Scenario_10_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_10_beta_0.25[[rep]], ncomp=25, mode="canonical")
})

#redundancy for PLS
redundancy.pls.Scenario3_beta_0.25 = lapply(1:1000, function(rep) {
  compute.redundancy.pls(list.scenarios.dense.Microbiotes$Scenario_3_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_3_beta_0.25[[rep]],pls.regression.Scenario_3_beta_0.25[[rep]])
})

redundancy.pls.Scenario1_beta_0.25 = lapply(1:1000, function(rep) {
  compute.redundancy.pls(list.scenarios.dense.Microbiotes$Scenario_1_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_1_beta_0.25[[rep]],pls.regression.Scenario_1_beta_0.25[[rep]])
})

redundancy.pls.Scenario10_beta_0.25 = lapply(1:1000, function(rep) {
  compute.redundancy.pls(list.scenarios.dense.Microbiotes$Scenario_10_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_10_beta_0.25[[rep]],pls.regression.Scenario_10_beta_0.25[[rep]])
})


conditional.redundancy.pls.X.Scenario_1_beta_0.25 = sapply(redundancy.pls.Scenario1_beta_0.25, function(x) sum(x$ConditionalRedundancy$X))
conditional.redundancy.pls.X.Scenario_3_beta_0.25 = sapply(redundancy.pls.Scenario3_beta_0.25, function(x) sum(x$ConditionalRedundancy$X))
conditional.redundancy.pls.X.Scenario_10_beta_0.25 = sapply(redundancy.pls.Scenario10_beta_0.25, function(x) sum(x$ConditionalRedundancy$X))

conditional.redundancy.pls.Y.Scenario_1_beta_0.25 = sapply(redundancy.pls.Scenario1_beta_0.25, function(x) sum(x$ConditionalRedundancy$Y))
conditional.redundancy.pls.Y.Scenario_3_beta_0.25 = sapply(redundancy.pls.Scenario3_beta_0.25, function(x) sum(x$ConditionalRedundancy$Y))
conditional.redundancy.pls.Y.Scenario_10_beta_0.25 = sapply(redundancy.pls.Scenario10_beta_0.25, function(x) sum(x$ConditionalRedundancy$Y))


boxplot(conditional.redundancy.pls.X.Scenario_1_beta_0.25,conditional.redundancy.pls.X.Scenario_3_beta_0.25,
        conditional.redundancy.pls.X.Scenario_10_beta_0.25)

boxplot(conditional.redundancy.pls.Y.Scenario_1_beta_0.25,conditional.redundancy.pls.Y.Scenario_3_beta_0.25,
        conditional.redundancy.pls.Y.Scenario_10_beta_0.25, ylim=c(0,1))


#RDA
rda.Scenario_1_beta_0.25 = lapply(1:1000, function(rep){
  vegan::rda(list.scenarios.dense.Microbiotes$Scenario_1_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_1_beta_0.25[[rep]])
})
rda.Scenario_3_beta_0.25 = lapply(1:1000, function(rep){
  vegan::rda(list.scenarios.dense.Microbiotes$Scenario_3_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_3_beta_0.25[[rep]])
})
rda.Scenario_10_beta_0.25 = lapply(1:1000, function(rep){
  vegan::rda(list.scenarios.dense.Microbiotes$Scenario_10_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_10_beta_0.25[[rep]])
})

boxplot(sapply(1:1000, function(rep) compute.explained.variance.rda(rda.Scenario_1_beta_0.25[[rep]])),
sapply(1:1000, function(rep) compute.explained.variance.rda(rda.Scenario_3_beta_0.25[[rep]])),
sapply(1:1000, function(rep) compute.explained.variance.rda(rda.Scenario_10_beta_0.25[[rep]])), ylim=c(0,1))



trw.X = sapply(1:1000, function(rep) {
  compute.jaccard.pls(pls.regression.Scenario_10_beta_0.25[[rep]],
                    list.scenarios.dense.index$Scenario_10_beta_0.25[[rep]]$Microbiotes, type="X", index.component=1:25)})

trw.Y = sapply(1:1000, function(rep) {
  compute.jaccard.pls(pls.regression.Scenario_10_beta_0.25[[rep]],
                      list.scenarios.dense.index$Scenario_10_beta_0.25[[rep]]$Metabolites, type="Y", index.component=1:25)})

boxplot(t(trw.X), ylim=c(0,1))
boxplot(t(trw.Y), ylim=c(0,1))


dft = sapply(1:1000, function(rep) compute.jaccard.cca(cca.Scenario_10_beta_0.25[[rep]],list.scenarios.dense.index$Scenario_10_beta_0.25[[rep]]$Microbiotes, type="X", index.component=1:25))
boxplot(t(dft))

dft.y = sapply(1:1000, function(rep) compute.jaccard.cca(cca.Scenario_10_beta_0.25[[rep]],list.scenarios.dense.index$Scenario_10_beta_0.25[[rep]]$Metabolites, type="Y", index.component=1:25))
boxplot(t(dft.y))
