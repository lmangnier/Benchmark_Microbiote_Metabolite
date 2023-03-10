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
        conditional.redundancy.Y.Scenario_10_beta_0.25)

#PLS



pls.regression.Scenario_3_beta_0.25 = lapply(1:1000, function(rep){
  mixOmics::pls(list.scenarios.dense.Microbiotes$Scenario_3_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_3_beta_0.25[[rep]], ncomp=25, mode="regression")
})

pls.regression.Scenario_1_beta_0.25 = lapply(1:1000, function(rep){
  mixOmics::pls(list.scenarios.dense.Microbiotes$Scenario_1_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_1_beta_0.25[[rep]], ncomp=25, mode="regression")
})

pls.regression.Scenario_10_beta_0.25 = lapply(1:1000, function(rep){
  mixOmics::pls(list.scenarios.dense.Microbiotes$Scenario_10_beta_0.25[[rep]],list.scenarios.dense.Metabolites$Scenario_10_beta_0.25[[rep]], ncomp=25, mode="regression")
})

#RDA

