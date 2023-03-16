list.Scenarios_25_25_M2M_index = list()
list.Scenarios_25_25_M2M_Microbiotes = list()
list.Scenarios_25_25_M2M_Metabolites = list()
list.Scenarios_25_25_M2M_effects = list()

for(rep in 1:100){
  print(rep)
  L0 = matrix(0, ncol=25, nrow=25)
  #Variance for microbiotes throughout samples
  #We can play on this parameter in order to assess biological variability
  diag(L0) = runif(25,1.5,2.5)
  
  #Off-diagonal elements are randomly selected to have either 0 covariance or a positive or negative covariance based on uniform distribution
  L0[lower.tri(L0)] = sapply(1:length(L0[lower.tri(L0)]), function(x) sample(c(0,runif(1,-1.5,1.5)),1, prob = c(0.7,0.3)))
  
  Precision0 = L0%*%t(L0)
  
  #We obtain Covariance matrix based on Cholesky decomposition of lower triangular matrix
  Sigma0 = solve(Precision0)
  Cor0 = cov2cor(Sigma0)
  
  multi.norm = MASS::mvrnorm(100, rep(0,25), Cor0)
  
  L1 = matrix(0, ncol=25, nrow=25)
  #Variance for microbiotes throughout samples
  #We can play on this parameter in order to assess biological variability
  diag(L1) = runif(25,1.5,2.5)
  
  #Off-diagonal elements are randomly selected to have either 0 covariance or a positive or negative covariance based on uniform distribution
  L1[lower.tri(L1)] = sapply(1:length(L1[lower.tri(L1)]), function(x) sample(c(0,runif(1,1,1.5)),1, prob = c(0.7,0.3)))
  
  Precision1 = L1%*%t(L1)
  
  #We obtain Covariance matrix based on Cholesky decomposition of lower triangular matrix
  Sigma1 = solve(Precision1)
  Cor1 = cov2cor(Sigma1)
  
  #Adding random noise  in case to obtain semi-positive matrix
  #Cor0 = Cor0 + diag(ncol(Cor0))*0.01
  
  #The multivariate normal distribution is generated for 100 individuals with mean 0 and the 
  #Correlation structure 
  multi.norm1 = MASS::mvrnorm(100, rep(0,25), Cor1)
  
  
  index.associated.microbiote = sample(1:25, round(25*round(runif(1,0.03,0.2), 2)))
  index.impacted.metabolite = sample(1:25, round(25*round(runif(1,0.03,0.2), 2)))
  
  
  index.associated.metabolite = sample(c(1:25)[-index.impacted.metabolite], round(25*round(runif(1,0.03,0.2), 2)))
  index.impacted.microbiote = sample(c(1:25)[-index.associated.microbiote], round(25*round(runif(1,0.03,0.2), 2)))
  
  n.associated.by.metabolite = sapply(index.impacted.metabolite, function(x) sample(1:length(index.associated.microbiote),1))
  lhg = lapply(1:length(index.impacted.metabolite), function(x) sample(index.associated.microbiote,n.associated.by.metabolite[x]))
  lhg.beta = lapply(lhg, function(x) rnorm(length(x), 0,0.2))
  
  
  n.associated.by.microbiote = sapply(index.impacted.microbiote, function(x) sample(1:length(index.associated.metabolite),1))
  lhh = lapply(1:length(index.impacted.microbiote), function(x) sample(index.associated.metabolite,n.associated.by.microbiote[x]))
  lhh.beta = lapply(lhh, function(x) rnorm(length(x), 0,0.2))
  
  mu.microbiotes = replicate(25,runif(100,0,1))
  mu.microbiotes[,index.associated.microbiote] = replicate(length(index.associated.microbiote), runif(100,1.5,2.5))
  
  mu.metabolites = replicate(25, runif(100,0,1))
  mu.metabolites[,index.associated.metabolite] = replicate(length(index.associated.metabolite), runif(100,1.5,2.5))
  
  size.microbiotes = runif(25,0.5,1)
  size.metabolites = runif(25,0.3,1)
  
  
  simulated.microbiotes = matrix(qnbinom(pnorm(multi.norm), mu=mu.microbiotes, size=size.microbiotes), ncol=25,nrow=100)
  simulated.metabolites = matrix(qnbinom(pnorm(multi.norm1), mu=mu.metabolites, size=size.metabolites), ncol=25,nrow=100)
  
  mu.microbiotes[,index.impacted.microbiote]=sapply(1:length(index.impacted.microbiote) , function(x) exp(simulated.metabolites[,lhh[[x]], drop=F]%*%lhh.beta[[x]]))
  mu.metabolites[,index.impacted.metabolite]=sapply(1:length(index.impacted.metabolite) , function(x) exp(simulated.microbiotes[,lhg[[x]], drop=F]%*%lhg.beta[[x]]))
  
  simulated.microbiotes = matrix(qnbinom(pnorm(multi.norm), mu=mu.microbiotes, size=size.microbiotes), ncol=25,nrow=100)
  simulated.metabolites = matrix(qnbinom(pnorm(multi.norm1), mu=mu.metabolites, size=size.metabolites), ncol=25,nrow=100)
  
  
  names(lhg) = index.impacted.metabolite
  names(lhg.beta) = index.impacted.metabolite
  
  names(lhh) = index.impacted.microbiote
  names(lhh.beta) = index.impacted.microbiote
  
  
  list.Scenarios_25_25_M2M_index[[rep]] = list("Microbiote" = lhh, "Metabolite" = lhg)
  list.Scenarios_25_25_M2M_effects[[rep]] = list("Microbiote" = lhh.beta, "Metabolite" = lhg.beta)
  list.Scenarios_25_25_M2M_Microbiotes[[rep]] = simulated.microbiotes
  list.Scenarios_25_25_M2M_Metabolites[[rep]] = simulated.metabolites
  
}

summary(pscl::zeroinfl(list.Scenarios_25_25_M2M_Metabolites[[1]][,3]~list.Scenarios_25_25_M2M_Microbiotes[[1]][,25]+list.Scenarios_25_25_M2M_Microbiotes[[1]][,13]+list.Scenarios_25_25_M2M_Microbiotes[[1]][,21],dist = "negbin"))
summary(pscl::zeroinfl(list.Scenarios_25_25_M2M_Microbiotes[[1]][,14]~list.Scenarios_25_25_M2M_Metabolites[[1]][,25]+
                         list.Scenarios_25_25_M2M_Metabolites[[1]][,11]+list.Scenarios_25_25_M2M_Metabolites[[1]][,20]+list.Scenarios_25_25_M2M_Metabolites[[1]][,22],dist = "negbin"))

cc.test = CCA::cc(list.Scenarios_25_25_M2M_Metabolites[[24]], list.Scenarios_25_25_M2M_Microbiotes[[24]])

redundancy.test = compute.redundancy.cca(list.Scenarios_25_25_M2M_Metabolites[[24]], list.Scenarios_25_25_M2M_Microbiotes[[24]],
                       cc.test)


sum(redundancy.test$ConditionalRedundancy$X)
