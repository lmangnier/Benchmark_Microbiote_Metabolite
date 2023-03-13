#Main scenarios for methods when P+Q>>N
#Beta= [-0.5,-0.25,0.25,0.5], 100/100 microbiotes-metabolites
#1000 replicates are generated under each scenario
set.seed(123456)

list.scenarios.sparse.Microbiotes = list()
list.scenarios.sparse.Metabolites = list()
list.scenarios.sparse.index = list()

for(m in c(1,3,10)){
  print(m)
  for(beta in c(-0.5, -0.25, 0.25, 0.5)){
    print(beta)
    list.microbiotes = list()
    list.metabolites = list()
    list.index = list()
    for(i in 1:1000){
      L0 = matrix(0, ncol=100, nrow=100)
      #Variance for microbiotes throughout samples
      #We can play on this parameter in order to assess biological variability
      diag(L0) = runif(100,1.5,2.5)
      
      #Off-diagonal elements are randomly selected to have either 0 covariance or a positive or negative covariance based on uniform distribution
      L0[lower.tri(L0)] = sapply(1:length(L0[lower.tri(L0)]), function(x) sample(c(0,runif(1,-1.5,1.5)),1, prob = c(0.7,0.3)))
      
      Precision0 = L0%*%t(L0)
      
      #We obtain Covariance matrix based on Cholesky decomposition of lower triangular matrix
      Sigma0 = solve(Precision0)
      Cor0 = cov2cor(Sigma0)
      
      multi.norm = MASS::mvrnorm(100, rep(0,100), Cor0)
      
      index.associated.microbiote = sample(1:100, m)
      index.associated.metabolite = sample(1:100,m)
      
      # beta.coef = rnorm(3,0,0.2)
      # print(beta.coef)
      list.index[[i]] =  list("Microbiotes"=index.associated.microbiote,"Metabolites"=index.associated.metabolite)
      
      mu.microbiotes = replicate(100,runif(100,0,1))
      mu.microbiotes[,index.associated.microbiote] = replicate(length(index.associated.microbiote), runif(100,1,2))
      
      # matrix.mu.microbiotes = MASS::mvrnorm(1000,mu.microbiotes, Sigma = diag(1,100))
      size.microbiotes = runif(100,0.5,1)
      simulated.microbiotes = matrix(qnbinom(pnorm(multi.norm), mu=mu.microbiotes, size=size.microbiotes), ncol=100,nrow=100)
      
      list.microbiotes[[i]] = simulated.microbiotes
      
      matrix.mu.metabolites = exp(MASS::mvrnorm(100,rep(0,100), Sigma = diag(0.2,100)))
      
      #For each associated microbiote we assume a log-linear model of the mean
      vec.mu.for.associated.metabolites = sapply(1:length(index.associated.microbiote), function(x){
        
        exp(beta*simulated.microbiotes[,index.associated.microbiote[x]])*matrix.mu.metabolites[,index.associated.metabolite[x]]
      })
      
      matrix.mu.metabolites[,index.associated.metabolite] = vec.mu.for.associated.metabolites
      
      size.metabolites = runif(100,0.5,1)
      
      L1 = matrix(0, ncol=100, nrow=100)
      #Variance for microbiotes throughout samples
      #We can play on this parameter in order to assess biological variability
      diag(L1) = runif(100,1.5,2.5)
      
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
      multi.norm1 = MASS::mvrnorm(100, rep(0,100), Cor1)
      
      #Metabolites are simulated with common overdispersion and mean depending on microbiotes
      simulated.metabolites = matrix(qnbinom(pnorm(multi.norm1), mu=matrix.mu.metabolites, size=size.metabolites), ncol=100,nrow=100)
      list.metabolites[[i]] = simulated.metabolites
    }
    list.scenarios.sparse.index[[paste0("Scenario_",m,"_beta_",beta)]] = list.index
    list.scenarios.sparse.Microbiotes[[paste0("Scenario_",m,"_beta_",beta)]] = list.microbiotes
    list.scenarios.sparse.Metabolites[[paste0("Scenario_",m,"_beta_",beta)]] = list.metabolites
  }
  
}


saveRDS(list.scenarios.sparse.index, "C:\\Users\\loicm\\Desktop\\Projets\\Mimint\\data\\data_scenarios_sparse_fixed_index.RDS")
saveRDS(list.scenarios.sparse.Microbiotes, "C:\\Users\\loicm\\Desktop\\Projets\\Mimint\\data\\data_scenarios_sparse_fixed_Microbiotes.RDS")
saveRDS(list.scenarios.sparse.Metabolites, "C:\\Users\\loicm\\Desktop\\Projets\\Mimint\\data\\data_scenarios_sparse_fixed_Metabolites.RDS")

list.scenarios.dense.random.index = list()
list.scenarios.dense.random.Microbiotes = list()
list.scenarios.dense.random.Metabolites = list()

for(m in c(1,3,10)){
  print(m)
  list.index = list()
  list.metabolites = list()
  list.microbiotes = list()
  
  for(i in 1:1000){
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
    
    index.associated.microbiote = sample(1:25, m)
    index.associated.metabolite = sample(1:25,m)
    
    beta.coef = rnorm(m,0,0.2)
    
    list.index[[i]] =  list("Microbiotes"=index.associated.microbiote,"Metabolites"=index.associated.metabolite)
    
    mu.microbiotes = replicate(25,runif(100,0,1))
    mu.microbiotes[,index.associated.microbiote] = replicate(length(index.associated.microbiote), runif(25,1,2))
    
    # matrix.mu.microbiotes = MASS::mvrnorm(1000,mu.microbiotes, Sigma = diag(1,100))
    size.microbiotes = runif(25,0.5,1)
    simulated.microbiotes = matrix(qnbinom(pnorm(multi.norm), mu=mu.microbiotes, size=size.microbiotes), ncol=25,nrow=100)
    
    list.microbiotes[[i]] = simulated.microbiotes
    
    matrix.mu.metabolites = exp(MASS::mvrnorm(100,rep(0,25), Sigma = diag(0.2,25)))
    
    #For each associated microbiote we assume a log-linear model of the mean
    vec.mu.for.associated.metabolites = sapply(1:length(index.associated.microbiote), function(x){
      
      exp(beta*simulated.microbiotes[,index.associated.microbiote[x]])*matrix.mu.metabolites[,index.associated.metabolite[x]]
    })
    
    matrix.mu.metabolites[,index.associated.metabolite] = vec.mu.for.associated.metabolites
    
    size.metabolites = runif(25,0.5,1)
    
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
    
    #Metabolites are simulated with common overdispersion and mean depending on microbiotes
    simulated.metabolites = matrix(qnbinom(pnorm(multi.norm1), mu=matrix.mu.metabolites, size=size.metabolites), ncol=25,nrow=100)
    list.metabolites[[i]] = simulated.metabolites
  }
  list.scenarios.dense.random.index[[paste0("Scenario_",m,"_beta_random")]] = list.index
  list.scenarios.dense.random.Microbiotes[[paste0("Scenario_",m,"_beta_random")]] = list.microbiotes
  list.scenarios.dense.random.Metabolites[[paste0("Scenario_",m,"_beta_random")]] = list.metabolites
}

saveRDS(list.scenarios.dense.random.index, "../data/data_scenarios_denses_random_index.RDS")
saveRDS(list.scenarios.dense.random.Microbiotes, "../data/data_scenarios_denses_random_Microbiotes.RDS")
saveRDS(list.scenarios.dense.random.Metabolites, "../data/data_scenarios_denses_random_Metabolites.RDS")
