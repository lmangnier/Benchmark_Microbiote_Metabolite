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


size.microbiotes = runif(25,0.5,1)
size.metabolites = runif(25,0.5,1)

#Microbiotes & Metabolites  are simulated with common overdispersion and mean depending on microbiotes
simulated.metabolites = matrix(qnbinom(pnorm(multi.norm1), mu=matrix.mu.metabolites, size=size.metabolites), ncol=100,nrow=100)
simulated.microbiotes = matrix(qnbinom(pnorm(multi.norm), mu=mu.microbiotes, size=size.microbiotes), ncol=100,nrow=100)


index.associated.microbiote = sample(1:25, round(25*round(runif(1,0.03,0.4), 2)))
index.associated.metabolite = sample(1:25,round(25*round(runif(1,0.03,0.4), 2)))


#Many to many with variable effect
n.associated.by.microbiote = sapply(index.associated.microbiote, function(x) sample(1:length(index.associated.metabolite), 1))
n.associated.by.metabolite = sapply(index.associated.metabolite, function(x) sample(1:length(index.associated.microbiote), 1))


lki=lapply(1:length(index.associated.microbiote), function(x) sample(index.associated.metabolite, n.associated.by.microbiote[x]))
lkj=lapply(1:length(index.associated.metabolite), function(x) sample(index.associated.microbiote, n.associated.by.metabolite[x]))

#If coefficients are drawn from a gaussian 
beta.coef.microbiote = rnorm(length(index.associated.microbiote),0,0.2)
beta.coef.metabolite = rnorm(length(index.associated.metabolite), 0, 0.2)

mu.microbiotes = replicate(25,runif(100,0,1))
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

