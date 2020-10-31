# Quantitative model for the genetic releasing model
#
# Siyang Xia
# 2017.07.07
#
# For cluster analysis:
# 1. Remove all plotting and interactive analysis
# 2. Store all results in csv files
# 3. Conduct parallel computing
#
# Change a few things to try to improve efficiency:
# 1. Stop cycles to initial field / immigrant population when the new generation is less than 0.01% different from
#    previous generation;
# 2. Separate relapse function and initial function
# 3. Adapt the code to parallel computation.
#
# Modification from last version:
# 1._ Add external gene flow into the model - migration
# 2. Add a new outcome variables: proportional reduction of integrated vector competence
# 3. (Estimate the frequency of releasing after intense initial releasing to maintain the outcome)
#
# Assumptions:
# a. Infinitesimal model of quantitative traits: 
#    infinite number of small-effect, unlinked loci determine the genetic components of the phenotype, i.e. vector competence.
# b. Vector competence is determined by both genetic and environmental components;
# c. All loci are additive with equal effect (i.e. no dominant effects);
# d. Populations have non-overlapping generations;
# e. Populations have equal sex ratio and mating is random;
# f. Mutation rate is low enough to be ignored;
# g. Vector competence is under stabilizing natural selection in the field;
# h. The original population has reached the optimal vector competence (observed VC = optimal VC by natural selection)
# i. Constant releasing population: genotype of the releasing population is normally distributed.



###################################################################################################################

# Load the packages useful for constructing genotype distribution
# if (!require("truncnorm")) install.packages("truncnorm")
# if (!require("randomForest")) install.packages("randomForest")
# if (!require("rpart")) install.packages("rpart")
# if (!require("compiler")) install.packages("compiler")
# if (!require("doParallel")) install.packages("doParallel")
# if (!require("animation")) install.packages("animation")
# if (!require("grDevices")) install.packages("grDevices")

library(truncnorm)
# library(randomForest)
# library(rpart)
library(compiler)
library(doParallel)
# library(animation)
# library(grDevices)

packages_list <- c("truncnorm", "compiler", "doParallel")


#----------------------------------------------- Functions-------------------------------------
###### Toolbox functions

# Make Simpson's sequence for integration (h in the level of 10e-5):
Simpson <- function (dx, Nx){
  # Input:  - dx: the interval of the genotype sequence
  #         - Nx: the total number of genotype sequence (=length(vz))
  
  # Output: - sim: the Simpson's sequence with the same length as genotype sequence, for numerical integration
  #                Integration of f(x) = sum(sim * f(x)) 
  
  if ((Nx %% 2==1) && (Nx>2)){
    sim <- (dx/3)*c(1, rep(c(4,2), (Nx-1-2)/2), 4,1)
  }else if((Nx %% 3==1) && (Nx>3)){
    sim <- (dx/8)*c(3, rep(c(9,9,6),(Nx-1-3)/3),9,9,3)
  }else if (Nx>3){
    sim <- c((dx/3)*c(1, rep(c(4,2), (Nx-1-2-1)/2),4), dx*(1/3+1/2), dx/2)
  }else{
    sim <- (dx/2)*rep(1, Nx)
  }  
  return (sim)
}

# Create genotype distribution for the field population and the releasing population
Create <- function(pop, para){
  if(pop=="field"){
    return(para$N0 * dtruncnorm(x=para$vz, mean=para$fm, sd=para$fsd, a=para$zmin, b=para$zmax))
  }else if(pop=="releasing"){
    return(para$Nr * dtruncnorm(x=para$vz, mean=para$rm, sd=para$rsd, a=para$zmin, b=para$zmax))
  }
}


###### Life cycle functions:

# 1. Reproduction 
# 1.1 Do not consider the difference between males and females:
#    Use Fast Fourier Transform to calculate the reproduction outcome (provided by Prof.Baskett).

reproduction1 = function(phi, para) {
  # calculate reproduction
  # Input:  - phi: genotype population distribution at time t
  #         - para: parameter list
  
  # Output: - phi2: size/phenotye/genotype distribution for offspring
  
  Ntot = sum(para$Sz*phi) # total population size
  rphi = phi/Ntot # normalized genotype distribution
  rphi = c(rphi[para$midl],rep(0,para$Nz/2)) # zero pad
  phi2 = Re(fft(fft(rphi)*fft(rphi)*fft(para$Trans),inverse=T)*para$dz) 
          # integrate to get offspring genotypes
  phi2 = para$R*Ntot/2*phi2/sum(para$Sz*phi2) 
          # normalize then multiply by total production (female population size * reproductive rate per female)
          # and egg harching density independent survivalship
  return(phi2)
}

# 1.2 Consider the possibility that males and females may have different genotype/phenotype distribution:
# Because in practice, it is common that only males are released, resulting in the distribution of males being
# a mix of releasing population and field population, while females are only from field population.

# Superior function to calculate reprodution
reproduction2 <- function(densitym, densityf, para){
  Sz <- para$Sz
  size <- length(densitym)
  if (size!=length(densityf) | size!=length(para$vz)){
    return ("Error: Male distribution does not compatable with femal distribution!")
  }
  m <- densitym/Popsize(Sz, densitym)
  f <- densityf/Popsize(Sz, densityf)
  parents2 <- (m / sqrt(pi * para$vle)) %*% t(f)                    # The first three items in the formula      
  t2 <- 0.5 * (para$vz + rep(para$vz, each=size))                   # Matrix of (g1+g2)/2
  a <- double(size)
  for (i in seq_len(size)) {
    a[i] <- sum(parents2 * exp(-(para$vz[i] - t2)^2/para$vle) * (Sz%*%t(Sz)))
  }
  offdist <- (a/sum(a*para$Sz)) * Popsize(Sz, densityf) * para$R   # Standardize the distribution and multiply by population size and reproductive rate
  return(offdist)                                         
}
## Changes made on the newest version: the offspring population size is only determined by mother but not father pop size,
## as only female will lay eggs.
## -> releasing male will not influence pop size of the next generation



# 2. Density-dependent survival:
Den.dependent <- function(density, para){
  # Input: - density: population distribution of genotype
  #        - para: parameter list
  
  # Output: population genotype distribution after density dependence
  
  N <- Popsize(para$Sz, density)
  return (density/(1+para$alpha*N))
}

# 3. Density-independent survival for pupae emerging:
Den.independent <- function(density, para){
  # Input: - density: population distribution of genotype
  #        - para: parameter list
  
  # Output: population genotype distribution after density independence
  
  return(density*para$ind2)
}

# 4. Stabilizing selection:
Selection <- function (density, para){
  # Input: - density: population distribution of genotype
  #        - para: parameter list
  
  # Output: population genotype distribution after stabilizing selection
  
  return (sqrt(para$vs/(para$vs+para$ve)) * exp((-1/2)*(para$vz-para$fm)^2/(para$vs+para$ve)) * density)
}

# 5. Releasing (account for the relative fitness loss of the releasing strains):
# 5.1 Releasing both sex with 1:1 sex ratio (corresponding to reproduction1)
Releasing1 <- function(density, para){
  # Input: - density: population distribution of genotype
  #        - para: parameter list
  
  # Output: mixed population genotype distribution after releasing the lab population
  
  return (density + para$rel*para$fit)
}

# 5.2 Releasing only males (corresponding to reproduction2)
Releasing2 <- function(density, para){
  # Input: - density: population distribution of genotype
  #        - para: parameter list
  
  # Output: mixed population genotype distribution after releasing the lab population
  
  sex <- para$sex
  male <- (sex/(1+sex))*density + para$rel*para$fit
  female <- (1/(1+sex))*density
  return(list(male=male, female=female))
}

# 6. Migration
Migration <- function(density, para){
  # Input: - density: population distribution of genotype
  #        - para: parameter list
  
  # Output: mixed population genotype distribution after migration from outside of the study area
  
  return (density + para$migrant*para$Nm)
}


###### Output variables:
# Calculate population size:
Popsize <- function(simpson, density){
  # Input: - simpson: Simpson's sequence with the same length of genotype sequence
  #        - density: population distribution of genotype
  
  # Output: Population size
  
  return(sum(simpson*density))
}

# Calculate the phenotypic distribution 
GtoPh <- function(genotype, para){
  l <- length(para$vz)
  f_g <- (rep(1,l) %o% para$vz) - (para$vz %o% rep(1, l))
  prop <- (1/sqrt(2*pi*para$ve)) * exp(-0.5*(f_g^2)/para$ve)
  # prop <- as.matrix(dtruncnorm(x = as.vector(f_g), a = para$zmin, b = para$zmax, mean = 0, sd = sqrt(para$ve)), nrow=l)
  ng <- genotype %o% rep(1,l)
  phenotype <- apply(prop*ng*para$Sz, MARGIN = 2, FUN = sum)
  return(phenotype/sum(phenotype*para$Sz))      # Return standardized phenotype density distribution (sum=1)
}



# Calculate population mean genotype:
Genotype <- function (density, para){
  N <- Popsize(para$Sz, density)
  return(sum(para$Sz*para$vz*density/N))
}

# Calculate population genotype variance:
GVariance <- function (density, para){
  N <- Popsize(para$Sz, density)
  mu <- Genotype(density, para)
  return(sum((para$vz-mu)^2*density/N*para$Sz))
}

# Calculate population phenotype mean (mean vector competence):
Phenotype <- function(density, para){
  ph <- GtoPh(density, para)
  return(sum(para$vz*ph*para$Sz))
}

# Calculate population phenotype variation (variance of vector competence):
PhVariance <- function(density, para){
  ph <- GtoPh(density, para)
  ph_mean <- Phenotype(density, para)
  return(sum((para$vz-ph_mean)^2*ph*para$Sz))
}

# Calculate population mean fitness:
Fitness <- function (density, para){
  N <- Popsize(para$Sz,density)
  d <- density/N
  return(sum(Selection(d, para)*para$Sz))
}


##### Implement funtions:

## Initialize the field population and iterate it until reaching equilibrium

# a) Selection happen only before migration:
Initial_before <- function(density0, para, alter){
  # Input: - density0: initial genotype distribution of the field population
  #        - para: parameter list
  #        - alter: indicator about which reproduction assumption is used
  
  # Output: (equilibrium) genotype distribution of the field population before releasing
  
  # Initialization
  Sz <- para$Sz
  density <- density0
  
  # Changes over generations
  init <- 200                 # Maximum number of generations to run for the field population to reach equilibrium
  
  inipop <- rep(0,(init+1))
  inimean <- rep(0,(init+1))
  inivar <- rep(0,(init+1))
  inimean_vc <- rep(0,(init+1))
  inivar_vc <- rep(0,(init+1))
  inifit <- rep(0,(init+1))
  
  # The initial status:
  inipop[1] <- Popsize(Sz, density)
  inimean[1] <- Genotype(density, para)
  inivar[1] <- GVariance(density, para)
  if(para$phenotype){
    inimean_vc[1] <- Phenotype(density, para)
    inivar_vc[1] <- PhVariance(density, para)
  }
  inifit[1] <- Fitness(density, para)
  
  # start value of counter
  i <- 2
  stable <- FALSE
  
  while (i <= (init+1) & stable == FALSE){
    if(para$mig==1){
      density <- Migration(density, para)                          # Immigration
    }
    
    if (alter==1){
      density1 <- reproduction1(density, para)   
    }else if (alter==2){
      densitym <- (para$sex/(1+para$sex))*density
      densityf <- (1/(1+para$sex))*density
      density1 <- reproduction2(densitym, densityf, para)
    }                                                              # Reproduction (including density independence of egg hatching)
    
    density2 <- Den.dependent(density1, para)                      # Density dependent survival
    density3 <- Den.independent(density2, para)                    # Density independent survival
    density4 <- Selection(density3, para)                          # Stabilizing selection
    
    inipop[i] <- Popsize(Sz, density4)
    inimean[i] <- Genotype(density4, para)
    inivar[i] <- GVariance(density4, para)
    if(para$phenotype){
      inimean_vc[i] <- Phenotype(density4, para)
      inivar_vc[i] <- PhVariance(density4, para)
    }
    inifit[i] <- Fitness(density4, para)
    
    # Comparison between this cycle and previous one:
    c.change <- 0.00001                             # criteria: 0.01% of the previous value
    
    s1 <- (step.change(inipop, i) < c.change)
    s2 <- (step.change(inimean, i) < c.change)
    s3 <- (step.change(inivar, i) < c.change)
    if(para$phenotype){
      s4 <- (step.change(inimean_vc, i) < c.change)
      s5 <- (step.change(inivar_vc, i) < c.change)
    }else{
      s4 <- TRUE
      s5 <- TRUE
    }
    s6 <- all(abs(density4 - density) < c.change * density)
    
    # update:
    i <- i+1
    stable <- all(s1, s2, s3, s4, s5, s6)
    density <- density4
  }
  
  print(paste("The population reach stable in", i-1, "generations."))
  
  return(list(init=1:(i-1), density=density, 
              inipop=inipop, inimean=inimean, inivar=inivar, inisum=inipop * inimean,
              inimean_vc=inimean_vc, inivar_vc=inivar_vc, inifit=inifit, generation=i-1))
}


# b) Selection happen only after migration:
Initial_after <- function(density0, para, alter){
  # Input: - density0: initial genotype distribution of the field population
  #        - para: parameter list
  #        - alter: indicator about which reproduction assumption is used
  
  # Output: (equilibrium) genotype distribution of the field population before releasing
  
  # Initialization
  Sz <- para$Sz
  density <- density0
  
  # Changes over generations
  init <- 200                 # Maximum number of generations to run for the field population to reach equilibrium
  
  inipop <- rep(0,(init+1))
  inimean <- rep(0,(init+1))
  inivar <- rep(0,(init+1))
  inimean_vc <- rep(0,(init+1))
  inivar_vc <- rep(0,(init+1))
  inifit <- rep(0,(init+1))
  
  # The initial status:
  inipop[1] <- Popsize(Sz, density)
  inimean[1] <- Genotype(density, para)
  inivar[1] <- GVariance(density, para)
  if(para$phenotype){
    inimean_vc[1] <- Phenotype(density, para)
    inivar_vc[1] <- PhVariance(density, para)
  }
  inifit[1] <- Fitness(density, para)
  
  # start value of counter
  i <- 2
  stable <- FALSE
  
  while (i <= (init+1) & stable == FALSE){
    if(para$mig==1){
      density <- Migration(density, para)                         # Immigration
    }
    density1 <- Selection(density, para)                          # Stabilizing selection
    
    if (alter==1){
      density2 <- reproduction1(density1, para)   
    }else if (alter==2){
      densitym <- (para$sex/(1+para$sex))*density1
      densityf <- (1/(1+para$sex))*density1
      density2 <- reproduction2(densitym, densityf, para)
    }                                                              # Reproduction (including density independence of egg hatching)
    
    density3 <- Den.dependent(density2, para)                      # Density dependent survival
    density4 <- Den.independent(density3, para)                    # Density independent survival
    
    
    inipop[i] <- Popsize(Sz, density4)
    inimean[i] <- Genotype(density4, para)
    inivar[i] <- GVariance(density4, para)
    if(para$phenotype){
      inimean_vc[i] <- Phenotype(density4, para)
      inivar_vc[i] <- PhVariance(density4, para)
    }
    inifit[i] <- Fitness(density4, para)
    
    # Comparison between this cycle and previous one:
    c.change <- 0.00001                             # criteria: 0.001% of the previous value
    
    s1 <- (step.change(inipop, i) < c.change)
    s2 <- (step.change(inimean, i) < c.change)
    s3 <- (step.change(inivar, i) < c.change)
    if(para$phenotype){
      s4 <- (step.change(inimean_vc, i) < c.change)
      s5 <- (step.change(inivar_vc, i) < c.change)
    }else{
      s4 <- TRUE
      s5 <- TRUE
    }
    s6 <- all(abs(density4 - density) < c.change * density)
    
    # update:
    i <- i+1
    stable <- all(s1, s2, s3, s4, s5, s6)
    density <- density4
  }
  
  print(paste("The population reach stable in", i-1, "generations."))
  
  return(list(init=1:(i-1), density=density, 
              inipop=inipop, inimean=inimean, inivar=inivar, inisum=inipop * inimean,
              inimean_vc=inimean_vc, inivar_vc=inivar_vc, inifit=inifit, generation=i-1))
}


# c) Selection happen both before and after migration:
Initial_both <- function(density0, para, alter){
  # Input: - density0: initial genotype distribution of the field population
  #        - para: parameter list
  #        - alter: indicator about which reproduction assumption is used
  
  # Output: (equilibrium) genotype distribution of the field population before releasing
  
  # Initialization
  Sz <- para$Sz
  density <- density0
  
  # Changes over generations
  init <- 200                 # Maximum number of generations to run for the field population to reach equilibrium
  
  inipop <- rep(0,(init+1))
  inimean <- rep(0,(init+1))
  inivar <- rep(0,(init+1))
  inimean_vc <- rep(0,(init+1))
  inivar_vc <- rep(0,(init+1))
  inifit <- rep(0,(init+1))
  
  # The initial status:
  inipop[1] <- Popsize(Sz, density)
  inimean[1] <- Genotype(density, para)
  inivar[1] <- GVariance(density, para)
  if(para$phenotype){
    inimean_vc[1] <- Phenotype(density, para)
    inivar_vc[1] <- PhVariance(density, para)
  }
  inifit[1] <- Fitness(density, para)
  
  # start value of counter
  i <- 2
  stable <- FALSE
  
  while (i <= (init+1) & stable == FALSE){
    if(para$mig==1){
      density <- Migration(density, para)                         # Immigration
    }
    density1 <- Selection(density, para)                          # Stabilizing selection after migration
    
    if (alter==1){
      density2 <- reproduction1(density1, para)   
    }else if (alter==2){
      densitym <- (para$sex/(1+para$sex))*density1
      densityf <- (1/(1+para$sex))*density1
      density2 <- reproduction2(densitym, densityf, para)
    }                                                              # Reproduction (including density independence of egg hatching)
    
    density3 <- Den.dependent(density2, para)                      # Density dependent survival
    density4 <- Den.independent(density3, para)                    # Density independent survival
    density5 <- Selection(density4, para)                          # Selection before migration (second selection)
    
    inipop[i] <- Popsize(Sz, density5)
    inimean[i] <- Genotype(density5, para)
    inivar[i] <- GVariance(density5, para)
    if(para$phenotype){
      inimean_vc[i] <- Phenotype(density5, para)
      inivar_vc[i] <- PhVariance(density5, para)
    }
    inifit[i] <- Fitness(density5, para)
    
    # Comparison between this cycle and previous one:
    c.change <- 0.00001                             # criteria: 0.01% of the previous value
    
    s1 <- (step.change(inipop, i) < c.change)
    s2 <- (step.change(inimean, i) < c.change)
    s3 <- (step.change(inivar, i) < c.change)
    if(para$phenotype){
      s4 <- (step.change(inimean_vc, i) < c.change)
      s5 <- (step.change(inivar_vc, i) < c.change)
    }else{
      s4 <- TRUE
      s5 <- TRUE
    }
    s6 <- all(abs(density5 - density) < c.change * density)
    
    # update:
    i <- i+1
    stable <- all(s1, s2, s3, s4, s5, s6)
    density <- density5
  }
  
  print(paste("The population reach stable in", i-1, "generations."))
  
  return(list(init=1:(i-1), density=density, 
              inipop=inipop, inimean=inimean, inivar=inivar, inisum=inipop * inimean,
              inimean_vc=inimean_vc, inivar_vc=inivar_vc, inifit=inifit, generation=i-1))
}


step.change <- function(v, c){
  return(abs((v[c] - v[c-1]) / v[c]))
}


## Iterate the life cycle with releasing (main function!!):
#-----#
# a) Release both sex, selection happens only before migration and releasing
Iterate_before_both <- function(equaldensity, para, alter){
  # Input: - equaldensity: genotype distribution of the field population before releasing
  #        - para: parameter list
  #        - alter: indicator about which reproduction assumption is used
  
  # Output: genotype distribution of the field population after "gen" generations of releasing
  
  Sz <- para$Sz
  vz <- para$vz
  density <- equaldensity
  
  # Create vectors to restore the population characteristic for each generations
  postpop <- rep(0,(para$gen+1))
  postmean <- rep(0,(para$gen+1))
  postvar <- rep(0,(para$gen+1))
  postmean_vc <- rep(0,(para$gen+1))
  postvar_vc <- rep(0,(para$gen+1))
  postfit <- rep(0, (para$gen+1))
  
  # Input the values for the fiest generation
  postpop[1] <- Popsize(Sz, density)
  postmean[1] <- Genotype(density, para)
  postvar[1] <- GVariance(density, para)
  if(para$phenotype){
    postmean_vc[1] <- Phenotype(density, para)
    postvar_vc[1] <- PhVariance(density, para)
  }
  postfit[1] <- Fitness(density, para)
  
  # Releasing for "gen" generations
  for (i in 2:(para$gen+1)){
    density1 <- Releasing1(density, para)    # Releasing both sex, generate a combined distribution
    
    for(f in 1:para$fq){
      if(para$mig==1){
        density1 <- Migration(density1, para)                   # Immigration
      }
      
      if (alter==1){
        density2 <- reproduction1(density1, para)   
      }else if (alter==2){
        density_male <- (para$sex/(1+para$sex))*density1
        density_female <- (1/(1+para$sex))*density1
        density2 <- reproduction2(density_male, density_female, para)
      }                                                         # Releasing and reproduction
      
      density3 <- Den.dependent(density2, para)                 # Density dependent survival
      density4 <- Den.independent(density3, para)               # Density independent survival  
      density5 <- Selection(density4, para)                     # Stabilizing selection
      
      density1 <- density5                                      # Next cycle
    }
    
    # Assign the population from the last time step of the current generation as the starting point of the next generation
    density <- density5     ## Next cycle
    
    # Record outcomes of the releasing
    postpop[i] <- Popsize(Sz, density)
    postmean[i] <- Genotype(density, para)
    postvar[i] <- GVariance(density, para)
    postfit[i] <- Fitness(density, para)
    if(para$phenotype){
      postmean_vc[i] <- Phenotype(density, para)
      postvar_vc[i] <- PhVariance(density, para)
    }
    
  }
  return (list(time=0:para$gen, density=density, 
               postpop=postpop, postmean=postmean, postvar=postvar, postfit=postfit, postsum = postmean * postpop,
               postmean_vc=postmean_vc, postvar_vc=postvar_vc))
}


# b) Release only male, selection happens only before migration and releasing
# When only males are released, alter can only be 2 (use matrix calculation, instead of FFT)
Iterate_before_male <- function(equaldensity, para, alter=2){
  # Input: - equaldensity: genotype distribution of the field population before releasing
  #        - para: parameter list
  #        - alter: indicator about which reproduction assumption is used
  
  # Output: genotype distribution of the field population after "gen" generations of releasing
  
  Sz <- para$Sz
  vz <- para$vz
  density <- equaldensity
  
  # Create vectors to restore the population characteristic for each generations
  postpop <- rep(0,(para$gen+1))
  postmean <- rep(0,(para$gen+1))
  postvar <- rep(0,(para$gen+1))
  postmean_vc <- rep(0,(para$gen+1))
  postvar_vc <- rep(0,(para$gen+1))
  postfit <- rep(0, (para$gen+1))
  
  # Input the values for the fiest generation
  postpop[1] <- Popsize(Sz, density)
  postmean[1] <- Genotype(density, para)
  postvar[1] <- GVariance(density, para)
  if(para$phenotype){
    postmean_vc[1] <- Phenotype(density, para)
    postvar_vc[1] <- PhVariance(density, para)
  }
  postfit[1] <- Fitness(density, para)
  
  # Releasing only males for "gen" generations
  for (i in 2:(para$gen+1)){
    
    if(para$mig==1){
      density <- Migration(density, para)                     # Immigration
    }
    density1 <- Releasing2(density, para)                     # Releasing only male: two distributions (male & female)
    
    density_male <- density1$male
    density_female <- density1$female
    density2 <- reproduction2(density_male, density_female, para)   # Reproduction
    
    density3 <- Den.dependent(density2, para)                 # Density dependent survival
    density4 <- Den.independent(density3, para)               # Density independent survival  
    density5 <- Selection(density4, para)                     # Stabilizing selection
    
    # Assign the population from the last time step of the current generation as the starting point of the next generation
    density <- density5                                       # Next cycle
    
    if(para$fq>=2){           # If not releasing every generation
      for(f in 1:(para$fq-1)){
        if(para$mig==1){
          density <- Migration(density, para)                     # Immigration
        }
        
        density_male <- (para$sex/(1+para$sex))*density
        density_female <- (1/(1+para$sex))*density
        density1 <- reproduction2(density_male, density_female, para)
        
        density2 <- Den.dependent(density1, para)                 # Density dependent survival
        density3 <- Den.independent(density2, para)               # Density independent survival  
        density4 <- Selection(density3, para)                     # Stabilizing selection
        
        density <- density4                                       # Next generation
      }
    }
    
    # Record the outcome of releasing
    postpop[i] <- Popsize(Sz, density)
    postmean[i] <- Genotype(density, para)
    postvar[i] <- GVariance(density, para)
    postfit[i] <- Fitness(density, para)
    if(para$phenotype){
      postmean_vc[i] <- Phenotype(density, para)
      postvar_vc[i] <- PhVariance(density, para)
    }
    
  }
  return (list(time=0:para$gen, density=density, 
               postpop=postpop, postmean=postmean, postvar=postvar, postfit=postfit, postsum = postmean * postpop,
               postmean_vc=postmean_vc, postvar_vc=postvar_vc))
}




#-----#
# c) Release both sex, selection happens only after migration and releasing
Iterate_after_both <- function(equaldensity, para, alter){
  # Input: - equaldensity: genotype distribution of the field population before releasing
  #        - para: parameter list
  #        - alter: indicator about which reproduction assumption is used
  
  # Output: genotype distribution of the field population after "gen" generations of releasing
  
  Sz <- para$Sz
  vz <- para$vz
  density <- equaldensity
  
  # Create vectors to restore the population characteristic for each generations
  postpop <- rep(0,(para$gen+1))
  postmean <- rep(0,(para$gen+1))
  postvar <- rep(0,(para$gen+1))
  postmean_vc <- rep(0,(para$gen+1))
  postvar_vc <- rep(0,(para$gen+1))
  postfit <- rep(0, (para$gen+1))
  
  # Input the values for the fiest generation
  postpop[1] <- Popsize(Sz, density)
  postmean[1] <- Genotype(density, para)
  postvar[1] <- GVariance(density, para)
  if(para$phenotype){
    postmean_vc[1] <- Phenotype(density, para)
    postvar_vc[1] <- PhVariance(density, para)
  }
  postfit[1] <- Fitness(density, para)
  
  # Releasing for "gen" generations
  for (i in 2:(para$gen+1)){
    density1 <- Releasing1(density, para)    # Releasing both sex, generate a combined distribution
    
    for(f in 1:para$fq){
      if(para$mig==1){
        density1 <- Migration(density1, para)                   # Immigration
      }
      density2 <- Selection(density1, para)                     # Stabilizing selection
      
      if (alter==1){
        density3 <- reproduction1(density2, para)   
      }else if (alter==2){
        density_male <- (para$sex/(1+para$sex))*density2
        density_female <- (1/(1+para$sex))*density2
        density3 <- reproduction2(density_male, density_female, para)
      }                                                         # Releasing and reproduction
      
      density4 <- Den.dependent(density3, para)                 # Density dependent survival
      density5 <- Den.independent(density4, para)               # Density independent survival  
      
      density1 <- density5                                      # Next cycle
    }
    
    # Assign the population from the last time step of the current generation as the starting point of the next generation
    density <- density5     ## Next cycle
    
    # Record outcomes of the releasing
    postpop[i] <- Popsize(Sz, density)
    postmean[i] <- Genotype(density, para)
    postvar[i] <- GVariance(density, para)
    postfit[i] <- Fitness(density, para)
    if(para$phenotype){
      postmean_vc[i] <- Phenotype(density, para)
      postvar_vc[i] <- PhVariance(density, para)
    }
    
  }
  return (list(time=0:para$gen, density=density, 
               postpop=postpop, postmean=postmean, postvar=postvar, postfit=postfit, postsum = postmean * postpop,
               postmean_vc=postmean_vc, postvar_vc=postvar_vc))
}


# d) Release only male, selection happens only after migration and releasing
# When only males are released, alter can only be 2 (use matrix calculation, instead of FFT)
Iterate_after_male <- function(equaldensity, para, alter=2){
  # Input: - equaldensity: genotype distribution of the field population before releasing
  #        - para: parameter list
  #        - alter: indicator about which reproduction assumption is used
  
  # Output: genotype distribution of the field population after "gen" generations of releasing
  
  Sz <- para$Sz
  vz <- para$vz
  density <- equaldensity
  
  # Create vectors to restore the population characteristic for each generations
  postpop <- rep(0,(para$gen+1))
  postmean <- rep(0,(para$gen+1))
  postvar <- rep(0,(para$gen+1))
  postmean_vc <- rep(0,(para$gen+1))
  postvar_vc <- rep(0,(para$gen+1))
  postfit <- rep(0, (para$gen+1))
  
  # Input the values for the fiest generation
  postpop[1] <- Popsize(Sz, density)
  postmean[1] <- Genotype(density, para)
  postvar[1] <- GVariance(density, para)
  if(para$phenotype){
    postmean_vc[1] <- Phenotype(density, para)
    postvar_vc[1] <- PhVariance(density, para)
  }
  postfit[1] <- Fitness(density, para)
  
  # Releasing only males for "gen" generations
  for (i in 2:(para$gen+1)){
    
    if(para$mig==1){
      density <- Migration(density, para)                     # Immigration
    }
    density1 <- Releasing2(density, para)                     # Releasing only male: two distributions (male & female)
    
    density_male <- Selection(density1$male, para)
    density_female <- Selection(density1$female, para)
    density3 <- reproduction2(density_male, density_female, para)   # Reproduction
    
    density4 <- Den.dependent(density3, para)                 # Density dependent survival
    density5 <- Den.independent(density4, para)               # Density independent survival
    
    # Assign the population from the last time step of the current generation as the starting point of the next generation
    density <- density5                                       # Next cycle
    
    if(para$fq>=2){           # If not releasing every generation
      for(f in 1:(para$fq-1)){
        if(para$mig==1){
          density <- Migration(density, para)                     # Immigration
        }
        density1 <- Selection(density, para)                            # Stabilizing selection
        
        density_male <- (para$sex/(1+para$sex))*density1
        density_female <- (1/(1+para$sex))*density1
        density2 <- reproduction2(density_male, density_female, para)
        
        density3 <- Den.dependent(density2, para)                 # Density dependent survival
        density4 <- Den.independent(density3, para)               # Density independent survival  
        
        density <- density4                                       # Next generation
      }
    }
    
    # Record the outcome of releasing
    postpop[i] <- Popsize(Sz, density)
    postmean[i] <- Genotype(density, para)
    postvar[i] <- GVariance(density, para)
    postfit[i] <- Fitness(density, para)
    if(para$phenotype){
      postmean_vc[i] <- Phenotype(density, para)
      postvar_vc[i] <- PhVariance(density, para)
    }
    
  }
  return (list(time=0:para$gen, density=density, 
               postpop=postpop, postmean=postmean, postvar=postvar, postfit=postfit, postsum = postmean * postpop,
               postmean_vc=postmean_vc, postvar_vc=postvar_vc))
}


#-----#
# e) Release both sex, selection happens both before and after migration and releasing
Iterate_both_both <- function(equaldensity, para, alter){
  # Input: - equaldensity: genotype distribution of the field population before releasing
  #        - para: parameter list
  #        - alter: indicator about which reproduction assumption is used
  
  # Output: genotype distribution of the field population after "gen" generations of releasing
  
  Sz <- para$Sz
  vz <- para$vz
  density <- equaldensity
  
  # Create vectors to restore the population characteristic for each generations
  postpop <- rep(0,(para$gen+1))
  postmean <- rep(0,(para$gen+1))
  postvar <- rep(0,(para$gen+1))
  postmean_vc <- rep(0,(para$gen+1))
  postvar_vc <- rep(0,(para$gen+1))
  postfit <- rep(0, (para$gen+1))
  
  # Input the values for the fiest generation
  postpop[1] <- Popsize(Sz, density)
  postmean[1] <- Genotype(density, para)
  postvar[1] <- GVariance(density, para)
  if(para$phenotype){
    postmean_vc[1] <- Phenotype(density, para)
    postvar_vc[1] <- PhVariance(density, para)
  }
  postfit[1] <- Fitness(density, para)
  
  # Releasing for "gen" generations
  for (i in 2:(para$gen+1)){
    density1 <- Releasing1(density, para)    # Releasing both sex, generate a combined distribution
    
    for(f in 1:para$fq){
      if(para$mig==1){
        density1 <- Migration(density1, para)                   # Immigration
      }
      density2 <- Selection(density1, para)                     # Stabilizing selection
      
      if (alter==1){
        density3 <- reproduction1(density2, para)   
      }else if (alter==2){
        density_male <- (para$sex/(1+para$sex))*density2
        density_female <- (1/(1+para$sex))*density2
        density3 <- reproduction2(density_male, density_female, para)
      }                                                         # Releasing and reproduction
      
      density4 <- Den.dependent(density3, para)                 # Density dependent survival
      density5 <- Den.independent(density4, para)               # Density independent survival  
      density6 <- Selection(density5, para)                     # Stabilizing selection (again)
      
      density1 <- density6                                      # Next cycle
    }
    
    # Assign the population from the last time step of the current generation as the starting point of the next generation
    density <- density6     ## Next cycle
    
    # Record outcomes of the releasing
    postpop[i] <- Popsize(Sz, density)
    postmean[i] <- Genotype(density, para)
    postvar[i] <- GVariance(density, para)
    postfit[i] <- Fitness(density, para)
    if(para$phenotype){
      postmean_vc[i] <- Phenotype(density, para)
      postvar_vc[i] <- PhVariance(density, para)
    }
    
  }
  return (list(time=0:para$gen, density=density, 
               postpop=postpop, postmean=postmean, postvar=postvar, postfit=postfit, postsum = postmean * postpop,
               postmean_vc=postmean_vc, postvar_vc=postvar_vc))
}


# f) Release only male, selection happens both before and after migration and releasing
# When only males are released, alter can only be 2 (use matrix calculation, instead of FFT)
Iterate_both_male <- function(equaldensity, para, alter=2){
  # Input: - equaldensity: genotype distribution of the field population before releasing
  #        - para: parameter list
  #        - alter: indicator about which reproduction assumption is used
  
  # Output: genotype distribution of the field population after "gen" generations of releasing
  
  Sz <- para$Sz
  vz <- para$vz
  density <- equaldensity
  
  # Create vectors to restore the population characteristic for each generations
  postpop <- rep(0,(para$gen+1))
  postmean <- rep(0,(para$gen+1))
  postvar <- rep(0,(para$gen+1))
  postmean_vc <- rep(0,(para$gen+1))
  postvar_vc <- rep(0,(para$gen+1))
  postfit <- rep(0, (para$gen+1))
  
  # Input the values for the fiest generation
  postpop[1] <- Popsize(Sz, density)
  postmean[1] <- Genotype(density, para)
  postvar[1] <- GVariance(density, para)
  if(para$phenotype){
    postmean_vc[1] <- Phenotype(density, para)
    postvar_vc[1] <- PhVariance(density, para)
  }
  postfit[1] <- Fitness(density, para)
  
  # Releasing only males for "gen" generations
  for (i in 2:(para$gen+1)){
    
    if(para$mig==1){
      density <- Migration(density, para)                     # Immigration
    }
    density1 <- Releasing2(density, para)                     # Releasing only male: two distributions (male & female)
    
    density_male <- Selection(density1$male, para)
    density_female <- Selection(density1$female, para)
    density3 <- reproduction2(density_male, density_female, para)   # Reproduction
    
    density4 <- Den.dependent(density3, para)                 # Density dependent survival
    density5 <- Den.independent(density4, para)               # Density independent survival
    density6 <- Selection(density5, para)                     # Stabilizing selection (again)
    
    # Assign the population from the last time step of the current generation as the starting point of the next generation
    density <- density6                                       # Next cycle
    
    if(para$fq>=2){           # If not releasing every generation
      for(f in 1:(para$fq-1)){
        if(para$mig==1){
          density <- Migration(density, para)                     # Immigration
        }
        density1 <- Selection(density, para)                            # Stabilizing selection
        
        density_male <- (para$sex/(1+para$sex))*density1
        density_female <- (1/(1+para$sex))*density1
        density2 <- reproduction2(density_male, density_female, para)
        
        density3 <- Den.dependent(density2, para)                 # Density dependent survival
        density4 <- Den.independent(density3, para)               # Density independent survival  
        density5 <- Selection(density4, para)                     # Stabilizing selection (again)
        
        density <- density5                                       # Next generation
      }
    }
    
    # Record the outcome of releasing
    postpop[i] <- Popsize(Sz, density)
    postmean[i] <- Genotype(density, para)
    postvar[i] <- GVariance(density, para)
    postfit[i] <- Fitness(density, para)
    if(para$phenotype){
      postmean_vc[i] <- Phenotype(density, para)
      postvar_vc[i] <- PhVariance(density, para)
    }
    
  }
  return (list(time=0:para$gen, density=density, 
               postpop=postpop, postmean=postmean, postvar=postvar, postfit=postfit, postsum = postmean * postpop,
               postmean_vc=postmean_vc, postvar_vc=postvar_vc))
}


# g) Release males and blood-fed females, selection happens only before migration and releasing
# Males enter the cycle before reproduction while females enter after reproduction (equals to releasing eggs)
# When only males are released, alter can only be 2 (use matrix calculation, instead of FFT)
Iterate_before_bloodfed <- function(equaldensity, para, alter=2){
  # Input: - equaldensity: genotype distribution of the field population before releasing
  #        - para: parameter list
  #        - alter: indicator about which reproduction assumption is used
  
  # Output: genotype distribution of the field population after "gen" generations of releasing
  
  Sz <- para$Sz
  vz <- para$vz
  density <- equaldensity
  
  # Create vectors to restore the population characteristic for each generations
  postpop <- rep(0,(para$gen+1))
  postmean <- rep(0,(para$gen+1))
  postvar <- rep(0,(para$gen+1))
  postmean_vc <- rep(0,(para$gen+1))
  postvar_vc <- rep(0,(para$gen+1))
  postfit <- rep(0, (para$gen+1))
  
  # Input the values for the fiest generation
  postpop[1] <- Popsize(Sz, density)
  postmean[1] <- Genotype(density, para)
  postvar[1] <- GVariance(density, para)
  if(para$phenotype){
    postmean_vc[1] <- Phenotype(density, para)
    postvar_vc[1] <- PhVariance(density, para)
  }
  postfit[1] <- Fitness(density, para)
  
  # Change the size of the releasing population as males and females will enter the cycle in different places
  para.male <- para
  para.female <- para
  
  para.male$rel <- 0.5 * para$rel
  para.female$rel <- 0.5 * para$rel * para$R.rel
  
  # Releasing only males for "gen" generations
  for (i in 2:(para$gen+1)){
    
    if(para$mig==1){
      density <- Migration(density, para)                     # Immigration
    }
    density1 <- Releasing2(density, para.male)                # Releasing only male: two distributions (male & female)
    
    density_male <- density1$male
    density_female <- density1$female
    density2 <- reproduction2(density_male, density_female, para)   # Reproduction
    
    density2.2 <- Releasing1(density2, para.female)           # Releasing blood-fed females: equal to releasing eggs
    
    density3 <- Den.dependent(density2.2, para)                 # Density dependent survival
    density4 <- Den.independent(density3, para)               # Density independent survival  
    density5 <- Selection(density4, para)                     # Stabilizing selection
    
    # Assign the population from the last time step of the current generation as the starting point of the next generation
    density <- density5                                       # Next cycle
    
    if(para$fq>=2){           # If not releasing every generation
      for(f in 1:(para$fq-1)){
        if(para$mig==1){
          density <- Migration(density, para)                     # Immigration
        }
        
        density_male <- (para$sex/(1+para$sex))*density
        density_female <- (1/(1+para$sex))*density
        density1 <- reproduction2(density_male, density_female, para)
        
        density2 <- Den.dependent(density1, para)                 # Density dependent survival
        density3 <- Den.independent(density2, para)               # Density independent survival  
        density4 <- Selection(density3, para)                     # Stabilizing selection
        
        density <- density4                                       # Next generation
      }
    }
    
    # Record the outcome of releasing
    postpop[i] <- Popsize(Sz, density)
    postmean[i] <- Genotype(density, para)
    postvar[i] <- GVariance(density, para)
    postfit[i] <- Fitness(density, para)
    if(para$phenotype){
      postmean_vc[i] <- Phenotype(density, para)
      postvar_vc[i] <- PhVariance(density, para)
    }
    
  }
  return (list(time=0:para$gen, density=density, 
               postpop=postpop, postmean=postmean, postvar=postvar, postfit=postfit, postsum = postmean * postpop,
               postmean_vc=postmean_vc, postvar_vc=postvar_vc))
}



# h) Release males and blood-fed females, selection happens only after migration and releasing
# Males enter the cycle before reproduction while females enter after reproduction (equals to releasing eggs)
# When only males are released, alter can only be 2 (use matrix calculation, instead of FFT)
Iterate_after_bloodfed <- function(equaldensity, para, alter=2){
  # Input: - equaldensity: genotype distribution of the field population before releasing
  #        - para: parameter list
  #        - alter: indicator about which reproduction assumption is used
  
  # Output: genotype distribution of the field population after "gen" generations of releasing
  
  Sz <- para$Sz
  vz <- para$vz
  density <- equaldensity
  
  # Create vectors to restore the population characteristic for each generations
  postpop <- rep(0,(para$gen+1))
  postmean <- rep(0,(para$gen+1))
  postvar <- rep(0,(para$gen+1))
  postmean_vc <- rep(0,(para$gen+1))
  postvar_vc <- rep(0,(para$gen+1))
  postfit <- rep(0, (para$gen+1))
  
  # Input the values for the fiest generation
  postpop[1] <- Popsize(Sz, density)
  postmean[1] <- Genotype(density, para)
  postvar[1] <- GVariance(density, para)
  if(para$phenotype){
    postmean_vc[1] <- Phenotype(density, para)
    postvar_vc[1] <- PhVariance(density, para)
  }
  postfit[1] <- Fitness(density, para)
  
  # Change the size of the releasing population as males and females will enter the cycle in different places
  para.male <- para
  para.female <- para
  
  para.male$rel <- 0.5 * para$rel
  para.female$rel <- 0.5 * para$rel * para$R.rel
  
  # Releasing only males for "gen" generations
  for (i in 2:(para$gen+1)){
    
    if(para$mig==1){
      density <- Migration(density, para)                     # Immigration
    }
    density1 <- Releasing2(density, para.male)                     # Releasing only male: two distributions (male & female)
    
    density_male <- Selection(density1$male, para)
    density_female <- Selection(density1$female, para)
    density3 <- reproduction2(density_male, density_female, para)   # Reproduction
    
    density3.2 <- Releasing1(density3, para.female)           # Releasing blood-fed females: equal to releasing eggs
    
    density4 <- Den.dependent(density3.2, para)               # Density dependent survival
    density5 <- Den.independent(density4, para)               # Density independent survival
    
    # Assign the population from the last time step of the current generation as the starting point of the next generation
    density <- density5                                       # Next cycle
    
    if(para$fq>=2){           # If not releasing every generation
      for(f in 1:(para$fq-1)){
        if(para$mig==1){
          density <- Migration(density, para)                     # Immigration
        }
        density1 <- Selection(density, para)                            # Stabilizing selection
        
        density_male <- (para$sex/(1+para$sex))*density1
        density_female <- (1/(1+para$sex))*density1
        density2 <- reproduction2(density_male, density_female, para)
        
        density3 <- Den.dependent(density2, para)                 # Density dependent survival
        density4 <- Den.independent(density3, para)               # Density independent survival  
        
        density <- density4                                       # Next generation
      }
    }
    
    # Record the outcome of releasing
    postpop[i] <- Popsize(Sz, density)
    postmean[i] <- Genotype(density, para)
    postvar[i] <- GVariance(density, para)
    postfit[i] <- Fitness(density, para)
    if(para$phenotype){
      postmean_vc[i] <- Phenotype(density, para)
      postvar_vc[i] <- PhVariance(density, para)
    }
    
  }
  return (list(time=0:para$gen, density=density, 
               postpop=postpop, postmean=postmean, postvar=postvar, postfit=postfit, postsum = postmean * postpop,
               postmean_vc=postmean_vc, postvar_vc=postvar_vc))
}




# i) Release males and blood-fed females, selection happens both before and after migration and releasing
# Males enter the cycle before reproduction while females enter after reproduction (equals to releasing eggs)
# When only males are released, alter can only be 2 (use matrix calculation, instead of FFT)
Iterate_both_bloodfed <- function(equaldensity, para, alter=2){
  # Input: - equaldensity: genotype distribution of the field population before releasing
  #        - para: parameter list
  #        - alter: indicator about which reproduction assumption is used
  
  # Output: genotype distribution of the field population after "gen" generations of releasing
  
  Sz <- para$Sz
  vz <- para$vz
  density <- equaldensity
  
  # Create vectors to restore the population characteristic for each generations
  postpop <- rep(0,(para$gen+1))
  postmean <- rep(0,(para$gen+1))
  postvar <- rep(0,(para$gen+1))
  postmean_vc <- rep(0,(para$gen+1))
  postvar_vc <- rep(0,(para$gen+1))
  postfit <- rep(0, (para$gen+1))
  
  # Input the values for the fiest generation
  postpop[1] <- Popsize(Sz, density)
  postmean[1] <- Genotype(density, para)
  postvar[1] <- GVariance(density, para)
  if(para$phenotype){
    postmean_vc[1] <- Phenotype(density, para)
    postvar_vc[1] <- PhVariance(density, para)
  }
  postfit[1] <- Fitness(density, para)
  
  # Change the size of the releasing population as males and females will enter the cycle in different places
  para.male <- para
  para.female <- para
  
  para.male$rel <- 0.5 * para$rel
  para.female$rel <- 0.5 * para$rel * para$R.rel
  
  # Releasing only males for "gen" generations
  for (i in 2:(para$gen+1)){
    
    if(para$mig==1){
      density <- Migration(density, para)                     # Immigration
    }
    density1 <- Releasing2(density, para.male)                     # Releasing only male: two distributions (male & female)
    
    density_male <- Selection(density1$male, para)
    density_female <- Selection(density1$female, para)
    density3 <- reproduction2(density_male, density_female, para)   # Reproduction
    
    density3.2 <- Releasing1(density3, para.female)           # Releasing blood-fed females: equal to releasing eggs
    
    density4 <- Den.dependent(density3.2, para)               # Density dependent survival
    density5 <- Den.independent(density4, para)               # Density independent survival
    density6 <- Selection(density5, para)                     # Stabilizing selection (again)
    
    # Assign the population from the last time step of the current generation as the starting point of the next generation
    density <- density6                                       # Next cycle
    
    if(para$fq>=2){           # If not releasing every generation
      for(f in 1:(para$fq-1)){
        if(para$mig==1){
          density <- Migration(density, para)                     # Immigration
        }
        density1 <- Selection(density, para)                            # Stabilizing selection
        
        density_male <- (para$sex/(1+para$sex))*density1
        density_female <- (1/(1+para$sex))*density1
        density2 <- reproduction2(density_male, density_female, para)
        
        density3 <- Den.dependent(density2, para)                 # Density dependent survival
        density4 <- Den.independent(density3, para)               # Density independent survival  
        density5 <- Selection(density4, para)                     # Stabilizing selection (again)
        
        density <- density5                                       # Next generation
      }
    }
    
    # Record the outcome of releasing
    postpop[i] <- Popsize(Sz, density)
    postmean[i] <- Genotype(density, para)
    postvar[i] <- GVariance(density, para)
    postfit[i] <- Fitness(density, para)
    if(para$phenotype){
      postmean_vc[i] <- Phenotype(density, para)
      postvar_vc[i] <- PhVariance(density, para)
    }
    
  }
  return (list(time=0:para$gen, density=density, 
               postpop=postpop, postmean=postmean, postvar=postvar, postfit=postfit, postsum = postmean * postpop,
               postmean_vc=postmean_vc, postvar_vc=postvar_vc))
}




### Function for relapse: medified from previous "Initial" function

# a) Selection happen only before migration:
Relapse_before <- function(post_releasing_density, para, alter, g=20){
  # Input: - density0: initial genotype distribution of the field population
  #        - para: parameter list
  #        - alter: indicator about which reproduction assumption is used
  
  # Output: (equilibrium) genotype distribution of the field population before releasing
  
  # Initialization
  Sz <- para$Sz
  density <- post_releasing_density
  
  # Changes over generations
  init <- g                      # How many generations of relapse
  
  # Initialize result vectors
  inipop <- rep(0,(init+1))
  inimean <- rep(0,(init+1))
  inivar <- rep(0,(init+1))
  inimean_vc <- rep(0,(init+1))
  inivar_vc <- rep(0,(init+1))
  inifit <- rep(0,(init+1))
  
  # The initial status:
  inipop[1] <- Popsize(Sz, density)
  inimean[1] <- Genotype(density, para)
  inivar[1] <- GVariance(density, para)
  if(para$phenotype){
    inimean_vc[1] <- Phenotype(density, para)
    inivar_vc[1] <- PhVariance(density, para)
  }
  inifit[1] <- Fitness(density, para)
  
  for (i in 2:(init+1)){
    if(para$mig==1){
      density <- Migration(density, para)    # Immigration
    }
    
    if (alter==1){
      density1 <- reproduction1(density, para)   
    }else if (alter==2){
      densitym <- (para$sex/(1+para$sex))*density
      densityf <- (1/(1+para$sex))*density
      density1 <- reproduction2(densitym, densityf, para)
    }   # Reproduction (including density independence of egg hatching)
    
    density2 <- Den.dependent(density1, para)                      # Density dependent survival
    density3 <- Den.independent(density2, para)                    # Density independent survival
    density4 <- Selection(density3, para)                          # Stabilizing selection
    
    inipop[i] <- Popsize(Sz, density4)
    inimean[i] <- Genotype(density4, para)
    inivar[i] <- GVariance(density4, para)
    if(para$phenotype){
      inimean_vc[i] <- Phenotype(density4, para)
      inivar_vc[i] <- PhVariance(density4, para)
    }
    inifit[i] <- Fitness(density4, para)
    density <- density4     # Next cycle
  }
  return(list(init=0:init, density=density, 
              inipop=inipop, inimean=inimean, inivar=inivar, inisum = inipop * inimean,
              inimean_vc=inimean_vc, inivar_vc=inivar_vc, inifit=inifit))
}


# b) Selection happen only after migration:
Relapse_after <- function(post_releasing_density, para, alter, g=20){
  # Input: - density0: initial genotype distribution of the field population
  #        - para: parameter list
  #        - alter: indicator about which reproduction assumption is used
  
  # Output: (equilibrium) genotype distribution of the field population before releasing
  
  # Initialization
  Sz <- para$Sz
  density <- post_releasing_density
  
  # Changes over generations
  init <- g                      # How many generations of relapse
  
  # Initialize result vectors
  inipop <- rep(0,(init+1))
  inimean <- rep(0,(init+1))
  inivar <- rep(0,(init+1))
  inimean_vc <- rep(0,(init+1))
  inivar_vc <- rep(0,(init+1))
  inifit <- rep(0,(init+1))
  
  # The initial status:
  inipop[1] <- Popsize(Sz, density)
  inimean[1] <- Genotype(density, para)
  inivar[1] <- GVariance(density, para)
  if(para$phenotype){
    inimean_vc[1] <- Phenotype(density, para)
    inivar_vc[1] <- PhVariance(density, para)
  }
  inifit[1] <- Fitness(density, para)
  
  for (i in 2:(init+1)){
    if(para$mig==1){
      density <- Migration(density, para)                         # Immigration
    }
    density1 <- Selection(density, para)                          # Stabilizing selection
    
    if (alter==1){
      density2 <- reproduction1(density1, para)   
    }else if (alter==2){
      densitym <- (para$sex/(1+para$sex))*density1
      densityf <- (1/(1+para$sex))*density1
      density2 <- reproduction2(densitym, densityf, para)
    }                                 # Reproduction (including density independence of egg hatching)
    
    density3 <- Den.dependent(density2, para)                      # Density dependent survival
    density4 <- Den.independent(density3, para)                    # Density independent survival
    
    inipop[i] <- Popsize(Sz, density4)
    inimean[i] <- Genotype(density4, para)
    inivar[i] <- GVariance(density4, para)
    if(para$phenotype){
      inimean_vc[i] <- Phenotype(density4, para)
      inivar_vc[i] <- PhVariance(density4, para)
    }
    inifit[i] <- Fitness(density4, para)
    
    density <- density4     # Next cycle
  }
  return(list(init=0:init, density=density, 
              inipop=inipop, inimean=inimean, inivar=inivar, inisum = inipop * inimean,
              inimean_vc=inimean_vc, inivar_vc=inivar_vc, inifit=inifit))
}


# c) Selection happen both before and after migration:
Relapse_both <- function(post_releasing_density, para, alter, g=20){
  # Input: - density0: initial genotype distribution of the field population
  #        - para: parameter list
  #        - alter: indicator about which reproduction assumption is used
  
  # Output: (equilibrium) genotype distribution of the field population before releasing
  
  # Initialization
  Sz <- para$Sz
  density <- post_releasing_density
  
  # Changes over generations
  init <- g                      # How many generations of relapse
  
  # Initialize result vectors
  inipop <- rep(0,(init+1))
  inimean <- rep(0,(init+1))
  inivar <- rep(0,(init+1))
  inimean_vc <- rep(0,(init+1))
  inivar_vc <- rep(0,(init+1))
  inifit <- rep(0,(init+1))
  
  # The initial status:
  inipop[1] <- Popsize(Sz, density)
  inimean[1] <- Genotype(density, para)
  inivar[1] <- GVariance(density, para)
  if(para$phenotype){
    inimean_vc[1] <- Phenotype(density, para)
    inivar_vc[1] <- PhVariance(density, para)
  }
  inifit[1] <- Fitness(density, para)
  
  for (i in 2:(init+1)){
    if(para$mig==1){
      density <- Migration(density, para)                         # Immigration
    }
    density1 <- Selection(density, para)                          # Stabilizing selection
    
    if (alter==1){
      density2 <- reproduction1(density1, para)   
    }else if (alter==2){
      densitym <- (para$sex/(1+para$sex))*density1
      densityf <- (1/(1+para$sex))*density1
      density2 <- reproduction2(densitym, densityf, para)
    }                                 # Reproduction (including density independence of egg hatching)
    
    density3 <- Den.dependent(density2, para)                      # Density dependent survival
    density4 <- Den.independent(density3, para)                    # Density independent survival
    density5 <- Selection(density4, para)                          # Stablizing selection (again)
    
    inipop[i] <- Popsize(Sz, density5)
    inimean[i] <- Genotype(density5, para)
    inivar[i] <- GVariance(density5, para)
    if(para$phenotype){
      inimean_vc[i] <- Phenotype(density5, para)
      inivar_vc[i] <- PhVariance(density5, para)
    }
    inifit[i] <- Fitness(density5, para)
    
    density <- density5      # Next cycle
  }
  return(list(init=0:init, density=density, 
              inipop=inipop, inimean=inimean, inivar=inivar, inisum = inipop * inimean,
              inimean_vc=inimean_vc, inivar_vc=inivar_vc, inifit=inifit))
}


#-----# Summarizing functions
# Summarize a density: mean, var, population size, fitness
Densummary <- function (status, para){
  # Input:  - status: output from "Initial" or "Iterate"
  #         - para: parameter list
  
  # Output: a list containing the summary of the population genotype distribution
  
  d <- status$density
  m <- Genotype(d, para)
  v <- GVariance(d, para)
  p <- Popsize(para$Sz, d)
  f <- Fitness(d, para)
  inte <- m * p
  vc_m <- Phenotype(d, para)
  vc_var <- PhVariance(d, para)
  vc_sd <- sqrt(vc_var)
  vc_inte <- vc_m * p
  return (list(full=status, density=d, 
               dmean=m, dvar=v, dpop=p, dfit=f, dsum=inte, 
               vcmean=vc_m, vcvar=vc_var, vcsd=vc_sd, vcsum=vc_inte))
}


# Examine how frequent should the release being kept to maintain the outcome of the releasing
Maintain <- function(relap, out){
  # input: - relap: outcome of the "Initial" function used to calculate relapse
  #        - out: outcome of the "outcome" function
  #
  # Output: how many generations per each releasing are needed to maintain the outcome of intense releasing
  crit <- out$postmean_vc[length(out$postmean_vc)-1]
  return(max(which(relap$inimean_vc<=crit)))
}


# Compile three main functions to increase computing speed:
reproduction2 <- cmpfun(reproduction2)

Initial_before <- cmpfun(Initial_before)
Initial_after <- cmpfun(Initial_after)
Initial_both <- cmpfun(Initial_both)

Iterate_before_both <- cmpfun(Iterate_before_both)
Iterate_before_male <- cmpfun(Iterate_before_male)
Iterate_after_both <- cmpfun(Iterate_after_both)
Iterate_after_male <- cmpfun(Iterate_after_male)
Iterate_both_both <- cmpfun(Iterate_both_both)
Iterate_both_male <- cmpfun(Iterate_both_male)

Relapse_before <- cmpfun(Relapse_before)
Relapse_after <- cmpfun(Relapse_after)
Relapse_both <- cmpfun(Relapse_both)


#-----# Local sensitivity analysis:
# Modified the "para" variable used in the "Locals" function to implement local sensitivity analysis
Update <- function(newpara, sen, newvalue){
  # Input:  - newpara: the para that need modified
  #         - sen: the variable of interest
  #         - newvalue: the new value of the variable that need to change in "newpara"
  
  # Output: - a list containing the modified parameter set and a indicator of the variable of interest
  
  newpara[[which(names(newpara)==sen)]]=newvalue
  
  # Update components of variance if necessary:
  if (sen=="vp" | sen=="h"){
    newpara$vle <- newpara$vp*newpara$h
    newpara$ve  <- newpara$vp*(1-newpara$h)
    newpara$fsd <- sqrt(newpara$vle)
    newpara$Trans = dnorm(newpara$vz,mean=newpara$vz[newpara$Nz/2+1], sd=sqrt(newpara$vle/2))
    newpara$rotate = c((newpara$Nz/2+1):newpara$Nz, 1:(newpara$Nz/2))
    newpara$Trans[newpara$rotate] = newpara$Trans
  }
  return(newpara)
}


Locals <- function(sen, sen.range, para, alter, equaldensity=equ.density, default=outcome){
  
  # Input: - sen: the name of the variable
  #        - max, min, intv: the maximum, minimum and the step size of the variable (setting the range)
  #        - para: parameter list
  #        - alter: indicator about which reproduction assumption is used
  #        - density: initial genotype distribution of the field population
  #                    default = density0 (initial genotype distribution)
  #        - equaldensity: genotype distribution of the field population before releasing
  #                        default = equ.density (genotype distribution before releasing)
  #        - default: the result of "Itearate" function using the default parameters
  #                   default = outcome
  
  # Output: summary statistics (mean, var, pop size and "speed") changing with the variable
  
  l <- length(sen.range)
  
  relapse_gen <- 20                           # How many generations to keep track of after stop releasing
  gen_total <- para$gen + relapse_gen + 1     # Total number of generations to record
  
  sen_pop <- matrix(0, nrow = gen_total, ncol = l)      # Population size of all 31 generations under each variable value
  sen_mean <- matrix(0, nrow = gen_total, ncol = l)     # Population genotype mean of all 31 generations under each variable
  sen_var <- matrix(0, nrow = gen_total, ncol = l)      # Population genotype variance of all 31 generations under each variable
  sen_mean_vc <- matrix(0, nrow = gen_total, ncol = l)  # Population mean VC of all 31 generations under each variable
  sen_var_vc <- matrix(0, nrow = gen_total, ncol = l)   # Population variance of VC of all 31 generations under each variable
  sen_fit <- matrix(0, nrow = gen_total, ncol = l)      # Population mean fitness of all 31 generations under each variable value
  sen_speed <- rep(-1, l)                                  # The number of steps required to reach the "cre" 
  sen_integrateVC <- sen_pop * sen_mean                    # Integrated VC (genotype)
  
  parallel_results <- foreach(i = 1:l, .packages="truncnorm", .export = function_list) %dopar% {
    para.sen <- para
    para.sen <- Update(para.sen, sen, sen.range[i])
    
    # Run the initiation and iterations
    if(sen %in% c("rm", "rsd", "size", "gen", "fit", "fq")){
      equ.sen <- equaldensity
    }else{
      # Update the initial field population and then let it reach equilibrium
      # Firest updatae migrating population
      if(para.sen$mig){
        para.sen.ext <- para.sen
        para.sen.ext$mig <- 0
        ext.density0.sen <- Create(pop = "field", para.sen.ext)
        ext.sen <- Initial(ext.density0.sen, para.sen.ext, alter)
        ext.density.sen <- ext.sen$density
        para.sen$migrant <- ext.density.sen/Popsize(simpson = para.sen.ext$Sz, density = ext.density.sen)
      }
      # Then update the pre-releasing target population
      density0.sen <- Create(pop = "field", para.sen)
      equ.sen.list <- Initial(density0.sen, para.sen, alter)
      equ.sen <- equ.sen.list$density
    }
    
    # Update the releasing population if necessary
    para.sen$Nr <- para.sen$size * Popsize(para.sen$Sz, density = equ.sen)
    para.sen$rel <- Create(pop = "releasing", para.sen)
    
    # Set a criteria to calculate the rate of changes of mean genotype: reduce two sd of the pre-releasing population
    para.sen$cre <- Genotype(equ.sen, para.sen) - 2*sqrt(GVariance(equ.sen, para.sen))

    # Run the releasing process
    out.sen <- Iterate(equ.sen, para.sen, alter)
    
    # Relapse
    relapse.sen <- Relapse(out.sen$density, para.sen, alter, g = relapse_gen)
    
    return(list(o=out.sen, r=relapse.sen, g=para.sen$gen, 
                speed=which(out.sen$postmean < para.sen$cre)[1]))
  }
  
  # Record the outcome of the post-releasing population
  for(i in 1:l){
    out.sen.result <- parallel_results[[i]]$o
    relapse.sen.result <- parallel_results[[i]]$r
    
    # Combine the releasing and relapsing results:
    sen_pop[,i] <- c(out.sen.result$postpop, relapse.sen.result$inipop[2:(relapse_gen+1)])
    sen_mean[,i] <- c(out.sen.result$postmean, relapse.sen.result$inimean[2:(relapse_gen+1)])
    sen_var[,i] <- c(out.sen.result$postvar, relapse.sen.result$inivar[2:(relapse_gen+1)])
    sen_mean_vc[,i] <- c(out.sen.result$postmean_vc, relapse.sen.result$inimean_vc[2:(relapse_gen+1)])
    sen_var_vc[,i] <- c(out.sen.result$postvar_vc, relapse.sen.result$inivar_vc[2:(relapse_gen+1)])
    sen_fit[,i] <- c(out.sen.result$postfit, relapse.sen.result$inifit[2:(relapse_gen+1)])
    sen_integrateVC[,i] <- c(out.sen.result$postsum, relapse.sen.result$inisum[2:(relapse_gen+1)])
    sen_speed[i] <- parallel_results[[i]]$speed
  }
  
  return(list(sen.range=sen.range, sen.gen=para$gen,
              sen_pop=sen_pop, sen_mean=sen_mean, sen_var=sen_var, 
              sen_mean_vc=sen_mean_vc, sen_var_vc=sen_var_vc, 
              sen_fit=sen_fit, sen_integrateVC = sen_integrateVC,
              sen_speed=sen_speed))
}


# Function to extract key information from local sensitivity analysis
Sensummary <- function(out, toi=c(-1,-1)){
  # Input:  - out: the output from function "Locals".
  #         - sen: the parameter of interest. e.g. "rm", "R" etc.
  #         - equsum: the equilibrium population size before releasing BY DEFAULT SETTING
  #         - toi: indicator of which generation to look at (default = looking at the last generation)
  
  # Output: a pdf file containing four figures
  
  # Results after releasing for "sen.gen" generations
  if (toi[1]==-1){
    o <- out$sen.gen+1  # default: the "out$sen.gen+1" row in the result matrices
  }else{
    o <- toi[1]         # Look at specific generations ("toi") in the process
  }
  
  # Results after relapsing
  if (toi[2]==-1){
    r <- nrow(out$sen_pop)  # default: the last row in the result matrices
  }else{
    r <- toi[2]         # Look at specific generations ("toi") in the process
  }
  
  sen.summary <- data.frame(out.range = out$sen.range,
                            equ.pop = out$sen_pop[1,],
                            equ.mean = out$sen_mean[1,],
                            equ.var = out$sen_var[1,],
                            equ.mean.vc = out$sen_mean_vc[1,],
                            equ.var.vc = out$sen_var_vc[1,],
                            equ.fit = out$sen_fit[1,],
                            equ.integrateVC = out$sen_integrateVC[1,],
                            out.pop = out$sen_pop[o,],
                            out.mean = out$sen_mean[o,],
                            out.var = out$sen_var[o,],
                            out.mean.vc = out$sen_mean_vc[o,],
                            out.var.vc = out$sen_var_vc[o,],
                            out.fit = out$sen_fit[o,],
                            out.integrateVC = out$sen_integrateVC[o,],
                            speed = out$sen_speed,
                            relapse.pop = out$sen_pop[r,],
                            relapse.mean = out$sen_mean[r,],
                            relapse.var = out$sen_var[r,],
                            relapse.mean.vc = out$sen_mean_vc[r,],
                            relapse.var.vc = out$sen_var_vc[r,],
                            relapse.fit = out$sen_fit[r,],
                            relapse.integrateVC = out$sen_integrateVC[r,])
  return(sen.summary)
}


# Translate the parameter name (e.g. "rm") into descriptions (e.g. "VC mean of the releasing population"):
Annotation <- function(sen){
  if (sen=="rm"){
    return ("Mean genotype of the releasing population")
  }else if(sen=="rm.prop"){
    return ("Ratio of mean genotype of releasing population against field population")
  }else if(sen=="rsd"){
    return ("Standard deviation(SD) of genotype of the releasing population")
  }else if(sen=="rsd.prop"){
    return ("Ratio of genotype sd of the releasing population against field population")
  }else if(sen=="size"){
    return ("Relative size of the releasing population")
  }else if(sen=="fit"){
    return ("Relative fitness of the releasing population")
  }else if(sen=="fm"){
    return ("Optimal genotype in the field")
  }else if(sen=="fsd"){
    return ("SD of the initial field population")
  }else if(sen=="N0"){
    return ("Size of the initial field population")
  }else if(sen=="R"){
    return ("Reproductive rate per female")
  }else if(sen=="vp"){
    return ("Total phenotype variance")
  }else if(sen=="h"){
    return ("Heritability of VC")
  }else if(sen=="alpha"){
    return ("Beverton-Holt density-dependent parameter")
  }else if(sen=="vs"){
    return ("Variance in selection surface")
  }else if(sen=="ind2"){
    return ("Density independent survival in pupae emerging")
  }else if(sen=="Nm"){
    return ("Absolute number of immigrants per generation")
  }else if(sen=="post_mean"){
    return ("Mean genotype after releasing")
  }else if(sen=="prop"){
    return ("Proportion of integrated VC")
  }
}


#----# Selectively summarize important results
# MOST IMPORTANT VARIABLES:
# a. mean VC (phenotype) after releasing
# b. Proportional remaing of integrated VC (phenotype) after releasing relative to the
#    pre-releasing equilibrium population.

Report <- function(out){
  mean_VC <- out$postmean_vc[length(out$postmean_vc)]
  integrate_VC_pre <- out$postmean_vc[1] * out$postpop[1]
  integrate_VC_post <- mean_VC * out$postpop[length(out$postpop)]
  prop <- integrate_VC_post / integrate_VC_pre
  return(m=mean_VC, prop=prop)
}





#-------------------------------- Model alternatives --------------------------------

# 1) Which reproduction function is used
alter <- 2  
# alter=1: Corresponding to Reproduction 1, use FFT to calculate convolution of the parental distribution.
#          Can only be used when both sexes are released.          
#          (Numerical error is too large)
# alter=2: Corresponding to Reproduction 2, use basic matrix calculation.
#          The two sexes can have different genotype distributions.
#          Can be used when either sex or both sexes are released.
#          (default methods, less numerical error but much slower)


# 4) Order of the life cycle events
selection_time <- "before_releasing"
# selection_time="before_releasing": selection happens only before migration and/or releasing
# selection_time="after_releasing": selection happens only after migration and/or releasing
# selection_time="both": selection happens both before and after migration and/or releasing

# 5) Releasing sex:
releasing_sex <- "male"
# releasing_sex="both": releasing both males and females
# releasing_sex="male": releasing only males (Only alter=2 is appliable)
# releasing_sex="bloodfed": releasing males and blood-fed females

#-----# Select the function "Initial", "Iterate" and "Relapse" based on the order of life cycle events
# and sex of releasing individuals
if(selection_time=="before_releasing"){
  Initial <- Initial_before
  Relapse <- Relapse_before
  if(releasing_sex=="both"){
    Iterate <- Iterate_before_both
  }else if(releasing_sex=="male"){
    Iterate <- Iterate_before_male
  }else if(releasing_sex=="bloodfed"){
    Iterate <- Iterate_before_bloodfed
  }else{
    Iterate <- Iterate_before_both
    print("No valid input of releasing_sex, default method used: releasing both sexes.", immediate. = TRUE)
  }
}else if(selection_time=="after_releasing"){
  Initial <- Initial_after
  Relapse <- Relapse_after
  if(releasing_sex=="both"){
    Iterate <- Iterate_after_both
  }else if(releasing_sex=="male"){
    Iterate <- Iterate_after_male
  }else if(releasing_sex=="bloodfed"){
    Iterate <- Iterate_after_bloodfed
  }else{
    Iterate <- Iterate_after_both
    print("No valid input of releasing_sex, default method used: releasing both sexes.", immediate. = TRUE)
  }
}else if(selection_time=="both"){
  Initial <- Initial_both
  Relapse <- Relapse_both
  if(releasing_sex=="both"){
    Iterate <- Iterate_both_both
  }else if(releasing_sex=="male"){
    Iterate <- Iterate_both_male
  }else if(releasing_sex=="bloodfed"){
    Iterate <- Iterate_both_bloodfed
  }else{
    Iterate <- Iterate_both_both
    print("No valid input of releasing_sex, default method used: releasing both sexes.", immediate. = TRUE)
  }
}else{
  Initial <- Initial_before
  Relapse <- Relapse_before
  print("No valid input of selection_time, default method used: selection before releasing.", immediate. = TRUE)
  
  if(releasing_sex=="both"){
    Iterate <- Iterate_before_both
  }else if(releasing_sex=="male"){
    Iterate <- Iterate_before_male
  }else if(releasing_sex=="bloodfed"){
    Iterate <- Iterate_before_bloodfed
  }else{
    Iterate <- Iterate_before_both
    print("No valid input of releasing_sex, default method used: releasing both sexes.", immediate. = TRUE)
  }
}


# Locals <- cmpfun(Locals)







#---------------------------------- Parameters---------------------------------------

par.input <- read.csv("GSM_parameters_manual.csv", 
                      header = TRUE, stringsAsFactors = FALSE)

#----# Register cores and clusters
ncore <- 4          # Number of CPUs to use!

if(getDoParWorkers()==1){
  cl <- makeCluster(ncore)
  registerDoParallel(cl)
}
getDoParWorkers()

# Record function number in a list: required input for parallel computing
function_list <- c("Annotation", "Create", "Den.dependent", "Den.independent", "Densummary", 
                   "Fitness", "Genotype", "GtoPh", "GVariance", "Initial", "Iterate", "Relapse",
                   "Maintain", "Migration", "Phenotype", "PhVariance", 
                   "Popsize", "Releasing1", "Releasing2", "Report", "reproduction1", "reproduction2", 
                   "Selection", "Simpson", "Update", "step.change")


#--- Model implementation
d <- foreach(i = 1:nrow(par.input), .packages="truncnorm", .export = function_list) %dopar% {
  
  #----# Modeling paprameters
  # Observed optimal VC in the field
  fm <- par.input$fm[i]     # initial (optimal) mean vector competence of field population
  N0 <- 100000   # initial size of field population
  
  # LIfe history parameters
  R <- 40            # mean reproductive rate per female (taking ind1 into consideration)
  vp <- par.input$vp[i]         # Total phenotypic variance
  h <- 0.4           # heritability
  vle <- vp*h        # genetic variance at LE
  ve <- (1-h)*vp     # environmental variance (genotype->phenotype)
  fsd <- sqrt(vle)   # initial sd of vector competence of the field population: determined by genetic variance at LE (vle)
  alpha <- 0.0001   # Beverton-Holt density-dependent parameter
  vs <- 1            # Variance in selection surface (1/selection strength)
  ind2 <- 0.7        # Density independent survival in pupae emerging
  sex <- 1           # Sex ratio: male:female = sex:1
  # ind1 <- 0.15  # Density independent survival in egg hatching (now included into R)
  
  # Immigration
  Nm <- 0      # absolute number of immigrating individuals per generation (default: no migration)
  
  # Releasing population
  rm <- par.input$rm[i]        # mean vector competence of lab strain (releasing population)
  rsd <- par.input$rsd[i]      # sd of vector competence of the lab strain (releasing population)
  size <- par.input$size[i]      # size of releasing population relative to the field population
  fit <- 0.75      # relative fitness of the releasing population (relative to field population)
  R.rel <- 50      # mean reproductive rate per releasing female after blood-fed before releasing
                   #  (slightly larger than females in the field as they are less likely to get a full blood meal)
  
  # Length of the releaing program
  gen <- 10    # how many generations of releasing: 20 generations approximately equals one year
  fq <- 1      # Releasing frequency: how many generations per release. fq = {1, 2, 3, 4, 5}
  
  
  #----# auxiliary variables (from Prof. Marissa Baskett)
  # basic matrices, vectors, and values for integration etc.
  zmin = 0                # minimum genotype -- should be much less than expected minimum
  zmax = 1                # maximum genotype -- should be greater than expected max
  Nz = 2^floor(log((zmax-zmin)/0.001,base=2))     # number of genotype entries -- for fft algorithm efficiency
  vz = seq(zmin,zmax,length=Nz)                   # vector of genotypes
  dz = diff(vz)[1]                                # step for genotype
  Sz = Simpson(dz,Nz)                             # Simpson vector for integrating over genotypes
  
  # For FFT implements:
  midl = seq(1,Nz-1,by=2) # middle values to extract
  Trans = dnorm(vz,mean=vz[Nz/2+1], sd=sqrt(vle/2)) # transmission function
  rotate = c((Nz/2+1):Nz,1:(Nz/2))
  Trans[rotate] = Trans
  
  # Create a vector to store the density distribution of the immigrant (distribution scaled by pop size)
  migrant <- rep(0, length(vz))      # Density distribution of migrants
  
  
  #----# Organizing all parameters into a list:
  para <- list(fm=fm, fsd=fsd,                                                      # Field population (target population)
               rm=rm, rsd=rsd, size=size, R.rel=R.rel,                              # Releasing population
               vp=vp, h=h, vle=vle, ve=ve, N0=N0,                                   # Variance composition
               R=R, fit=fit, sex=sex, alpha=alpha, vs=vs, ind2=ind2,                # Life cycle parameters
               gen=gen, fq=fq,                                                      # Time steps
               Nm=Nm, migrant=migrant,                                              # Immigrants
               zmin=zmin,zmax=zmax,Nz=Nz,vz=vz,dz=dz,Sz=Sz,                         # Numerically implements
               midl=midl, Trans=Trans)                                              # FFT
  
  # Whether to calculate phenotype (vector competence):
  para$phenotype <- FALSE
  
  # Whether there are immigrants in the population
  para$mig <- 0
    # para$mig=1: Nm Immigrants are moving into the population every generation
    # para$mig=0: No immigrants
  
  #-------------------------------- Check all model options -------------------------------
  print(paste0("Parameter combination: ", par.input$Index[i]))
  print(paste0("Releasing strategy: ", releasing_sex))
  # print(para$mig)
  # print(para$phenotype)
  # print(selection_time)
  # print(releasing_sex)
  
  
  
  #--------------------------------- Model implements ----------------------------
  
  
  #----# initialize the external field population if considering migration
  if(para$mig){
    # parameters for the external population
    para_external <- para
    para_external$mig <- 0
    
    # Initial field population:
    density_external <- Create(pop = "field", para_external)
    
    # Reaching equilibrium status before releasing
    ext <- Initial(density_external, para_external, alter)
    
    # Record the equilibrium density distribution
    ext.density <- ext$density
    
    # Summary of the equilibrium population status before releasing.
    ext.summary <- Densummary(status = ext, para = para)
    # Checked analytically based on normal approximation 
    # (See the pdf file "Iterations" provided by Prof. Baskett)   
    
    # modify the parameter set to include immigrants:
    para$migrant <- ext.density/Popsize(simpson = para$Sz, density = ext.density)
  }
  
  
  
  #----# initialize the TARGET field population
  # Initial field population:
  density0 <- Create(pop = "field", para)
  
  # Reaching equilibrium status before releasing
  equ <- Initial(density0, para, alter)
  
  # Record the equilibrium density distribution
  equ.density <- equ$density
  
  # Summary of the equilibrium population status before releasing.
  equ.summary <- Densummary(status = equ, para = para)
  
  # Put all metrics into a data frame for output
  equ.subset <- 1:length(equ$init)
  equ.bind <- cbind(equ$init, as.data.frame(equ[3:9])[equ.subset,])
  names(equ.bind)[1] <- "init"
  
  # wirte out the data:
  # write.csv(equ.bind, paste0("Pre-releasing_Index_", par.input$Index[i], ".csv"), 
  #           row.names = FALSE)
  
       
  
  #----# Set up the releasing population:
  # Distribution of the Lab strain / Releasing population:
  Nr <- para$size*equ.summary$dpop          # Size of the releasing pop
  para$Nr <- Nr
  
  rel <- Create(pop = "releasing", para)    # Distribution of the releasing pop
  para$rel <- rel
  
  
  #----# Iterate releasing procuss:
  # Outcome of iterate the releasing program for XXX generations using the default settings
  outcome <- Iterate(equ.density, para, alter)
  outcome.summary <- Densummary(outcome, para)
  outcome.density <- outcome$density            # Density distribution after the last releasing
  
  # Put all metrics into a data frame for output
  outcome.bind <- as.data.frame(outcome[c(1,3:9)])
  
  # wirte out the data:
  write.csv(outcome.bind, paste0("Releasing_Index_", par.input$Index[i], ".csv"),
            row.names = FALSE)
  
  return(outcome.bind)
}

save.image("Genetic_shifting_manual_male.RData")


# load("Genetic_shifting_bloodfed.RData")

for(j in 1:nrow(par.input)){
  par.input[j, paste0("VC_mean.", 1:10)] <- d[[j]]$postmean[1:10]
  par.input[j, paste0("VC_sd.", 1:10)] <- sqrt(d[[j]]$postvar[1:10] / 0.4)
}

write.csv(x = par.input, file = paste0("GSM_manual_", releasing_sex, ".csv"), 
          quote = FALSE, row.names = FALSE)

# #----# After releasing stops:
# # Changes for 20 generations after the releasing stops
# relapse <- Relapse(outcome.density, para, alter)     
# relapse.summary <- Densummary(relapse, para)         # Summarizing the results
# 
# # Put all metrics into a data frame for output
# relapse.bind <- as.data.frame(relapse[c(1,3:9)])
# 
# # wirte out the data:
# write.csv(relapse.bind, "Relapse_default.csv", row.names = FALSE)
# 
# 
# #----# Output some key density distribution
# savepoint <- as.data.frame(cbind(equ.summary$density, outcome.summary$density, relapse.summary$density))
# names(savepoint) <- c("Pre-releasing", "Post-releasing", "Relapse")
# write.csv(savepoint, "Genotype distribution pre-releasing, post-releasing and post-relapsing.csv")



#################################################################################################################
#################################################################################################################

# #--------------------------------- Sensitivity analysis preparation ----------------------------------
# 
# #----# 0. Parameter ranges
# rm.prop.range <- seq(from = 0, to = 1, by = 0.01)      # rm.range should be proportional to the equilibrium population mean
# rm.range <- rm.prop.range * fm
# rsd.prop.range <- seq(from = 0.01, to = 1, by = 0.01)  # rsd.range should be proportional to the equilibrium population sd
# rsd.range <- rsd.prop.range * fsd
# gen.range <- seq(from = 1, to = 50, by = 1)
# fit.range <- seq(from = 0.01, to = 1, by = 0.01)
# size.range <- seq(from = 0.01, to = 0.5, by = 0.01)
# fq.range <- 1:5                                        # Releasing every 1, 2, 3, 4 or 5 generations
# gen.range <- 1:50                                      # Releasing for 1-50 generations
# R.rel.range <- seq(from = 20, to = 150, by = 5)        # Reproductive rate per females of the releasing population
#                                                        # (only apply for the scenario where blood-fed females were released)
# 
# fm.range <- seq(from = 0.1, to = 0.9, by = 0.01)
# R.range <- seq(from = 5, to = 150, by = 5)
# vp.range <- seq(from = 0.01, to = 0.5, by = 0.005)^2                   # square scale
# h.range <- seq(from = 0.01, to = 1, by = 0.01)
# alpha.range <- seq(0.00001, 0.001, length.out = 100)
# vs.range <- seq(0.1, 10, length.out = 100)^2                           # square scale
# ind2.range <- seq(from = 0.2, to = 1, by = 0.01)
# Nm.range <- seq(from = 0, to = 500, by = 10)                          # Population size of the immigrations
#
#
# #----# Register cores and clusters
# ncore <- 20          # Number of CPUs to use!
# 
# if(getDoParWorkers()==1){
#   cl <- makeCluster(ncore)
#   registerDoParallel(cl)
# }
# getDoParWorkers()
#
# #--------------------------------- Local sensitivity analysis ----------------------------------------
# 
# #----# 1. Parameters of the releasing population
# # In this part, the equilibrium field population before releasing do not change, only change parameters
# # in the function "Iterate()".
# sen.alter <- 2
# 
# # 1.1 mean vector competence of the releasing population: rm
# rm.systime <- system.time(rm.out <- Locals(sen="rm", sen.range = rm.range, para=para, alter=sen.alter))
# rm.summary <- Sensummary(rm.out)      # Each row represent the final results of 20 generations of releasing given a variable value
# write.csv(rm.summary, "LSA_rm.csv", row.names = FALSE)
# 
# # 1.2 Standard deviation of vector competence of the releasing population: rsd
# rsd.systime <- system.time(rsd.out <- Locals(sen="rsd", sen.range = rsd.range, para=para, alter=sen.alter))
# rsd.summary <- Sensummary(rsd.out)      # Each row represent the final results of 20 generations of releasing given a variable value
# write.csv(rsd.summary, "LSA_rsd.csv", row.names = FALSE)
# 
# # 1.3 Size of the releasing population: Nr
# size.systime <- system.time(size.out <- Locals(sen="size", sen.range = size.range, para=para, alter=sen.alter))
# size.summary <- Sensummary(size.out)      # Each row represent the final results of 20 generations of releasing given a variable value
# write.csv(size.summary, "LSA_size.csv", row.names = FALSE)
# 
# # 1.4 The relative fitness of releasing population: fit
# fit.systime <- system.time(fit.out <- Locals(sen="fit", sen.range = fit.range, para=para, alter=sen.alter))
# fit.summary <- Sensummary(fit.out)      # Each row represent the final results of 20 generations of releasing given a variable value
# write.csv(fit.summary, "LSA_fit.csv", row.names = FALSE)
# 
# # 1.5 releasing frequency: fq
# fq.systime <- system.time(fq.out <- Locals(sen="fq", sen.range = fq.range, para=para, alter=sen.alter))
# fq.summary <- Sensummary(fq.out)      # Each row represent the final results of 20 generations of releasing given a variable value
# write.csv(fq.summary, "LSA_fq.csv", row.names = FALSE)
# 
# 
# ###### 2. Life history parameters
# # The change of life history parameters will influence both the initial process as well as the releasing
# # 2.1 The reproductive rate: R
# R.systime <- system.time(R.out <- Locals(sen="R", sen.range = R.range, para=para, alter=sen.alter))
# R.summary <- Sensummary(R.out)      # Each row represent the final results of 20 generations of releasing given a variable value
# write.csv(R.summary, "LSA_R.csv", row.names = FALSE)
# 
# # 2.2 The total phenotypic variance: vp
# vp.systime <- system.time(vp.out <- Locals(sen="vp", sen.range = vp.range, para=para, alter=sen.alter))
# vp.summary <- Sensummary(vp.out)      # Each row represent the final results of 20 generations of releasing given a variable value
# write.csv(vp.summary, "LSA_vp.csv", row.names = FALSE)
# 
# # 2.3 Heritability: h
# h.systime <- system.time(h.out <- Locals(sen="h", sen.range = h.range, para=para, alter=sen.alter))
# h.summary <- Sensummary(h.out)      # Each row represent the final results of 20 generations of releasing given a variable value
# write.csv(h.summary, "LSA_h.csv", row.names = FALSE)
# 
# # 2.4 Beverton-Holt density-dependent parameter: alpha
# alpha.systime <- system.time(alpha.out <- Locals(sen="alpha", sen.range = alpha.range, para=para, alter=sen.alter))
# alpha.summary <- Sensummary(alpha.out)      # Each row represent the final results of 20 generations of releasing given a variable value
# write.csv(alpha.summary, "LSA_alpha.csv", row.names = FALSE)
# 
# # 2.5 Variance of selection surface: vs (1/selection strength)
# vs.systime <- system.time(vs.out <- Locals(sen="vs", sen.range = vs.range, para=para, alter=sen.alter))
# vs.summary <- Sensummary(vs.out)      # Each row represent the final results of 20 generations of releasing given a variable value
# write.csv(vs.summary, "LSA_vs.csv", row.names = FALSE)
# 
# # 2.6 Density independence of pupa emerging: ind2
# ind2.systime <- system.time(ind2.out <- Locals(sen="ind2", sen.range = ind2.range, para=para, alter=sen.alter))
# ind2.summary <- Sensummary(ind2.out)      # Each row represent the final results of 20 generations of releasing given a variable value
# write.csv(ind2.summary, "LSA_ind2.csv", row.names = FALSE)
# 
# # 2.7 External immigrants: Nm
# Nm.systime <- system.time(Nm.out <- Locals(sen="Nm", sen.range = Nm.range, para=para, alter=sen.alter))
# Nm.summary <- Sensummary(Nm.out)      # Each row represent the final results of 20 generations of releasing given a variable value
# write.csv(Nm.summary, "LSA_Nm.csv", row.names = FALSE)
# 
# # 3. Optimal genotype in the field: fm
# fm.systime <- system.time(fm.out <- Locals(sen="fm", sen.range = fm.range, para=para, alter=sen.alter))
# fm.summary <- Sensummary(fm.out)      # Each row represent the final results of 20 generations of releasing given a variable value
# write.csv(fm.summary, "LSA_fm.csv", row.names = FALSE)
# 
# # 4. Number of generations of releasing: gen (no relapse)
# para.gen <- para
# para.gen$gen <- 50
# out.gen <- Iterate(equaldensity = equ.density, para = para.gen, alter = sen.alter)
# # rel.gen <- Relapse(post_releasing_density = out.gen$density, para = para.gen, alter = sen.alter)
# 
# gen.summary <- data.frame(out.range = out.gen$time,
#                           equ.pop = out.gen$postpop[1],
#                           equ.mean = out.gen$postmean[1],
#                           equ.var = out.gen$postvar[1],
#                           equ.mean.vc = out.gen$postmean_vc[1],
#                           equ.var.vc = out.gen$postvar_vc[1],
#                           equ.fit = out.gen$postfit[1],
#                           equ.integrateVC = out.gen$postsum[1],
#                           out.pop = out.gen$postpop,
#                           out.mean = out.gen$postmean,
#                           out.var = out.gen$postvar,
#                           out.mean.vc = out.gen$postmean_vc,
#                           out.var.vc = out.gen$postvar_vc,
#                           out.fit = out.gen$postfit,
#                           out.integrateVC = out.gen$postsum)
# write.csv(gen.summary, "LSA_gen.csv", row.names = FALSE)
# 
# # 5. Reproduction rate per female of the releasing population
# if(releasing_sex == "bloodfed"){
#   R.rel.systime <- system.time(R.rel.out <- Locals(sen="R.rel", sen.range = R.rel.range, para=para, alter=sen.alter))
#   R.rel.summary <- Sensummary(R.rel.out)      # Each row represent the final results of 20 generations of releasing given a variable value
#   write.csv(R.rel.summary, "LSA_R.rel.csv", row.names = FALSE)
# }
# 
# 
# #----# Report system time for LSA
# LSA.systime <- as.data.frame(rbind(rm.systime, rsd.systime, size.systime, fit.systime,
#                                    R.systime, vp.systime, h.systime, alpha.systime, vs.systime, ind2.systime, Nm.systime))
# write.csv(LSA.systime, "systime of LSA.csv", row.names = TRUE)






########################################################################################
########################################################################################

# # -------------------------- Global Sensitivity Analysis -----------------------
# #-----# Preparation
# # MC sampling times
# N_sample <- 5000
# 
# #----# 0. Parameter ranges: in finer resolution than LSA
# rm.prop.candidate <- seq(from = 0, to = 1, by = 0.001)      # rm.candidate should be proportional to the equilibrium population mean
# rsd.prop.candidate <- seq(from = 0.01, to = 1, by = 0.001)  # rsd.candidate should be proportional to the equilibrium population sd
# gen.candidate <- seq(from = 1, to = 50, by = 1)
# fit.candidate <- seq(from = 0.01, to = 1, by = 0.001)
# size.candidate <- seq(from = 0.01, to = 0.5, by = 0.001)
# fq.candidate <- 1:5                                        # Releasing every 1, 2, 3, 4 or 5 generations
# R.rel.candidate <- seq(from = 20, to = 150, by = 1)              # Reproductive rate per females of the releasing population
# # (only apply for the scenario where blood-fed females were released)
# 
# 
# fm.candidate <- seq(from = 0.1, to = 0.9, by = 0.001)
# R.candidate <- seq(from = 5, to = 150, by = 1)
# vp.candidate <- seq(from = 0.01, to = 0.5, by = 0.001)^2                   # square scale
# h.candidate <- seq(from = 0.01, to = 1, by = 0.001)
# alpha.candidate <- seq(0.00001, 0.001, length.out = 500)
# vs.candidate <- seq(0.1, 10, by = 0.01)^2                           # square scale
# ind2.candidate <- seq(from = 0.2, to = 1, by = 0.001)
# Nm.candidate <- seq(from = 0, to = 500, by = 10)                          # Population size of the immigrations
# 
# # Use a list to record ranges of all variables
# candidate_list <- list(rm = rm.prop.candidate,
#                        rsd = rsd.prop.candidate,
#                        size = size.candidate,
#                        fit = fit.candidate,
#                        gen = gen.candidate,
#                        fq = fq.candidate,
#                        fm = fm.candidate,
#                        R = R.candidate,
#                        vp = vp.candidate,
#                        h = h.candidate,
#                        alpha = alpha.candidate,
#                        vs = vs.candidate,
#                        ind2 = ind2.candidate,
#                        Nm = Nm.candidate,
#                        R.rel = R.rel.candidate)
# 
# 
# #-----# Which variables are allowed to vary:
# # a. all variables
# list_para <- c("fm", "rm", "rsd", "size", "fit", "fq", "R.rel",
#                "R", "vp", "h", "alpha", "vs", "ind2", "gen", "Nm")
# # b. only manageable variables
# list_para_manageable <- c("rm", "rsd", "size", "fit", "gen", "fq", "R.rel")
# 
# # Default parameter settings:
# para.global <- para
# 
# 
# ########################################################################################################
# #-----# Function to implement global sensitivity analysis
# Global <- function(para.global, alter, list_para, N_sample = 5000, candidate=candidate_list, fix=FALSE, var_mat=NULL){
# 
#   #---# Calculate pre-defined initial population (used when only vary manageable variables):
#   if(all(list_para %in% c("rm", "rsd", "size", "fit", "gen", "fq"))){
#     # external population:
#     if(para.global$mig==1){
#       para.global.ext <- para.global
#       para.global.ext$mig <- 0
#       density.ext.global <- Create(pop = "field", para.global.ext)
#       ext.global <- Initial(density.ext.global, para.global.ext, alter)
#       para.global$migrant <- ext.global$density / Popsize(para.global$Sz, ext.global$density)
#     }else{
#       para.global$migrant <- rep(0, length(para.global$vz))
#     }
# 
#     # Pre-releasing population:
#     density0.global.predefine <- Create(pop = "field", para.global)
#     equ.predefine <- Initial(density0.global.predefine, para.global, alter)
#   }
# 
# 
#   #-----# Create a data frame to record all metrics of interest:
#   list_results <- c("pre_mean", "pre_var", "pre_pop", "pre_fitness",
#                     "post_mean", "post_var", "post_pop", "post_fitness", "speed")
#   list_phenotype <- c("pre_mean_phenotype", "pre_var_phenotype", "post_mean_phenotype", "post_var_phenotype")
# 
#   # Combine all lists together
#   list_all <- c(list_para, list_results)
# 
#   # Whether to consider phenotype (explicit VC)
#   if(para.global$phenotype){
#     list_all <- c(list_all, list_phenotype)      # include 4 more variables to record phenotype values
#   }
# 
#   n_sub <- 1000                   # Number of runs before each recording points
#   n_run <- floor(N_sample/n_sub)  # Number of recording points
# 
#   # Create an empty data frame to record all results from each checking points
#   global.results.sum <- NULL
# 
#   for(rec in 1:n_run){
# 
#     #-----# MC sampling and implement of releasing model
#     global.results <- foreach(i = 1:n_sub, .packages="truncnorm", .export = function_list, .combine = rbind, .multicombine = TRUE) %dopar% {
#       #-# create a vector to store results from each iteration
#       sample.results <- rep(NA, length(list_all))
#       names(sample.results) <- list_all
# 
#       #-# Sample parameter values and update parameter list:
#       if(fix){
#         for(p in list_para) para.global[[p]] <- var_mat[[p]][i]
#       }else{
#         for(p in list_para) para.global[[p]] <- sample(candidate[[p]], size = 1)
#       }
# 
#       # Some calculation to update para.global
#       para.global$vle <- para.global$vp * para.global$h
#       para.global$ve <- para.global$vp * (1 - para.global$h)
#       para.global$fsd <- sqrt(para.global$vle)
# 
#       para.global$Trans <- dnorm(para.global$vz,
#                                  mean=para.global$vz[para.global$Nz/2+1],
#                                  sd=sqrt(para.global$vle/2))
#       para.global$rotate = c((para.global$Nz/2+1):para.global$Nz,1:(para.global$Nz/2))
#       para.global$Trans[para.global$rotate] = para.global$Trans
# 
# 
#       #-# Pre-releasing population:
#       if(all(list_para %in% c("rm", "rsd", "size", "fit", "gen"))){
#         equ.global.list <- equ.predefine
#       }else{
#         # update the external population:
#         if(para.global$mig==1){
#           para.global.ext <- para.global
#           para.global.ext$mig <- 0
#           density.ext.global <- Create(pop = "field", para.global.ext)
#           ext.global <- Initial(density.ext.global, para.global.ext, alter)
#           para.global$migrant <- ext.global$density / Popsize(para.global$Sz, ext.global$density)
#         }else{
#           para.global$migrant <- rep(0, length(para.global$vz))
#         }
# 
#         # target population before releasing:
#         density0.global <- Create(pop = "field", para.global)
#         equ.global.list <- Initial(density0.global, para.global, alter)
#       }
# 
#       #-# Summary of the pre-releasing equilibrium population
#       equ.global.density <- equ.global.list$density
#       equ.global.summary <- Densummary(status = equ.global.list, para = para.global)
# 
#       #-# Update information about releasing population in the parameter list
#       para.global$rm <- para.global$rm * equ.global.summary$dmean         # mean genotype of the releasing pop (proportion to the pre-releasing pop)
#       para.global$rsd <- para.global$rsd * sqrt(equ.global.summary$dvar)  # genotye sd of the releasing pop (proportion to the pre-releaing pop)
#       para.global$Nr <- para.global$size * equ.global.summary$dpop        # Pop size of the releasing pop (proportion to the pre-releasing pop)
#       para.global$rel  <- Create(pop = "releasing", para.global)
# 
#       # Set a criteria to calculate the rate of changes of mean genotype: reduce two sd of the pre-releasing population
#       para.global$cre <- equ.global.summary$dmean  - 2*sqrt(equ.global.summary$dvar)
# 
#       #---# Post-releasing population
#       outcome.global <- Iterate(equ.global.density, para.global, alter)
#       outcome.global.summary <- Densummary(outcome.global, para.global)
#       #---#
# 
#       # Record all parameters and metrics of interest:
#       for(p in list_para) sample.results[p] <- para.global[[p]]
# 
#       sample.results["pre_mean"] <- equ.global.summary$dmean
#       sample.results["pre_var"] <- equ.global.summary$dvar
#       sample.results["pre_pop"] <- equ.global.summary$dpop
#       sample.results["pre_fitness"] <- equ.global.summary$dfit
# 
#       sample.results["post_mean"] <- outcome.global.summary$dmean
#       sample.results["post_var"] <- outcome.global.summary$dvar
#       sample.results["post_pop"] <- outcome.global.summary$dpop
#       sample.results["post_fitness"] <- outcome.global.summary$dfit
#       sample.results["speed"] <- which(outcome.global$postmean < para.global$cre)[1]
# 
#       if(para.global$phenotype){
#         sample.results["pre_mean_phenotype"] <- equ.global.summary$vcmean
#         sample.results["pre_var_phenotype"] <- equ.global.summary$vcvar
# 
#         sample.results["post_mean_phenotype"] <- outcome.global.summary$vcmean
#         sample.results["post_var_phenotype"] <- outcome.global.summary$vcvar
#       }
#       # print(paste("MC sample:", i))
#       return(sample.results)
#     }
# 
#     # Make a data frame from all parallel computing results
#     global.results <- as.data.frame(global.results)
# 
#     # Calculating remaining proportion of integrated VC:
#     global.results$prop <- (global.results$post_mean * global.results$post_pop) / (global.results$pre_mean * global.results$pre_pop)
# 
#     # Write out the results:
#     write.csv(global.results, paste0("matrix of global sensitivity analysis_all variables_N=", N_sample, "_subset=", rec, ".csv"), row.names=FALSE)
# 
#     # Storing all results in a single data frame
#     global.results.sum <- rbind(global.results.sum, global.results)
#   }
# 
#   return(list(global.results = global.results.sum, list_para = list_para, N_sample = N_sample))
# }
# 
# 
# #----# Register cores and clusters
# if(getDoParWorkers()==1){
#   cl <- makeCluster(ncore)
#   registerDoParallel(cl)
# }
# getDoParWorkers()
# 
# 
# #-----# a. all variables
# # Read in the variable matrix
# csv_finder <- list.files(pattern = "variable_matrix*")
# if(length(csv_finder)==1){
#   variable_matrix <- read.csv(csv_finder, header = TRUE, stringsAsFactors = FALSE)
#   print(csv_finder)
# }else{
#   print("More than one variable matrices exist!!")
# }
# 
# # Run GSA
# global.results.default.all <- Global(para.global = para, list_para = list_para, N_sample = N_sample, alter = 2,
#                                      fix = TRUE, var_mat = variable_matrix)
# global.results <- global.results.default.all$global.results
# write.csv(global.results, paste0("matrix of global sensitivity analysis_all variables_N=", 
#                                  global.results.default.all$N_sample, ".csv"), row.names=FALSE)
# 
# 
# # #-----# Read in the matrix generated by the MC sampling:
# # # global.results.readin <- read.csv("matrix of global sensitivity analysis_no migration_N=10000.csv", head=TRUE, stringsAsFactor=FALSE)
# # # global.results <- global.results.readin
# # 
# # 
# # #-----# Global sensitivity analysis results:
# # # 1. Assess relative importance of parameters:
# # # 1.1 Consider the mean genotyep after releasing
# # # Use default parameter of random forest
# # global.rf <- randomForest(x = global.results[,which(names(global.results) %in% list_para)], y = global.results$post_mean,
# #                                   importance = TRUE, na.action = na.omit)
# # 
# # # Output the importance ranking
# # global.importance <- as.data.frame(round(importance(global.rf), 3))
# # global.importance.MSE <- global.importance[order(global.importance$'%IncMSE'),]
# # write.csv(global.importance.MSE, paste0("GSA_VI_Postmean_NoScale_IncMSE_N=", N_sample, "default.csv"))
# # global.importance.NP <- global.importance[order(global.importance$IncNodePurity),]
# # write.csv(global.importance.NP, paste0("GSA_VI_Postmean_NoScale_IncNodePurity_N=", N_sample, "default.csv"))
# # 
# # global.importance.scale <- as.data.frame(round(importance(global.rf, scale = TRUE), 3))
# # global.importance.scale.MSE <- global.importance.scale[order(global.importance.scale$'%IncMSE'),]
# # write.csv(global.importance.scale.MSE, paste0("GSA_VI_Postmean_Scale_IncMSE_N=", N_sample, "default.csv"))
# # global.importance.scale.NP <- global.importance.scale[order(global.importance.scale$IncNodePurity),]
# # write.csv(global.importance.scale.NP, paste0("GSA_VI_Postmean_Scale_IncNodePurity_N=", N_sample, "default.csv"))
# # 
# # # Increase ntree and mtry?
# # # global.rf1 <- randomForest(x = global.results[,which(names(global.results) %in% list_para)], y = global.results$post_mean, ntree = 1000, mtry = 8,
# # #                            importance = TRUE, na.action = na.omit)
# # # global.importance1 <- as.data.frame(round(importance(global.rf1, scale=TRUE), 2))
# # # global.importance1.sort <- global.importance[order(global.importance1$IncNodePurity),]
# # # write.csv(global.importance.sort, paste0("global sensitivity_parameter importance_1_N=", N_sample, ".csv"))
# # 
# # 
# # # 1.2 Consider the relative integrated vector competence
# # global.rf.prop <- randomForest(x = global.results[,which(names(global.results) %in% list_para)], y = global.results$prop,
# #                           importance = TRUE, na.action = na.omit)
# # 
# # # Output the importance ranking
# # global.importance.prop <- as.data.frame(round(importance(global.rf.prop), 3))
# # global.importance.prop.MSE <- global.importance.prop[order(global.importance.prop$'%IncMSE'),]
# # write.csv(global.importance.prop.MSE, paste0("GSA_VI_InteVC_NoScale_IncMSE_N=", N_sample, "default.csv"))
# # global.importance.prop.NP <- global.importance.prop[order(global.importance.prop$IncNodePurity),]
# # write.csv(global.importance.prop.NP, paste0("GSA_VI_InteVC_NoScale_IncNodePurity_N=", N_sample, "default.csv"))
# # 
# # global.importance.prop.scale <- as.data.frame(round(importance(global.rf.prop, scale = TRUE), 3))
# # global.importance.prop.scale.MSE <- global.importance.prop.scale[order(global.importance.prop.scale$'%IncMSE'),]
# # write.csv(global.importance.prop.scale.MSE, paste0("GSA_VI_InteVC_Scale_IncMSE_N=", N_sample, "default.csv"))
# # global.importance.prop.scale.NP <- global.importance.prop.scale[order(global.importance.prop.scale$IncNodePurity),]
# # write.csv(global.importance.prop.scale.NP, paste0("GSA_VI_InteVC_Scale_IncNodePurity_N=", N_sample, "default.csv"))
# # 
# # 
# # 
# # # 2. Build a regression tree
# # # Post-releasing mean VC genotype value
# # subset.post_mean <- global.results[,c(which(names(global.results) %in% list_para),
# #                                       which(names(global.results)=="post_mean"))]
# # tree <- rpart(post_mean ~ ., subset.post_mean)
# # # par(mar=c(5, 4, 4, 2)-1)
# # pdf(paste0("regression tree_N=", N_sample, "_post-releasing mean VC.pdf"), width=12, height=12)
# # plot(tree)
# # text(tree, use.n = TRUE)
# # dev.off()
# # 
# # # Post-releasing integrated VC:
# # subset.prop <- global.results[,c(which(names(global.results) %in% list_para),
# #                                       which(names(global.results)=="prop"))]
# # tree <- rpart(prop ~ ., subset.prop)
# # # par(mar=c(5, 4, 4, 2)-1)
# # pdf(paste0("regression tree_N=", N_sample, "_post-releasing integrated VC.pdf"), width=12, height=12)
# # plot(tree)
# # text(tree, use.n = TRUE)
# # dev.off()
# #
# # 
# # 
# # #-----# 2. Only manageable variables
# # global.results.manageable.default.all <- Global(para.global = para, allvariable = FALSE, N_sample = N_sample, alter=2)
# # global.results.manageable <- global.results.manageable.default.all$global.results
# # 
# # write.csv(global.results.manageable, paste0("matrix of global sensitivity analysis_manageable variables_N=", 
# #                                             global.results.manageable.default.all$N_sample, ".csv"), row.names=FALSE)
# # 
# # 
# # #-----# Analysis of global sensitivity analysis results (only managaable variables):
# # # 1. Assess relative importance of only manageable variables:
# # # 1.1 Consider the mean genotype after releasing
# # # Use default parameter of random forest
# # global.rf.manageable <- randomForest(x = global.results.manageable[,which(names(global.results.manageable) %in% list_para)], 
# #                                      y = global.results.manageable$post_mean,
# #                                      importance = TRUE, na.action = na.omit)
# # 
# # # Output the importance ranking
# # global.importance.manageable <- as.data.frame(round(importance(global.rf.manageable), 3))
# # global.importance.manageable.MSE <- global.importance.manageable[order(global.importance.manageable$'%IncMSE'),]
# # write.csv(global.importance.manageable.MSE, paste0("GSA_manageable_VI_Postmean_NoScale_IncMSE_N=", N_sample, "default.csv"))
# # global.importance.manageable.NP <- global.importance.manageable[order(global.importance.manageable$IncNodePurity),]
# # write.csv(global.importance.manageable.NP, paste0("GSA_manageable_VI_Postmean_NoScale_IncNodePurity_N=", N_sample, "default.csv"))
# # 
# # global.importance.manageable.scale <- as.data.frame(round(importance(global.rf.manageable, scale = TRUE), 3))
# # global.importance.manageable.scale.MSE <- global.importance.manageable.scale[order(global.importance.manageable.scale$'%IncMSE'),]
# # write.csv(global.importance.manageable.scale.MSE, paste0("GSA_manageable_VI_Postmean_Scale_IncMSE_N=", N_sample, "default.csv"))
# # global.importance.manageable.scale.NP <- global.importance.manageable.scale[order(global.importance.manageable.scale$IncNodePurity),]
# # write.csv(global.importance.manageable.scale.NP, paste0("GSA_manageable_VI_Postmean_Scale_IncNodePurity_N=", N_sample, "default.csv"))
# # 
# # 
# # # 1.2 Consider the relative integrated vector competence
# # global.rf.manageable.prop <- randomForest(x = global.results.manageable[,which(names(global.results.manageable) %in% list_para)], y = global.results$prop,
# #                                           importance = TRUE, na.action = na.omit)
# # 
# # # Output the importance ranking
# # global.importance.manageable.prop <- as.data.frame(round(importance(global.rf.manageable.prop), 3))
# # global.importance.manageable.prop.MSE <- global.importance.manageable.prop[order(global.importance.manageable.prop$'%IncMSE'),]
# # write.csv(global.importance.manageable.prop.MSE, paste0("GSA_manageable_VI_InteVC_NoScale_IncMSE_N=", N_sample, "default.csv"))
# # global.importance.manageable.prop.NP <- global.importance.manageable.prop[order(global.importance.manageable.prop$IncNodePurity),]
# # write.csv(global.importance.manageable.prop.NP, paste0("GSA_manageable_VI_InteVC_NoScale_IncNodePurity_N=", N_sample, "default.csv"))
# # 
# # global.importance.manageable.prop.scale <- as.data.frame(round(importance(global.rf.manageable.prop, scale = TRUE), 3))
# # global.importance.manageable.prop.scale.MSE <- global.importance.manageable.prop.scale[order(global.importance.manageable.prop.scale$'%IncMSE'),]
# # write.csv(global.importance.manageable.prop.scale.MSE, paste0("GSA_manageable_VI_InteVC_Scale_IncMSE_N=", N_sample, "default.csv"))
# # global.importance.manageable.prop.scale.NP <- global.importance.manageable.prop.scale[order(global.importance.manageable.prop.scale$IncNodePurity),]
# # write.csv(global.importance.manageable.prop.scale.NP, paste0("GSA_manageable_VI_InteVC_Scale_IncNodePurity_N=", N_sample, "default.csv"))
# # 
# # 
# # 
# # # 2. Build a regression tree
# # # Post-releasing mean VC genotype value
# # subset.manageable.post_mean <- global.results.manageable[,c(which(names(global.results.manageable) %in% list_para), 
# #                                       which(names(global.results.manageable)=="post_mean"))]
# # tree.manageable <- rpart(post_mean ~ ., subset.manageable.post_mean)
# # # par(mar=c(5, 4, 4, 2)-1)
# # pdf(paste0("regression tree_manageable_N=", N_sample, "_post-releasing mean VC.pdf"), width=12, height=12)
# # plot(tree.manageable)
# # text(tree.manageable, use.n = TRUE)
# # dev.off()
# # 
# # # Post-releasing integrated VC:
# # subset.manageable.prop <- global.results.manageable[,c(which(names(global.results.manageable) %in% list_para), 
# #                                  which(names(global.results.manageable)=="prop"))]
# # tree.manageable <- rpart(prop ~ ., subset.manageable.prop)
# # # par(mar=c(5, 4, 4, 2)-1)
# # pdf(paste0("regression tree_manageable_N=", N_sample, "_post-releasing integrated VC.pdf"), width=12, height=12)
# # plot(tree.manageable)
# # text(tree.manageable, use.n = TRUE)
# # dev.off()
# # 
# # 
# 
# 
# ###########################################################################################################
# ###########################################################################################################

# #----------------------- Two variable local sensitivity analysis ---------------------------------
# # The two variables are selected based on the variable importance examined in the above global sensitivity
# # analysis.
# # Except for the two target parameters/variables, all other parameters are the same as in the default settings.
# 
# ### Predictor: post-releasing mean genotype & relative integrated VC after releasing
# 
# # Function to calculate local sensitivity of two variables
# Two_variable <- function(sen1, sen2, sen1.range, sen2.range, para.two, alter){
#   results.mean.two <- matrix(NA, nrow=length(sen1.range), ncol=length(sen2.range))
#   rownames(results.mean.two) <- sen1.range
#   colnames(results.mean.two) <- sen2.range
#   
#   results.integrate.two <- results.mean.two
#   
#   for(j in 1:length(sen1.range)){
#     para.two <- Update(para.two, sen1, sen1.range[j])
#     out.two <- Locals(sen = sen2, sen.range = sen2.range, para = para.two, alter)
#     
#     last <- nrow(out.two$sen_mean)
#     
#     results.mean.two[j,] <- out.two$sen_mean[last,]
#     results.integrate.two[j,] <- (out.two$sen_mean[last,]*out.two$sen_pop[last,]) / (out.two$sen_mean[1,]*out.two$sen_pop[1,])
#   }
#   return(list(mean.mat=results.mean.two, integrate.mat=results.integrate.two))
# }
# Two_variable <- cmpfun(Two_variable)
# 
# #----# Register cores and clusters
# if(getDoParWorkers()==1){
#   cl <- makeCluster(ncore)
#   registerDoParallel(cl)
# }
# getDoParWorkers()
# 
# 
# #-----# Two-variable sensitivity analysis on target variables (chosen based on the importance index from GSA)
# # 1. fm vs. vs
# allvariable.mean <- Two_variable("fm", "vs", fm.range, vs.range, para.two = para, alter)
# write.csv(allvariable.mean$mean.mat, "fm and vs_post-releasing mean.csv", row.names=FALSE)
# 
# # Make a filled contour plot!
# pdf("try.pdf", width=7, height = 6)
# filled.contour(x = fm.range, y = log(vs.range), z = allvariable.mean$mean.mat, nlevels=10,
#                color.palette = colorRampPalette(c("white", "gray10")), las=1, axes=TRUE,
#                xlab="field VC mean", ylab="log selection variance", cex.lab=1.5,
#                plot.axes= {contour(x = fm.range, y = log(vs.range), 
#                                    z = allvariable.mean$mean.mat, nlevels=10,
#                                    drawlabels=FALSE, add=TRUE);
#                            axis(1, cex.axis=1.5); axis(2, cex.axis=1.5)})
# dev.off()
# 
# # 2. rm vs. effective size
# size.enlarge <- seq(from = 0.01, to = 1, by = 0.01)
# allvariable.prop1 <- Two_variable("rm", "size", rm.range, size.enlarge, para.two = para, alter)
# write.csv(allvariable.prop$integrate.mat, "rm and size_proportion of integrated VC.csv", row.names=FALSE)
# 
# # Make a filled contour plot!
# pdf("rm_size.pdf", width=7, height = 6)
# filled.contour(x = rm.range[28:length(rm.range)], y = (size.enlarge * para$fit), 
#                z = allvariable.prop$integrate.mat[28:length(rm.range),], 
#                nlevels=10,
#                color.palette = colorRampPalette(c("white", "gray10")), las=0, axes=TRUE,
#                xlab="Releasing VC mean (rm)", ylab="Effective releasing size (size×fit)", cex.lab=1.5,
#                plot.axes= {contour(x = rm.range[28:length(rm.range)], 
#                                    y = size.enlarge * para$fit,
#                                    z = allvariable.prop$integrate.mat[28:length(rm.range),],
#                                    nlevels=10,
#                                    drawlabels=FALSE, add=TRUE);
#                            axis(1, cex.axis=1.5); axis(2, cex.axis=1.5)})
# dev.off()
