library(ggplot2)
library(ggpubr)
library(ANTs)
library(sjPlot)
###########################
#### Basic simulation######
###########################
#' @param N an integer indicating the total number of individuals
#' @param Obs.mean.bias amount of observation bias between sex (in favor of females)
#' @param diff.soc.pheno amount of interaction differences per observation between sex (in favor of females)
#' @param N.obs.pheno the base number of observation per individuals (will be define for males)
#' @param N.alters.per.obs the number of alters per observation (will be define for males)
#' @param interactions.biases the percent of interaction misses within females observations
sim <- function(N=30, Obs.mean.bias = 0, obs.sd.bias = 1, diff.soc.pheno = 50, N.obs.pheno = 10, N.alters.per.obs = 5,
                interactions.biases = -50, interactions.biases.sd = 0, females.interactions.biases = "F", print = T){
  ids = 1 : N
  sex = sample(c("F", "M"), N, replace = T, prob = c(0.5,0.5))
  
  # observations biases--------------
  obs.m = N.obs.pheno
  obs.f = obs.m + ((obs.m*Obs.mean.bias)/100) # add the amount of biases define by "Obs.mean.bias"
  
  # Interactions biases--------------
  N.m = N.obs.pheno
  N.f = N.m + ((N.m*diff.soc.pheno)/100) # add the amount of sociality differences define by "diff.soc.pheno"
  
  Result = NULL
  sampling.effort = NULL
  tmp = data.frame(ids,sex)
  tmp$interactions = NA
  
  # Observations biases--------------
  for (a in 1:N) {
    ego.obs = NULL
    ego = a
    ego.sex = sex[a]
    if(ego.sex == "M"){
      ego.obs = rnorm(1, mean = obs.m, sd =  abs((obs.m*obs.sd.bias)/100))
      ego.obs = round(ego.obs)
      sampling.effort = rbind(sampling.effort, data.frame(ego, "sampling.effort" = ego.obs, "sex" = ego.sex))
      if(ego.obs < 1 ){next()}
      for (b in 1:ego.obs) {
        cat("Processing individual ", a, " observation #", b, "\r")
        alters = sample(ids[!ids %in% ego], rnorm(1, mean = N.m, sd = 0), replace = T)
        Result = rbind(Result, data.frame(ego, "sex" = ego.sex, "focal" = b, "alters" = alters))
        
        tmp$interactions[a] = sum(tmp$interactions[a], length(alters), na.rm = T)
        
      }
      
    }else{
      ego.obs = rnorm(1, mean = obs.f, sd = abs((obs.f*obs.sd.bias)/100))
      ego.obs = round(ego.obs)
      sampling.effort = rbind(sampling.effort, data.frame(ego, "sampling.effort" = ego.obs, "sex" = ego.sex))
      if(ego.obs < 1 ){next()}
      for (b in 1:ego.obs) {
        cat("Processing individual ", a, " observation #", b, "\r")
        alters = sample(ids[!ids %in% ego], round(rnorm(1, mean = N.f, sd = 0)), replace = T)
        Result = rbind(Result, data.frame( ego, "sex" = ego.sex, "focal" = b, "alters" = alters))
        tmp$interactions[a] = sum(tmp$interactions[a], length(alters), na.rm = T)
        
      }
    }
  }
  tmp2 = split(Result, Result$ego)
  p1 = ggplot(data = NULL, aes(x = sex, y = (unlist(lapply(tmp2,function(x){nrow(x)}))), group = sex, color = sex))+geom_boxplot()+geom_point()
  
  
  # Interactions biases--------------
  if(interactions.biases != 0){
    tmp = Result[Result$sex %in% females.interactions.biases,]# selection of phenotype
    tmp$ctrl = paste(tmp$ego, tmp$focal)
    tmp = split(tmp, tmp$ctrl)
    
    tmp = lapply(tmp, function(x, interactions.biases){
      t = abs(round((nrow(x) * interactions.biases )/100))
      if(t == 0){}else{
        to.remove = round(rnorm(1,mean = t, sd = abs((t*interactions.biases.sd)/100)))
        x = x[-sample(1:nrow(x), to.remove, replace = F),]
      }
      x
      
    }, interactions.biases = interactions.biases)
    tmp = do.call("rbind", tmp)
    tmp = tmp[,-ncol(tmp)]
    Result = rbind(Result[!Result$sex %in% females.interactions.biases,],tmp)
  }
  
  tmp3 = split(Result, Result$ego)
  p2 = ggplot(data = NULL, aes(x = sex, y = (unlist(lapply(tmp3,function(x){nrow(x)}))), group = sex, color = sex))+geom_boxplot()+geom_point()
  
  if(print){
    print(ggarrange(p1,p2,nrow = 1, ncol = 2, common.legend = T))
  }
  return(list("interactions" = Result, "informations" = sampling.effort))
}

#' @name Test test simulation 
test.sim <- function(sim, correction = T){
  ## Observations per individuals
  tmp = sim$interactions
  tmp = split(tmp, tmp$ego)
  tmp2 = sim$informations
  tmp2$obs = (unlist(lapply(tmp,function(x){length(unique(x$focal))})))
  p1 = ggplot(tmp2, aes(x = sex, y = obs, group = sex, color = sex))+geom_boxplot()+geom_point()+ylab('Number of observations')
  
  ## Number of interactions within observations
  tmp2$Int.mean = (unlist(lapply(tmp,function(x){
    t = split(x, x$focal)
    mean(unlist(lapply(t, nrow)))
  })))
  p2 = ggplot(tmp2, aes(x = sex, y = Int.mean, group = sex, color = sex))+geom_boxplot()+geom_point()+ylab('Frequencies in interactions')
  
  ## Strength
  if(correction){
    m = df.to.mat(sim$interactions, 1,4,num.ids = T, tobs = sim$informations$sampling.effort)# correction
  }else{
    m = df.to.mat(sim$interactions, 1,4,num.ids = T)# No correction
  }
  
  tmp2 = met.strength(m, df = tmp2, dfid = "ego")
  p3 = ggplot(tmp2, aes(x = sex, y = strength, group = sex, color = sex))+geom_boxplot()+geom_point()
  
  print(
    ggarrange(p1,p2,p3, nrow = 1, ncol = 3, common.legend = T)
  )
  
  return(tmp2)
}

# No biases
test = sim(N = 30, Obs.mean.bias = 0,obs.sd.bias = 5, diff.soc.pheno = 0, N.obs.pheno= 50, N.alters.per.obs = 10, interactions.biases.sd = 0)
test.sim(test)
summary(lm(data= test.sim(test), formula = strength~sex))# No significant effect

# No biases but sd variation in observations
test = sim(N = 30, Obs.mean.bias = 0, diff.soc.pheno = 0, N.obs.pheno= 50, N.alters.per.obs = 10, obs.sd.bias = 5, interactions.biases.sd = 0)
test.sim(test)
summary(lm(data= test.sim(test), formula = strength~sex))# No significant effect

# No biases but sd variation in interactions
test = sim(N = 30, Obs.mean.bias = 0, diff.soc.pheno = 0, N.obs.pheno= 50, N.alters.per.obs = 10, obs.sd.bias = 0, interactions.biases.sd = 5)
test.sim(test)
summary(lm(data= test.sim(test), formula = strength~sex))# Only sd variation in number of observations or interaction biases can be enough to lead to false results

# Sampling effort biases
test = sim(N = 30, Obs.mean.bias = -50, diff.soc.pheno = 0, N.obs.pheno= 50, N.alters.per.obs = 10, obs.sd.bias = 5, interactions.biases.sd = 0,
           interactions.biases = 0)
test.sim(test)
summary(lm(data= test.sim(test), formula = strength~sex))# Correction work most of time for observation biases

# Interactions biases
test = sim(N = 30, Obs.mean.bias = 0, diff.soc.pheno = 0, N.obs.pheno= 50, N.alters.per.obs = 10, obs.sd.bias = 0, interactions.biases.sd = 5,
           interactions.biases = 5)
test.sim(test)
summary(lm(data= test.sim(test), formula = strength~sex))# But not when there is interaction biases

# Sampling effort and interactions biases
test = sim(N = 30, Obs.mean.bias = -50, diff.soc.pheno = 0, N.obs.pheno= 50, N.alters.per.obs = 10, obs.sd.bias = 0, interactions.biases.sd = 0,
           interactions.biases = 10)
test.sim(test)
summary(lm(data= test.sim(test), formula = strength~sex))


