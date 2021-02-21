#' @include interACT-encoding.R

NULL

library(dplyr)
library(tidyr)
library(ggplot2)
library(parallel)
library(LaplacesDemon)
library(mvtnorm)
library(truncnorm)


## some helper functions:
rmsd <- function (obs, pred) {
  sqrt(mean((obs - pred)^2, na.rm = TRUE))
}

compute_int_means <- function(d){
  int <- select(filter(d, Distractor=="Match"), -Condition, -Distractor)
  dim(int)
  int$int <- filter(d, Distractor=="Match")$latency - filter(d, Distractor=="Mismatch")$latency
  #
  # means <- group_by(int, Set, Target, lf, ans, mas, mp, rth, bll, lp, ldp, blc, dbl, ndistr) %>% summarise(Effect=mean(int), SE=sd(int)/sqrt(length(int))) %>% ungroup() %>% mutate(lower=Effect-SE, upper=Effect+SE)
  # means
  means <- group_by(int, Set, Target, lf, ans, mas, mp, rth, bll, psc, pic, qcf, qco, cuesim, tprom, dprom, lp, ldp, blc, dbl, ndistr, cueweighting) %>% summarise(Effect=mean(int), SE=sd(int)/sqrt(length(int)), Sji_neg=sum(Sji_neg)) %>% ungroup() %>% mutate(lower=Effect-SE, upper=Effect+SE)
  means
}

convert2log <- function(x){
  ifelse(x>=1, log(x), ifelse(x<=-1, -log(abs(x)), 0))
}

convert2log10 <- function(x){
  x <- ifelse(x>-1 & x<1, 0, x)
  x <- ifelse(x<=-1, -log10(abs(x)), x)
  x <- ifelse(x>=1, log10(abs(x)), x)
}


# Introduction


## Basic engine for generating predictions

#' @printcounts<-FALSE

iterate_lf <- function(values,iterations=1000){
  ## values is a dataframe contaning either 1 row of parameters (a combination of latency   factor and cue-weighting parameter, or multiple rows i.e., multiple combinations of two parameters)
  ## iterations is the number of iterations for that given value.
  ## We need multiple iterations as noise is non-zero and there will be some
  ## variability due to noise.
  le        <<- 1       # latency exponent
  rth 	    <<- -1.5    # retrieval threshold
  bll       <<- 0.5     # decay parameter
  ans       <<- 0.2    # activation noise
  mas       <<- 1       # maximum associative strength
  mp        <<- 1       # mismatch penalty
  ga        <<- 1       # goal source activation
  rand_time <<- 3       # latency variability
  blc       <<- 0       # base-level activation
  lp        <<- 1     # default time since last presentation (msec)
  ldp       <<- 1       # last distractor presentation (msec)
  ndistr    <<- 1				# number of distractors
  dbl       <<- 0       # distractor base-level
  normalizeWeights <<- TRUE
  qcf 		  <<- 1 			# match quality correction factor
  qco 		  <<- -2*rth 			# match quality correction offset
  psc       <<- 1       # prominence scaling constant C1
  pic       <<- 0       # prominence scaling constant C2
  tprom     <<- 0       # target prominence
  dprom     <<- 0       # distractor prominence
  # cl 				<<- 0 			# cue confusion level (0-100)
  cuesim    <<- -1    # cue-feature similarity [-1..0]
  #
  VERBOSE <<- FALSE

  maxset <- 0
  means <- NULL
  for(v in 1:nrow(values)){
    lf <<- values[v,]$latency_f
    cueweighting <<- values[v,]$cue_w
    perc <<- values[v,]$prob
    parameters <- list(lf,le,rth,bll,ans,mas,mp,ga,rand_time,lp,blc,ldp,dbl,ndistr,cueweighting,psc,pic,qcf, qco,cuesim,tprom,dprom)
    paramnames <- c("lf","le","rth","bll","ans","mas","mp","ga","rand_time","lp","blc","ldp","dbl","ndistr","cueweighting","psc","pic","qcf", "qco" ,"cuesim","tprom","dprom")
    names(parameters) <- paramnames
    pmatr <- create_param_matrix(perc,distortion(perc), parameters, iterations)
    results <- run(pmatr)
    means2 <- compute_int_means(results)
    means2$Set <- means2$Set+maxset
    means <- bind_rows(means, means2)
  }
  means
}




priorfn <- function(priorname,proposal=NA){
  priorname <- gsub(" ", "", priorname)
  priorname <- gsub("mean=", "", priorname)
  priorname <- gsub("sd=", "", priorname)
  priorname <- gsub("lb=", "", priorname)
  priorname <- gsub("ub=", "", priorname)
  p1 <- as.numeric(strsplit(strsplit(strsplit(strsplit(priorname,"[(]")[[1]]," ")[[2]],",")[[1]]," ")[[1]])
  p2 <- as.numeric(strsplit(strsplit(strsplit(strsplit(strsplit(priorname,"[(]")[[1]]," ")[[2]],",")[[1]]," ")[[2]],"[)]")[[1]])
  if(length(strsplit(strsplit(strsplit(priorname,"[(]")[[1]]," ")[[2]],",")[[1]])>2){
    lb <- as.numeric(strsplit(strsplit(strsplit(strsplit(strsplit(priorname,"[(]")[[1]]," ")[[2]],",")[[1]]," ")[[3]],"[)]")[[1]])
    ub <- as.numeric(strsplit(strsplit(strsplit(strsplit(strsplit(priorname,"[(]")[[1]]," ")[[2]],",")[[1]]," ")[[4]],"[)]")[[1]])
  }else{
    lb <- -Inf
    ub <- Inf
  }
  if(strsplit(strsplit(priorname,"[(]")[[1]]," ")[[1]]=="Normal"){
    priorsim <- rtruncnorm(1,a=lb,b=ub,mean=p1,sd=p2)
    priordens <- dtruncnorm(proposal,a=lb,b=ub,mean=p1,sd=p2)
  }else if(strsplit(strsplit(priorname,"[(]")[[1]]," ")[[1]]=="Beta"){
    priorsim <- rbeta(1,p1,p2)
    priordens <- dbeta(proposal,p1,p2)
  }else if(strsplit(strsplit(priorname,"[(]")[[1]]," ")[[1]]=="Uniform"){
    priorsim <- runif(1,p1,p2)
    priordens <- dunif(proposal,p1,p2)
  }else{
    return(print("Error: Only Normal, Beta, and Uniform priors are supported"))
    priordens <- 1
  }
  return(list(sim=priorsim,dens=priordens,p1=p1,p2=p2,lb=lb,ub=ub))
}

#' ## Prior predictions

priorpredict <- function(chain_id,priors,samples){
  nsamples <- samples
  df.pool <- data.frame(matrix(ncol = 6, nrow = nsamples))
  colnames(df.pool) <- c("sample_id","latency_f","cue_w","perc","sim_gram","sim_ungram")
  for(i in 1:nsamples){
    if("lf" %in% names(priors)){
      latency_f <- priorfn(priors$lf)$sim
    }else{
      latency_f <- 0.15
    }
    if("cw" %in% names(priors)){
      cue_w <- priorfn(priors$cw)$sim
    }else{
      cue_w <- 1
    }
    if("prob" %in% names(priors)){
      prob <- priorfn(priors$prob)$sim
    }else{
      prob <- 0
    }
    proposal <- data.frame(latency_f,cue_w,prob)
    ## get generated effect for ungrammatical conditions:
    simx <- iterate_lf(proposal)$Effect
    df.pool[i,]<-c(i,latency_f,cue_w,prob,simx[1],simx[2])
  }
  df.pool$chain_id <- chain_id
  return(df.pool)
}


#' ## Posterior predictions

posteriorpredict <- function(chain_id,posteriors,samples){
  nsamples <- samples
  df.pool <- data.frame(matrix(ncol = 6, nrow = nsamples))
  colnames(df.pool) <- c("sample_id","latency_f","cue_w","perc","sim_gram","sim_ungram")
  for(i in 1:nsamples){
    if("lf" %in% names(posteriors)){
      latency_f <- posteriors$lf[i]
    }else{
      latency_f <- 0.15
    }
    if("cw" %in% names(posteriors)){
      cue_w <- posteriors$cw[i]
    }else{
      cue_w <- 1
    }
    if("prob" %in% names(posteriors)){
      prob <- posteriors$prob[i]
    }else{
      prob <- 0
    }
    proposal <- data.frame(latency_f,cue_w,prob)
    ## get generated effect for ungrammatical conditions:
    simx <- iterate_lf(proposal)$Effect
    df.pool[i,]<-c(i,latency_f,cue_w,prob,simx[1],simx[2])
  }
  df.pool$chain_id <- chain_id
  df.pool
}


#' ##Parameter estimation function

estimate <- function(algorithm,freeparams,priors,data,samples,warmup,cores,chains){
  if(algorithm=="ABC"){
    df.posterior <- sampler_avgeff_meanSS_PRC(freeparams,priors,data,samples,warmup,cores,chains)
    return(df.posterior)
  }else{
    print("Algorithm not available!")
  }
}

