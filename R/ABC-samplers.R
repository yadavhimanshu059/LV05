sampler_avgeff_meanSS_PRC <- function(freeparams,priors,data,samples,warmup,cores,chains){
  priors <- priors
  freeparams <- freeparams
  ybar_gram <- mean(data$effect_gram)
  ybar_ungram <- -mean(data$effect_ungram)
  sdx <- sd(data$effect_gram)
  tol <- 10
  if(samples>=1000){
    niter <- 8
    npart <- samples/8
    burniter <- (warmup/samples)*niter
  }else{
    niter <- 2
    npart <- samples/2
    burniter <- (warmup/samples)*niter
  }

  sampler <- function(chain_id){
    reject <- 0
    df.pool <- data.frame(matrix(ncol = 11, nrow = 0))
    colnames(df.pool) <- c("pool_id","sample_id","latency_f","cue_w","perc_rate","prior.prob","weight","likelihood","imp_density","xsim.gram","xsim.ungram")
    for(t in 1:niter){
      if(t==1){
        for(i in 1:npart){
          tolw <- tol/(t+4)
          k <- 1
          while(k==1){
            if("lf" %in% freeparams){
              latency_f <- runif(1,0.1,0.6)
            }else{
              latency_f <- 0.15
            }
            if("cw" %in% freeparams){
              cue_w <- runif(1,1,4)
            }else{
              cue_w <- 1
            }
            if("prob" %in% freeparams){
              prob <- runif(1,0,0.6)
            }else{
              prob <- 0
            }
            proposal <- data.frame(latency_f,cue_w,prob)
            ## get generated effect for ungrammatical conditions:
            generated_effect <- iterate_lf(proposal)$Effect
            xsim_gram <- generated_effect[1]
            xsim_ungram <- generated_effect[2]
            disc_gram <- abs(xsim_gram-ybar_gram)
            disc_ungram <- abs(xsim_ungram-ybar_ungram)
            if(disc_gram<100&disc_ungram<100){
              if(all(c("lf","cw","prob") %in% freeparams)){
                prior_i <- priorfn(priors$lf,proposal=latency_f)$dens*priorfn(priors$prob,proposal=prob)$dens*priorfn(priors$cw,proposal=cue_w)$dens
              }else if(all(c("lf","cw") %in% freeparams)){
                prior_i <- priorfn(priors$lf,proposal=latency_f)$dens*priorfn(priors$cw,proposal=cue_w)$dens
              }else if(all(c("lf","prob") %in% freeparams)){
                prior_i <- priorfn(priors$lf,proposal=latency_f)$dens*priorfn(priors$prob,proposal=prob)$dens
              }else{
                prior_i <- priorfn(priors$lf,proposal=latency_f)$dens
              }
              weight_i <- dnorm(disc_gram,0,tolw)*dnorm(disc_ungram,0,tolw)*prior_i
              likelihood_i <- dnorm(disc_gram,0,tolw)*dnorm(disc_ungram,0,tolw)
              imp_density_i <- prior_i
              df.pool[length(df.pool$pool_id)+1,]<-c(t,i,latency_f,cue_w,prob,prior_i,weight_i,likelihood_i,imp_density_i,xsim_gram,xsim_ungram)
              k <- k+1
            }else{
              reject <- reject+1
            }
          }
        }
      }else{
        prev.pool.lf <- subset(df.pool,pool_id==(t-1))$latency_f
        prev.weight <- subset(df.pool,pool_id==(t-1))$weight
        prev.pool.prob <- subset(df.pool,pool_id==(t-1))$perc_rate
        prev.pool.cw <- subset(df.pool,pool_id==(t-1))$cue_w
        for(i in 1:npart){
          tolw <- tol/(t+4)
          samp.lf <- sample(prev.pool.lf,size=1,prob = prev.weight)
          laf <- -1
          while(laf<0.1|laf>1){
            latency_f <- rnorm(1,samp.lf,0.03)
            laf <- latency_f
          }
          if(all(c("lf","cw","prob") %in% freeparams)){
            samp.prob <- prev.pool.prob[which(prev.pool.lf==samp.lf)]
            samp.cw <- prev.pool.cw[which(prev.pool.lf==samp.lf)]
            cuw <- -1
            while(cuw<0|cuw>0.8){
              prob <- rnorm(1,samp.prob,0.03)
              cuw <- prob
            }
            cuw <- -1
            while(cuw<1|cuw>4){
              cue_w <- rnorm(1,samp.cw,0.03)
              cuw <- cue_w
            }
          }else if(all(c("lf","cw") %in% freeparams)){
            samp.cw <- prev.pool.cw[which(prev.pool.lf==samp.lf)]
            cuw <- -1
            while(cuw<1|cuw>4){
              cue_w <- rnorm(1,samp.cw,0.03)
              cuw <- cue_w
            }
            prob <- 0
          }else if(all(c("lf","prob") %in% freeparams)){
            samp.prob <- prev.pool.prob[which(prev.pool.lf==samp.lf)]
            cuw <- -1
            while(cuw<0|cuw>0.8){
              prob <- rnorm(1,samp.prob,0.03)
              cuw <- prob
            }
            cue_w <- 1
          }else{
            prob <- 0
            cue_w <- 1
          }
          # Sample a proposal from the previous pool
          proposal <- data.frame(latency_f,cue_w,prob)
          ## get generated effect for ungrammatical conditions:
          generated_effect <- iterate_lf(proposal)$Effect
          xsim_gram <- generated_effect[1]
          xsim_ungram <- generated_effect[2]
          disc_gram <- abs(xsim_gram-ybar_gram)
          disc_ungram <- abs(xsim_ungram-ybar_ungram)
          if(all(c("lf","cw","prob") %in% freeparams)){
            prior_i <- priorfn(priors$lf,proposal=latency_f)$dens*priorfn(priors$prob,proposal=prob)$dens*priorfn(priors$cw,proposal=cue_w)$dens
            imp_density_i <- dnorm(latency_f,samp.lf,0.03)*dnorm(prob,samp.prob,0.03)*dnorm(cue_w,samp.cw,0.03)
          }else if(all(c("lf","cw") %in% freeparams)){
            prior_i <- priorfn(priors$lf,proposal=latency_f)$dens*priorfn(priors$cw,proposal=cue_w)$dens
            imp_density_i <- dnorm(latency_f,samp.lf,0.03)*dnorm(cue_w,samp.cw,0.03)
          }else if(all(c("lf","prob") %in% freeparams)){
            prior_i <- priorfn(priors$lf,proposal=latency_f)$dens*priorfn(priors$prob,proposal=prob)$dens
            imp_density_i <- dnorm(latency_f,samp.lf,0.03)*dnorm(prob,samp.prob,0.03)
          }else{
            prior_i <- priorfn(priors$lf,proposal=latency_f)$dens
            imp_density_i <- dnorm(latency_f,samp.lf,0.03)
          }
          weight_i <- dnorm(disc_gram,0,tolw)*dnorm(disc_ungram,0,tolw)*prior_i
          likelihood_i <- dnorm(disc_gram,0,tolw)*dnorm(disc_ungram,0,tolw)
          df.pool[length(df.pool$pool_id)+1,]<-c(t,i,latency_f,cue_w,prob,prior_i,weight_i,likelihood_i,imp_density_i,xsim_gram,xsim_ungram)
        }
      }
    }
    df.pool$chain_id <- rep(chain_id,length(df.pool$pool_id))
    df.pool$rejected <- reject
    df.pool
  }

  if(.Platform$OS.type == "unix"){
    df.posterior <- mclapply(1:chains, sampler,mc.cores = cores) %>% bind_rows()
  }else{
    df.posterior <- data.frame(matrix(ncol = 11, nrow = 0))
    colnames(df.posterior) <- c("pool_id","sample_id","latency_f","cue_w","perc_rate","prior.prob","weight","likelihood","imp_density","xsim.gram","xsim.ungram")
    for(i in 1:chains){
      old <- Sys.time()
      df.posterior <- rbind(df.posterior, sampler(i))
      new <- Sys.time() - old
      print(paste("Sampling for chain ",as.character(i)," completed!"))
      print(new)
    }
  }
  df.posterior <- subset(as.data.frame(df.posterior),pool_id>burniter)
  return(df.posterior)
}
