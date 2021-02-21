#' @include API-encoding.R

NULL

priorPredictions <- function(model="LV05-encoding",
                              priors=list(lf="Normal(0.15,0.03)"),
                              samples=4000,
                              cores=1){
  if(model=="LV05-encoding"){
    predicted_effect <- priorpredict(1,priors,samples)
  }
  return(predicted_effect)
}

posteriorPredictions <- function(model="LV05-encoding",
                              posteriors=list(lf=rbeta(100,20,60)),
                              samples=100,
                              cores=1){
  if(model=="LV05-encoding"){
    predicted_effect <- posteriorpredict(1,posteriors,samples)
  }
  return(predicted_effect)
}

estimateParameters <- function(model="LV05-encoding",
                                algorithm="ABC",
                                freeparams=c("lf"),
                                priors=list(lf="Normal(0.15,0.03)"),
                                data=list(effect_gram=rnorm(1000,20,2),
                                          effect_ungram=rnorm(1000,-20,2)),
                                samples=100,
                                warmup=50,
                                cores=1,
                                chains=1){
  if(model=="LV05-encoding"){
    posterior_samples <- estimate(algorithm,freeparams,priors,data,samples,warmup,cores,chains)
  }
  return(posterior_samples)
}


bayesFactor <- function(fit.m1,fit.m2){
  ML_m1 <- mean(fit.m1$weight/fit.m1$imp_density)
  ML_m2 <- mean(fit.m2$weight/fit.m2$imp_density)
  return(ML_m1/ML_m2)
}
