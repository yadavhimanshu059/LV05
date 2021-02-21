\name{posteriorPredictions}
\alias{posteriorPredictions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate posterior predictions
}
\description{
The function generates posterior predictive distribution from the Lewis and Vasishth (2005) model
}
\usage{
posteriorPredictions(model = "LV05-encoding",
                      posteriors = list(lf = rbeta(100, 20, 60)),
                      samples = 100, cores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
    Version of the model
}
  \item{posteriors}{
    Posterior distributions to be used to generate predictions

    Posterior distributions should be defined as list of vectors, like this
    list(lf=lfpost, cw=cwpost) where lfpost and cwpost are vectors containing posterior distribution of latency factor and cue weighting respectively
}
  \item{samples}{
    Number of datapoints in the posterior distribution of each parameter
}
  \item{cores}{
    Number of cores to be used
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{

}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## Let's say we have estimated posterior distribution for latency factor

data = list(effect_gram=rnorm(1000,22,5),
            effect_ungram=rnorm(1000,-22,5))

fit.model <- estimateParameters(data=data,
                              freeparams=c("lf"),
                              priors=list(lf="Normal(0.15,0.03,lb=0.10,ub=1)"),
                              samples=2000,warmup=1000,chains=4,cores=4)

## We can extract posterior distribution for latency factor
lfpost <- fit.model$latency_f

## Now we can generate posterior predictive distribution
df.predictions <- posteriorPredictions(posteriors = list(lf=lfpost),
                                        samples = length(lf_post))
