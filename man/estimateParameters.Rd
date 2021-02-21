\name{estimateParameters}
\alias{estimateParameters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate model parameters
}
\description{
The function estimates posterior distributions for the free parameters of the Lewis and Vasishth (2005) model.
}
\usage{
estimateParameters(model = "LV05-encoding",
                    algorithm = "ABC",
                    freeparams = c("lf"),
                    priors = list(lf = "Normal(0.15,0.03)"),
                    data = list(effect_gram = rnorm(1000, 20, 2), effect_ungram = rnorm(1000, -20, 2)),
                    samples = 100, warmup = 50, cores = 1, chains = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
    Version of the model.
}
  \item{algorithm}{
    Algorithm to be used for parameter estimation.
}
  \item{freeparams}{
    List of free parameters to be estimated.

    The current version allows - latency factor "lf", Cue weighting "cw", and feature percolation rate "prob"
}
  \item{priors}{
    Prior distributions on free parameters

    Priors should be defined in this format: list(lf="Normal(.2,1)",cw="Beta(2,6)")

    The current version allows only Normal, Beta and Uniform priors.
    You can define lower and upper bound on the prior distiributions like this "Normal(.2,1,lb=0.1,ub=1)".
}
  \item{data}{
    Data used for estimating the parameters.

    Data should be defined in this format: list(effect_gram,effect_ungram), where effect_gram and effect_ungram are vectors containing inference effect distribution in grammatical and ungrammatical conditions respectively
}
  \item{samples}{
    Total number of samples in each chain
}
  \item{warmup}{
    Number of burn-in samples in each chain
}
  \item{cores}{
    Number of cores

    Multiple cores can be used if your are using Linux system.
}
  \item{chains}{
    Number of markov chains

    Make sure number of chains are equal to number of cores
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
Himanshu Yadav (hyadav@uni-potsdam.de)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## Estimate latency factor parameter
data = list(effect_gram=rnorm(1000,22,5),
            effect_ungram=rnorm(1000,-22,5))

fit.model <- estimateParameters(data=data,
                              freeparams=c("lf"),
                              priors=list(lf="Normal(0.15,0.03,lb=0.10,ub=1)"),
                              samples=2000,warmup=1000,chains=4,cores=4)

## Estimate both latency factor and percolation rate simultaneously
fit.model <- estimateParameters(data = data,
                              freeparams=c("lf","prob"),
                              priors=list(
                                      lf="Normal(mean=0.15,sd=0.03,lb=0.1,ub=1)",
                                      prob="Normal(mean=0,sd=0.2,lb=0,ub=1)"),
                              samples=2000,warmup=1000,chains=4,cores=4)
