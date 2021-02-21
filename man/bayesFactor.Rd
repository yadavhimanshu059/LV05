\name{bayesFactor}
\alias{bayesFactor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate Bayes factor
}
\description{
The function computes Bayes factor for the models fitted using estimateParameter function in the LV05 package.
}
\usage{
bayesFactor(fit.m1, fit.m2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit.m1}{
    Fitted model 1
}
  \item{fit.m2}{
    Fitted model 2
}
}
\details{
Estimated Bayes factor value indicates evidence in favor of model 1 over model 2
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Data
data = list(effect_gram=rnorm(1000,22,5),
            effect_ungram=rnorm(1000,-22,5))

## Fit model 1 (retrieval-only model)
fit.model1 <- estimateParameters(data=data,
                              freeparams=c("lf"),
                              priors=list(lf="Normal(0.15,0.03,lb=0.10,ub=1)"),
                              samples=2000,warmup=1000,chains=4,cores=4)

## Fit model 2 (encoding-plus-retrieval model)
fit.model2 <- estimateParameters(data = data,
                              freeparams=c("lf","prob"),
                              priors=list(
                                      lf="Normal(mean=0.15,sd=0.03,lb=0.1,ub=1)",
                                      prob="Normal(mean=0,sd=0.2,lb=0,ub=1)"),
                              samples=2000,warmup=1000,chains=4,cores=4)

## Estimate evidence in favor of model 1 over model 2
bf <- bayesFactor(fit.model1,fit.model2)

## The function is currently defined as
function (fit.m1, fit.m2)
{
    ML_m1 <- mean(fit.m1$weight/fit.m1$imp_density)
    ML_m2 <- mean(fit.m2$weight/fit.m2$imp_density)
    return(ML_m1/ML_m2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
