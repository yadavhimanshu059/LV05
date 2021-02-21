\name{priorPredictions}
\alias{priorPredictions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate prior predictions
}
\description{
The function generates prior predictive distribution from the Lewis and Vasishth (2005) model.
}
\usage{
priorPredictions(model="LV05-encoding",
                  priors=list(lf="Normal(0.15,0.03)"),
                  samples=4000,
                  cores=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
    Version of the model
  }
  \item{priors}{
    Define priors for each parameter

    For example, priors=list(lf="Normal(0.15,0.03)",cw="Uniform(1,4)")

    Default: list(lf="Normal(0.15,0.03)")
  }
  \item{samples}{
    Total number of prior predicted datapoints
  }
  \item{cores}{
    Number of cores to be used

    Preset to 1 for now.
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

## Priors on latency factor, lf and percolation rate, prob
priors = list(lf = "Normal(mean=.15, sd=.03, lb=.1, ub=1)",
                 prob = "Normal(mean=0, sd=.25, lb= 0, ub=1)")

predictions <- priorPredictions(priors=priors,
                                samples=2000)

## The function is currently defined as
function (model = "LV05-encoding", priors = list(lf = "Normal(0.15,0.03)"),
    samples = 4000, cores = 1)
{
    if (model == "LV05-encoding") {
        predicted_effect <- priorpredict(1, priors, samples)
    }
    return(predicted_effect)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
