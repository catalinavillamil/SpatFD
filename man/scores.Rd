\name{scores}
\alias{scores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ scores
     %%  ~~function to do ... ~~
}
\description{
Scalar associated with functional principal components   that are a scalar spatial random field     %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
     scores(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
     \item{X}{ An object of class SpatFD or scores_pred, which are output of functions SpatFD or KS_scores.}

     %%     ~~Describe \code{x} here~~
}


\value{
     \item{scores}{Scalar associated with functional principal components that are a scalar spatial random field.}
     %%  ~Describe the value returned
     %%  If it is a LIST, use
     %%  \item{comp1 }{Description of 'comp1'}
     %%  \item{comp2 }{Description of 'comp2'}
     %% ...
}

\author{ Angie Villamil, Martha Bohorquez
     %%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

     \seealso{\code{\link{SpatFD}},\code{\link{KS_scores}}
          %% ~~objects to See Also as \code{\link{help}}, ~~~
     }
\examples{
     ##---- Should be DIRECTLY executable !! ----
     ##-- ==>  Define data, use random,
     ##--	or do  help(data=index)  for the standard data sets.

     ## The function is currently defined as
     function (x)
     {
     }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

