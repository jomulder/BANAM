#'  @title Covariate data frame for the Alabama voter turnout data
#'
#'  @description This data.frame includes logarithmized data on population casting votes, of the 
#'  population that is 25 years and older who completed 12th grade or higher 
#'  education, on the number of owner-occupied housing units, and of aggregate 
#'  income per county in Alabama, US. This data.frame can be specified as covariate 
#'  matrix in the network autocorrelation model to predict the voter turnout in 
#'  Alabama.
#'
#'  \tabular{lll}{
#'  \strong{pop_eligible} \tab \code{numeric} \tab Logarithm of the population casting votes \cr
#'  \strong{pop_college} \tab \code{numeric} \tab Logarithm of the population that is 25 years and older who completed 12th grade or higher education \cr
#'  \strong{homeownership} \tab \code{numeric} \tab Logarithm of the number of owner-occupied housing units \cr
#'  \strong{income} \tab \code{numeric} \tab Logarithmized aggregate income per county \cr
#'  }
#'  @docType data
#'  @keywords datasets
#'  @name X_votes
#'  @usage data(X_votes)
#'  @references Pace, R. K. and R. Barry. 1997. Quick Computation of 
#'  Spatial Autoregressive Estimators. Geographical Analysis, 29, 232-47. 
#'  Data can be downloaded from 
#'  http://www.spatial-econometrics.com/html/jplv7.zip
#'  @format A data.frame with 67 rows and 4 columns
##
