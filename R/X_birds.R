#'  @title Covariate data frame for the threatened bird data
#'
#'  @description Generated data from a study on spatial effect on country's 
#'  threatened bird rates by McPherson and Nieswiadomy (2005). The data includes 
#'  the percentage of threatened birds in 113 countries in the year 2000
#'  along with 10  other predictors, which the authors employed in
#'  their analysis.
#'
#' \tabular{lll}{
#' \strong{intercept} \tab \code{integer} \tab vector of ones to model intercept \cr
#' \strong{ISLAND} \tab \code{integer} \tab Island: 1 = Yes, 0 = No \cr
#' \strong{BIRDENPR} \tab \code{numeric} \tab percentage of endemic birds \cr
#' \strong{POPD8100} \tab \code{numeric} \tab persons per square kilometer, 1981-2000 average \cr
#' \strong{PPP8100} \tab \code{numeric} \tab per capita income in 1995$ purchasing power parity, 1981-2000 average \cr
#' \strong{PPP81002} \tab \code{numeric} \tab Squared per capita income in 1995$ purchasing power parity, 1981-2000 average \cr
#' \strong{POLC8100} \tab \code{numeric} \tab political rights and civil liberties, 1981-2000 average \cr
#' \strong{DEMO8100} \tab \code{numeric} \tab antigovernment demonstrations, per year, 1981-2000 average \cr
#' \strong{CIVIL} \tab \code{integer} \tab Civil law: 1 = Yes, 0 = No \cr
#' \strong{MUSLIM} \tab \code{integer} \tab Muslim law: 1 = Yes, 0 = No \cr
#' \strong{COMMUN} \tab \code{integer} \tab Communist law: 1 = Yes, 0 = No \cr
#' }
#'  @docType data
#'  @keywords datasets
#'  @name X_birds
#'  @usage data(X_birds)
#'  @references McPherson, M. A., & Nieswiadomy, M. L. (2005).
#'  Environmental Kuznets curve: threatened species and spatial
#'  effects. Ecological Economics, 55(3), 395-407.
#'  @format A data.frame with 113 rows and 11 columns
##
