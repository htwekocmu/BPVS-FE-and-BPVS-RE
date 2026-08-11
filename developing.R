#' Growth determinants in 25 developing countries, 2000-2023
#'
#' Balanced panel of 25 developing countries over 24 years (2000-2023), giving
#' 600 observations. Used in the empirical application of the BPVS paper. All
#' variables are annual, country-level observations.
#'
#' @format A data frame with 600 rows and 13 columns:
#' \describe{
#'   \item{year}{Calendar year (2000-2023).}
#'   \item{id}{Country identifier, an integer 1-25.}
#'   \item{Growth}{Real GDP per capita growth rate (percent).}
#'   \item{IE}{Investment share (percent of GDP).}
#'   \item{IH}{Human capital proxy (percent).}
#'   \item{Debt}{External debt ratio (percent of GDP).}
#'   \item{FDI}{Foreign direct investment net inflows (percent of GDP).}
#'   \item{INF}{Inflation (percent).}
#'   \item{POP}{Population growth (percent).}
#'   \item{TOP}{Trade openness, (exports + imports)/GDP (percent).}
#'   \item{UR}{Unemployment rate (percent).}
#'   \item{EXR}{Nominal exchange rate index.}
#'   \item{GR}{Government revenue (percent of GDP).}
#' }
#'
#' @source Compiled from World Bank World Development Indicators.
#' @examples
#' data(developing)
#' dim(developing)
#' summary(developing$Growth)
"developing"
