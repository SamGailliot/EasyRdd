#' @title rd_plot
#'
#' @description rd_plot calls the powerful function rdplot from rdrobust but makes only the most important options available to the user to make the function easier to use.
#'
#' @param Y Nx1 vector containing dependent variable
#' @param Y Nx1 vector containing running variable (forcing variable)
#' @param cutpoint Location of discontinuity.
#' @param bin_method specifies the procedure to select the number of bins. Options are:
#'
#' \code{es}: IMSE-optimal evenly-spaced method using spacings estimators.
#'
#' \code{espr}: IMSE-optimal evenly-spaced method using polynomial regression.
#'
#' \code{esmv}: mimicking variance evenly-spaced method using spacings estimators. Default
#'
#' \code{esmvpr}: mimiching variance evenly-spaced method using  polynomial regression
#'
#' \code{qs}: IMSE-optimal quantile-spaced method using spacings estimators.
#'
#' \code{qspr}: IMSE-optimal quantile-spaced method using polynomial regression.
#'
#' \code{qsmv}: mimicking variance quantile-spaced method using spacings estimators.
#'
#' \code{qsmvpr}: mimicking variance quantile-spaced method using polynomial regression.
#'
#' @param kernel_method specifies the vernel function used to construct local-polynomial estimators. Options are: \code{triangular}, \code{epanechnikov}, and \code{uniform}
#' @param ci optional graphical option to display confidence intervals of selected level (ex. 95) for each bin.
#' @param title optional title for the plot
#' @param x.lab optional title for the x-axis
#' @param y.lab optional title for the y-axis
#' @param x.lim optional setting for the range of the x-axis
#' @param y.lim optional setting for the range of the y-axis
#'
#' @return \code{rdplot} a standard \code{ggplot} object that can be used for further customization.
#'
#' @author Samuel Gailliot, Texas A (and) M University, College Station, TX. samuel.gailliot@stat.tamu.edu
#'
#' @references Calonico, S., M. D. Cattaneo, M. H. Farrell, and R. Titiunik. 2017. rdrobust: Software for Regression Discontinuity Designs. Stata Journal 17(2): 372-404.
#'
#' Calonico, S., M. D. Cattaneo, and R. Titiunik. 2014. Robust Data-Driven Inference in the Regression-Discontinuity Design. Stata Journal 14(4): 909-946.
#'
#' Calonico, S., M. D. Cattaneo, and R. Titiunik. 2015a. Optimal Data-Driven Regression Discontinuity Plots. Journal of the American Statistical Association 110(512): 1753-1769.
#'
#' Calonico, S., M. D. Cattaneo, and R. Titiunik. 2015b. rdrobust: An R Package for Robust Nonparametric Inference in Regression-Discontinuity Designs. R Journal 7(1): 38-51.
#'
#' Cattaneo, M. D., B. Frandsen, and R. Titiunik. 2015. Randomization Inference in the Regression Discontinuity Design: An Application to the Study of Party Advantages in the U.S. Senate. Journal of Causal Inference 3(1): 1-24.
#'
#' @export

rd_plot <- function(Y, X, cutpoint = 0, bin_method = "esmv", kernel_method = "uni", ci = NULL,
                    title = NULL, x.lab = NULL, y.lab = NULL, x.lim = NULL, y.lim = NULL){

rdrobust::rdplot(y = Y, x = X, c = cutpoint, nbins = NULL, binselect = bin_method,
                 scale = NULL, kernel = kernel_method, weights = NULL, h = NULL,
                 covs = NULL, covs_eval = 0, covs_drop = TRUE, support = NULL,
                 subset = NULL, hide = FALSE, ci = ci, shade = FALSE, title = title,
                 x.label = x.lab, y.label = y.lab, x.lim = x.lim, y.lim = y.lim,
                 col.dots = NULL, col.lines = NULL)
}
