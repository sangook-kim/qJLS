#' qJLS (quantile regression based joint location and scale) Test
#'
#' This function performs a qJLS (quantile regression based joint location and scale) test (Kim and Strug 2023) using "quantreg" package (Koenker et al.).
#' @param model_G a formula object to adjust the genotype for covariates, with the genotype (G) on the left and a ~ operator separating the covariates of interest on the right, separated by + operators.
#' @param b_loc a vector of rankscores for the location.  This is required.
#' @param b_scale a vector of rankscores for the scale  This is required.
#' @param data a data frame containing the variables named in model and correlation arguments.  This is required.
#' @export
#' @author Sangook Kim
#' @import quantreg
#' @details No missing data are allowed - function will return an "error". Outcome (phenotype) must be quantitative and genotype may be discrete (categorical) or continuous.
#' @return pval the p-value for the qJLS test
#' @references Kim, S. and Strug, L.J. (2023). A robust association test leveraging unknown genetic interactions: Application to cystic fibrosis lung disease (Submitted)
#' @examples
#' #################################################################################
#' ## Example simulating data from Model 3 (Kim and Strug 2023 Plos Genetics)
#' #################################################################################
#'
#' #### Simulation parameters
#' n <- 2000  # sample size
#' maf <- 0.3 # minor allele frequency
#' pE <- 0.3  # exposure (E) probability
#' alpha <- 0.1 # correlation parameter between genotype (G) and covariate (Z)
#' bE <- 0.3  # effect of E on response (Y)
#' bZ <- 0.3  # effect of Z on Y
#' bG <- 0.01 # effect of G on Y
#' bGE <- 0.1  # GxE interaction effect
#'
#' p <- c((1-maf)^2, 2*(1-maf)*maf, maf*maf)
#' pcum <- cumsum(p)
#'
#' #### Generates E, G and Z
#' E <- sample(c(0,1), n, replace=TRUE, prob=c(1-pE, pE))
#' X <- mvrnorm(n, mu=c(0,0), Sigma=matrix(c(1,alpha,alpha,1), nrow=2))
#' G <- cut(X[,1], breaks=c(-Inf, qnorm(pcum[1:2]), +Inf), labels=0:2)
#' G <- as.numeric(levels(G))[G]
#' Z <- X[,2]
#'
#' #### Generates error based on normal distribution and Y
#' err <- rnorm(n)
#' mu <- bG*G + sign(bGE)*bE*E + bGE*G*E + bZ*Z
#' Y <- mu + err
#'
#' df <- data.frame(Y, G, E, Z)
#'
#' #### Generates rank
#' mod0 <- rq(Y ~ Z, tau=-1, data=df)
#' b_loc <- ranks(mod0, score='normal')$ranks
#' b_scale <- ranks(mod0, score='normalscale')$ranks
#'
#' qJLS_test(G ~ Z, b_loc, b_scale, df)
#'

qJLS_test <- function(model_G, b_loc, b_scale, data)
{
  B <- cbind(b_loc, b_scale)
  Gstar <- lm(model_G, data=data)$residuals
  A <- matrix(c(1,0,0,2), nrow=2)
  S <- t(Gstar) %*% B
  M <- A * (sum(Gstar^2))
  T_rs <- S^2/diag(M)
  return(pchisq(sum(T_rs), length(S), lower.tail=FALSE))
}






