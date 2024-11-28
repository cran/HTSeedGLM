

#' @title Distribution of base seed water potential

#' @description This function provides the estimates of stress and uniformity parameters along with respective variances using generalised linear model fitted to observed germination percentage of seed. The model can be fitted under logit, probit and cloglog transformations.


#' @param model Fitted model

#' @return
#' \itemize{
#'   \item Degrees of freedom
#'   \item p_value: For testing significance of water potential
#'   \item stress: Location parameter of the base seed water potential
#'   \item uniformity: Scale parameter of the base seed water potential
#'   \item var_stress: Variance of estimator of the location parameter
#'   \item var_uniformity: Variance of estimator of the scale parameter
#' }
#' @export
#' @import stats
#' @usage BaseWPDist(model)
#' @examples
#' X <- c(0,-0.3,-0.6,-0.9) # Various water potentials
#' y <- c(44,10,10,4) # Number of germinated seeds
#' n <- c(100,100,100,100) # Total number of viable seeds
#' n_y <- n-y
#' sg.mat <- cbind(y,n_y)
#' res.glm1 <- glm(sg.mat~ X,family=binomial(link=logit)) # Using logit transformation
#' my.bdl<- BaseWPDist(res.glm1)

#' res.glm2 <- glm(sg.mat~ X,family=binomial(link=probit)) # Using probit transformation
#' my.bdp<- BaseWPDist(res.glm2)

#' res.glm3 <- glm(sg.mat~ X,family=binomial(link=cloglog))# Using cloglog transformation
#' my.bdcl<- BaseWPDist(res.glm3)

#' @references
#' \itemize{
#' \item Bradford, K. J. (2002). Applications of Hydrothermal Time to Quantifying and Modeling Seed Germination and Dormancy. Weed Science, 50(2), 248–260. http://www.jstor.org/stable/4046371

#' \item Kebreab, E., & Murdoch, A. J. (1999). Modelling the effects of water stress and temperature on germination rate of Orobanche aegyptiaca seeds. Journal of Experimental Botany, 50(334), 655-664. doi:10.1093/jxb/50.334.655

#' \item Dobson, A. J., & Barnett, A. G. (2018). An introduction to generalized linear models. Chapman and Hall/CRC.
#' }

BaseWPDist <- function(model) {

  # F-Statistic and p-value#

  F_cal <- ((model$null.deviance - model$deviance)/(model$df.null - model$df.residual))/(model$deviance/model$df.residual)
  p_value <- pf(F_cal, (model$df.null - model$df.residual), model$df.residual, lower.tail = FALSE)

  # Stress(A) and Uniformity(B) from coefficients

  a <- unname(model$coefficients[1])
  b <- unname(model$coefficients[2])
  stress <- -(a/b)
  uniformity <- 1/b

  # Var(stress) and Var(uniformity)

  cov_mat <- vcov(model)
  var_stress <- ((-1/b)^2)*cov_mat[1,1] + ((a/(b^2))^2)*cov_mat[2,2] + 2*(-1/b)*(a/(b^2))*cov_mat[1,2]
  var_uniformity <- ((-1/(b^2))^2)*cov_mat[2,2]

  rm(F_cal, a, b, cov_mat) # Remove intermediaries

  # Output as named list #
  message(paste0("Degrees of freedom (", (model$df.null - model$df.residual), ",", model$df.residual, ")"))
  return(list(p_value = p_value, stress = stress, uniformity = uniformity, var_stress = var_stress, var_uniformity = var_uniformity))
}










#' @title F-test between two fitted models

#' @description This function considers two fitted models as inputs. Considering the first model as full model, it performs testing equality of uniformity parameters representing the model under null hypothesis and provides the p-value and degrees of freedom of the test statistic.


#' @param model1 First fitted model
#' @param model2 Second fitted model

#' @return
#' \itemize{
#'   \item Degrees of freedom and p-value
#' }
#' @export
#'
#' @usage FStat(model1, model2)
#' @examples
#' data1 <- data.frame(cbind(sg = c(rep(1, 95), rep(0, 5), rep(1, 87), rep(0, 13),
#'rep(1, 80), rep(0, 20), rep(1, 59), rep(0, 41),
#'rep(1, 50), rep(0, 50), rep(1, 79), rep(0, 21),
#'rep(1, 69), rep(0, 31), rep(1, 72), rep(0, 28),
#'rep(1, 44), rep(0, 56), rep(1, 14), rep(0, 86)),
#'v1 = c(rep(1, 500), rep(0, 500)),
#'v2 = c(rep(0, 500), rep(1, 500)),
#'wp1 = c(rep(0, 100), rep(-0.3, 100), rep(-0.6, 100),
#'        rep(-0.9, 100), rep(-1.2, 100), rep(0, 500)),
#'wp2 = c(rep(0, 600), rep(-0.3, 100), rep(-0.6, 100),
#'        rep(-0.9, 100), rep(-1.2, 100))))

#' data2 <- data.frame(cbind(sg = c(rep(1, 95), rep(0, 5), rep(1, 87), rep(0, 13),
#'rep(1, 80), rep(0, 20), rep(1, 59), rep(0, 41),
#'rep(1, 50), rep(0, 50), rep(1, 79), rep(0, 21),
#'rep(1, 69), rep(0, 31), rep(1, 72), rep(0, 28),
#'rep(1, 44), rep(0, 56), rep(1, 14), rep(0, 86)),
#'v1 = c(rep(1, 500), rep(0, 500)),
#'v2 = c(rep(0, 500), rep(1, 500)),
#'wp = c(rep(0, 100), rep(-0.3, 100), rep(-0.6, 100),
#'       rep(-0.9, 100), rep(-1.2, 100), rep(0, 100),
#'       rep(-0.3, 100), rep(-0.6, 100), rep(-0.9, 100),
#'       rep(-1.2, 100))))
#' myprobit1 <- glm(sg ~ v1 + v2 + wp1 + wp2 - 1, data = data1, family = binomial(link = probit))
#' myprobit2 <- glm(sg ~ v1 + v2 + wp - 1, data = data2,
#'family = binomial(link = probit))

#' my.f<- FStat(myprobit1, myprobit2)


#' @references
#' \itemize{
#' \item Bradford, K. J. (2002). Applications of Hydrothermal Time to Quantifying and Modeling Seed Germination and Dormancy. Weed Science, 50(2), 248–260. http://www.jstor.org/stable/4046371

#' \item Kebreab, E., & Murdoch, A. J. (1999). Modelling the effects of water stress and temperature on germination rate of Orobanche aegyptiaca seeds. Journal of Experimental Botany, 50(334), 655-664. doi:10.1093/jxb/50.334.655

#' \item Dobson, A. J., & Barnett, A. G. (2018). An introduction to generalized linear models. Chapman and Hall/CRC.
#' }

FStat <- function(model1, model2) {
  message("The first model is considered as full model")
  # F-Statistic and p-value #
  F_cal <- ((model2$deviance - model1$deviance)/(model2$df.residual - model1$df.residual))/(model1$deviance/model1$df.residual)
  p_value <- pf(F_cal, (model2$df.residual - model1$df.residual), model1$df.residual, lower.tail = FALSE)
  rm(F_cal)
  message(paste0("The p-value of the F test at (", (model2$df.residual - model1$df.residual), ",", model1$df.residual,") degrees of freedom is ", p_value))
}
