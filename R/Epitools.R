#' inter.coeff
#'
#' This function allows you to caculate coefficient and confidence interval of each stratum of interaction term.
#' @param mod a model include interaction fiited by coxme
#' @param gp a numeric indicate the level of interacted factor
#' @examples
#' inter.coeff(model1,2)

inter.coeff <- function (mod,gp){
  require(coxme)
  g <- vector()
  coff <- vector()
  se <- vector()
  p <- vector()


  for (i in 1:gp) {
    g[i] <- names(mod$coefficients)[length(mod$coefficients)-gp+i]
  }


  beta <- mod$coefficients
  rmain <- as.numeric(which(beta==mod$coefficients[g[1]]))
  cmain <- as.numeric(which(beta==mod$coefficients[g[1]]))
  coff[1] <- beta[g[1]]
  se[1] <- sqrt(vcov(mod)[rmain,cmain])
  p[1] <- NA


  for (i in 2:gp) {
    rinter <- as.numeric(which(beta==mod$coefficients[g[i]]))
    cinter <- as.numeric(which(beta==mod$coefficients[g[i]]))
    coff[i] <- beta[rinter]+beta[rmain]
    se[i] <- sqrt(vcov(mod)[rmain,cmain]+vcov(mod)[rinter,cinter]+2*vcov(mod)[rmain,cinter])
    p[i] <- signif(1 - pchisq((beta[rinter]/sqrt(vcov(mod)[rinter,cinter]))^2, 1), 2)
  }


  table <- as.data.frame(cbind(g,coff,se,p))
  names(table)[1] <- "group"
  for (i in 2:4) {
    table[,i] <- as.character(table[,i])
    table[,i] <- as.numeric(table[,i])
  }
  or <- exp(table$coff*10)
  low <- exp(table$coff*10-1.96*table$se*10)
  up <- exp(table$coff*10+1.96*table$se*10)
  table$or <- paste(round(or,digits = 2)," (",round(low,digits = 2),", ",round(up,digits = 2),")",sep="")
  return(table)
}




#' inter.reri
#'
#' This function allows you to caculate RERI, APAB and S with confidence interval of interaction term.
#' @param model a model include interaction fiited by coxme
#' @param coeff a vector indicate the location of joint risk AB, seperate risk A and seperate risk B
#' @param type a string indicate the needed statistics ("RERI", "APAB" or "S")
#' @param conf.level a numeric indicate the level of confidence interval
#' @examples
#' inter.reri(model1,c(2,3,4),"RERI",0.95)

inter.reri <- function (model, coeff, type = c("RERI", "APAB", "S"), conf.level = 0.95)
{
  require(coxme)
  N. <- 1 - ((1 - conf.level)/2)
  z <- qnorm(N., mean = 0, sd = 1)
  if (type == "RERI") {
    if (class(model)[1] != "glm" & class(model)[2] != "lm" &
        class(model)[1] != "clogit" & class(model)[1] !=
        "coxph" & class(model)[1] != "coxme")
      stop("Error: model must be either a glm or coxph object")
    if (class(model)[1] == "glm" & class(model)[2] == "lm") {
      theta1 <- as.numeric(model$coefficients[coeff[1]])
      theta2 <- as.numeric(model$coefficients[coeff[2]])
      theta3 <- as.numeric(model$coefficients[coeff[3]])
    }
    if (class(model)[1] == "clogit" | class(model)[1] == "coxme" | class(model)[1] ==
        "coxph") {
      theta1 <- as.numeric(model$coefficients[coeff[1]])
      theta2 <- as.numeric(model$coefficients[coeff[2]])
      theta3 <- as.numeric(model$coefficients[coeff[3]])
    }
    cov.mat <- vcov(model)
    h1 <- -exp(theta1)
    h2 <- -exp(theta2)
    h3 <- exp(theta3)
    reri.var <- (h1^2 * (cov.mat[coeff[1], coeff[1]])) +
      (h2^2 * (cov.mat[coeff[2], coeff[2]])) + (h3^2 *
                                                  (cov.mat[coeff[3], coeff[3]])) + (2 * h1 * h2 * cov.mat[coeff[1],
                                                                                                          coeff[2]]) + (2 * h1 * h3 * cov.mat[coeff[1], coeff[3]]) +
      (2 * h2 * h3 * cov.mat[coeff[2], coeff[3]])
    reri.se <- sqrt(reri.var)
    reri.p <- exp(theta3) - exp(theta1) - exp(theta2) + 1
    reri.l <- reri.p - (z * reri.se)
    reri.u <- reri.p + (z * reri.se)
    rval <- data.frame(reri.p, reri.l, reri.u)
    names(rval) <- c("est", "lower", "upper")
  }
  if (type == "APAB") {
    if (class(model)[1] != "glm" & class(model)[2] != "lm" &
        class(model)[1] != "clogit" & class(model)[1] != "coxme" & class(model)[1] !=
        "coxph")
      stop("Error: model must be either a glm or coxph object")
    if (class(model)[1] == "glm" & class(model)[2] == "lm") {
      theta1 <- as.numeric(model$coefficients[coeff[1]])
      theta2 <- as.numeric(model$coefficients[coeff[2]])
      theta3 <- as.numeric(model$coefficients[coeff[3]])
    }
    if (class(model)[1] == "clogit" | class(model)[1] == "coxme" | class(model)[1] ==
        "coxph") {
      theta1 <- as.numeric(model$coefficients[coeff[1]])
      theta2 <- as.numeric(model$coefficients[coeff[2]])
      theta3 <- as.numeric(model$coefficients[coeff[3]])
    }
    cov.mat <- vcov(model)
    h1 <- -exp(theta1 - theta3)
    h2 <- -exp(theta2 - theta3)
    h3 <- (exp(theta1) + exp(theta2) - 1)/exp(theta3)
    apab.var <- (h1^2 * (cov.mat[coeff[1], coeff[1]])) +
      (h2^2 * (cov.mat[coeff[2], coeff[2]])) + (h3^2 *
                                                  (cov.mat[coeff[3], coeff[3]])) + (2 * h1 * h2 * cov.mat[coeff[1],
                                                                                                          coeff[2]]) + (2 * h1 * h3 * cov.mat[coeff[1], coeff[3]]) +
      (2 * h2 * h3 * cov.mat[coeff[2], coeff[3]])
    apab.se <- sqrt(apab.var)
    apab.p <- (exp(theta3) - exp(theta1) - exp(theta2) +
                 1)/exp(theta3)
    apab.l <- apab.p - (z * apab.se)
    apab.u <- apab.p + (z * apab.se)
    rval <- data.frame(apab.p, apab.l, apab.u)
    names(rval) <- c("est", "lower", "upper")
  }
  return(rval)
}