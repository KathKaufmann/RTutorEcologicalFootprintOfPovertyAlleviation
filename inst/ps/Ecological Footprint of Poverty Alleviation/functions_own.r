#Load those results from exercises 1.1 and 1.2 which are relevant later
solve_1_1 = function(){
  data_deforest = readRDS("data_deforestation.rds")
  hist_pov_ind = ggplot(data_deforest, aes(x=pov_ind_95)) + 
    geom_histogram(binwidth=0.1, center=-1.17, aes(y=..density..), col="blue") +
    geom_vline(xintercept = -1.22, col="red")
  data_deforest = mutate(data_deforest, 
                         aux_bins_povind = round(pov_ind_95 - 0.03, digits = 1))
  dat = summarize(group_by(data_deforest, aux_bins_povind), 
                  prop_treated = mean(treat))
  p_proptreated = ggplot(data=dat, aes(x=aux_bins_povind +0.03, y=prop_treated)) +
    geom_point(col="blue") +
    geom_vline(xintercept = -1.22, col="red")+
    xlab("bins marginality index")
  return(list(p_proptreated=p_proptreated, hist_pov_ind=hist_pov_ind))
}

solve_1_2 = function(){
  data_deforest = readRDS("data_deforestation.rds")
  data_deforest = mutate(data_deforest, 
                         aux_bins_povind = round(pov_ind_95 - 0.03, digits = 1))
  dat_pct = summarize(group_by(data_deforest, aux_bins_povind), 
                      mean_pctdefor = mean(pct_deforested))
  p_pctdefor = ggplot(data=dat_pct, aes(x=aux_bins_povind+0.03, y=mean_pctdefor)) +
    geom_point(col="darkgreen") + 
    geom_vline(xintercept = -1.22, col="red")+
    xlab("bins marginality index")
  return(p_pctdefor)
}

#This function estimates the marginal effects of each explanatory variable
#on P(y>0|x) and E[y|x,y>0] for the Tobit model.

#input: - Tobit model estimated by command "tobit" from package 'AER'
#       - corresponding data frame
margEff_tobit_dM_AER = function(tobit,data){
  #load package 'car' for later use of delta method
  library(car)
  
  #extract estimates from estimated Tobit model for coefficients and
  #for the logarithmized variance of the error term epsilon
  tobit_coef = summary(tobit)$coef[,1]
  
  #define number of estimates
  n = length(tobit_coef)
  
  #extract estimates of the coefficients {beta_j} from Tobit
  tobit_beta = tobit_coef[-n]
  
  #concerning all estimates (tobit_coef):
  #rename "intercept" "beta0" and "log(Sigma)" "lnsigma"
  #save renamed coefficients in renamedcoefs
  renamedcoefs = tobit_coef
  names(renamedcoefs)[c(1,n)] = c("beta0", "lnsigma")
  
  #extract explanatory variables from data
  cols = names(tobit_beta)[-1]
  X = data[,cols]
  #calculate the mean of all explanatory variables and name them appropriately
  if ( length(cols)==1){
    X_mean=mean(X)
  } else{
    X_mean=colMeans(X, na.rm=TRUE)
  }
  names(X_mean) = paste(cols, "mean", sep="")
  
  #determine binary variables and continuous variables
  firstfound = FALSE
  for (i in 1:(n-2)) {
    if (length(cols)==1){
      num_contained = unique(X)
    } else{
      num_contained = unique(X[,i])
    }
    if (length(num_contained)==2){
      if ((num_contained[1]==0 & num_contained[2]==1) | (num_contained[1]==1 & num_contained[2]==0)){
        if (firstfound==FALSE){
          firstfound = TRUE
          bin = c(i)
        } else{
          bin = c(bin, i)
        }
      }
    }
    if (length(num_contained)==1){
      if (num_contained==0 | num_contained==1){
        if (firstfound==FALSE){
          firstfound = TRUE
          bin = c(i)
        } else{
          bin = c(bin, i)
        }
      }
    }
  }
  
  conti = seq(1:(n-2))[-bin]
  
  #for each type of marginal effect:
  #- define the formula as a string
  #- calculate estimate of marginal effect and standard error with delta method
  #- add t value, corresponding p value and significance level
  #  to resulting data frame
  #- name the resulting data frame
  
  #marginal effects for continuous variables
  firstconti = FALSE
  for (i in conti){
    #marginal effect on P(y>0|x)
    formula_P_conti_i = paste("(dnorm((beta0 + ", paste(cols, names(X_mean), sep="*", collapse=" + "),
                              ")/exp(lnsigma))/exp(lnsigma)) * ", cols[i], sep="")
    margEff_P_conti_i = deltaMethod(object = renamedcoefs, vcov. = vcov(tobit),
                                    g = formula_P_conti_i, constants = X_mean)
    margEff_P_conti_i = mutate(margEff_P_conti_i, tvalue = Estimate/SE, 
                               pvalue = 2*pt(-abs(tvalue), nrow(data) - (n-1)),
                               significance = as.factor(ifelse(pvalue<0.1, 
                                                               ifelse(pvalue<0.05,
                                                                      ifelse(pvalue<0.01, "***" , "**") , "*"), "")))
    rownames(margEff_P_conti_i)[1] = cols[i]
    colnames(margEff_P_conti_i)[1] = "Estimate P(y>0)"
    
    #marginal effect on E[y|x,y>0]
    formula_E_conti_i = paste("(1-(dnorm((beta0 + ", 
                              paste(cols, names(X_mean), sep="*", collapse= " + ") ,
                              ")/exp(lnsigma))/pnorm((beta0 + ", paste(cols, names(X_mean), sep="*", collapse= " + ") ,
                              ")/exp(lnsigma))) * (  ((beta0 + ", paste(cols, names(X_mean), sep="*", collapse = " + "),
                              ")/exp(lnsigma)) + (dnorm((beta0 + ", paste(cols, names(X_mean), sep="*", collapse= " + ") ,
                              ")/exp(lnsigma))/pnorm((beta0 + ", paste(cols, names(X_mean), sep="*", collapse= " + ") ,
                              ")/exp(lnsigma)))  )) * ", cols[i], sep="")
    margEff_E_conti_i = deltaMethod(object = renamedcoefs, vcov. = vcov(tobit),
                                    g = formula_E_conti_i, constants = X_mean)
    margEff_E_conti_i = mutate(margEff_E_conti_i, tvalue = Estimate/SE, 
                               pvalue = 2*pt(-abs(tvalue), nrow(data) - (n-1)),
                               significance = as.factor(ifelse(pvalue<0.1, 
                                                               ifelse(pvalue<0.05,
                                                                      ifelse(pvalue<0.01, "***" , "**") , "*"), "")))
    rownames(margEff_E_conti_i)[1] = cols[i]
    colnames(margEff_E_conti_i)[1] = "Estimate E[y|y>0]"
    
    #Save results in data frames for each kind of marginal effect
    #of all continuous variables
    if (firstconti==FALSE){
      firstconti = TRUE
      margEff_P_conti = margEff_P_conti_i
      margEff_E_conti = margEff_E_conti_i
    } else{
      margEff_P_conti = rbind(margEff_P_conti, margEff_P_conti_i)
      margEff_E_conti = rbind(margEff_E_conti, margEff_E_conti_i)
    }
  }
  
  #marginal effects for binary variables
  firstbin = FALSE
  for (i in bin){
    #adapt the mean of the explanatory variables for the calculations below
    #and save them in X_mean_bin
    X_mean_bin = X_mean
    X_mean_bin[i] = 1
    
    #marginal effect on P(y>0|x)
    formula_P_bin_i = paste("pnorm((beta0 + ", paste(cols, names(X_mean_bin), sep="*", collapse= " + ") ,
                            ")/exp(lnsigma))", " - ", 
                            "pnorm((beta0 + ", paste(cols, names(X_mean_bin), sep="*", collapse= " + "), 
                            " - ", cols[i] , ")/exp(lnsigma))", sep="")
    margEff_P_bin_i = deltaMethod(object = renamedcoefs, vcov. = vcov(tobit),
                                  g = formula_P_bin_i, constants = X_mean_bin)
    margEff_P_bin_i = mutate(margEff_P_bin_i, tvalue = Estimate/SE, 
                             pvalue = 2*pt(-abs(tvalue), nrow(data) - (n-1)),
                             significance = as.factor(ifelse(pvalue<0.1, 
                                                             ifelse(pvalue<0.05,
                                                                    ifelse(pvalue<0.01, "***" , "**") , "*"), "")))
    rownames(margEff_P_bin_i)[1] = cols[i]
    colnames(margEff_P_bin_i)[1] = "Estimate P(y>0)"
    
    #marginal effect on E[y|x,y>0]
    formula_E_bin_i = paste(cols[i], " + exp(lnsigma) * (dnorm((beta0 + ", 
                            paste(cols, names(X_mean_bin), sep="*", collapse= " + ") ,
                            ")/exp(lnsigma))/pnorm((beta0 + ",
                            paste(cols, names(X_mean_bin), sep="*", collapse= " + ") ,
                            ")/exp(lnsigma)))", " - ",
                            "exp(lnsigma) * (dnorm((beta0 + ",
                            paste(cols, names(X_mean_bin), sep="*", collapse= " + ") ,
                            " - ", cols[i], ")/exp(lnsigma))/pnorm((beta0 + ",
                            paste(cols, names(X_mean_bin), sep="*", collapse= " + ") ,
                            " - ", cols[i], ")/exp(lnsigma)))", sep="")
    margEff_E_bin_i = deltaMethod(object = renamedcoefs, vcov. = vcov(tobit),
                                  g = formula_E_bin_i, constants = X_mean_bin)
    margEff_E_bin_i = mutate(margEff_E_bin_i, tvalue = Estimate/SE, 
                             pvalue = 2*pt(-abs(tvalue), nrow(data) - (n-1)),
                             significance = as.factor(ifelse(pvalue<0.1, 
                                                             ifelse(pvalue<0.05,
                                                                    ifelse(pvalue<0.01, "***" , "**") , "*"), "")))
    rownames(margEff_E_bin_i)[1] = cols[i]
    colnames(margEff_E_bin_i)[1] = "Estimate E[y|y>0]"
    
    #Save results in data frames for each kind of marginal effect
    #of all binary variables
    if (firstbin==FALSE){
      firstbin = TRUE
      margEff_P_bin = margEff_P_bin_i
      margEff_E_bin = margEff_E_bin_i
    } else{
      margEff_P_bin = rbind(margEff_P_bin, margEff_P_bin_i)
      margEff_E_bin = rbind(margEff_E_bin, margEff_E_bin_i)
    }
    
  }
  
  #Create the data frame or warning to return
  if (firstbin ==TRUE) {
    if (firstconti == TRUE){
      return (list(margEff_P = rbind(margEff_P_bin, margEff_P_conti),
                   margEff_E = rbind(margEff_E_bin, margEff_E_conti)))
    } else {
      return (list(margEff_P = margEff_P_bin,
                   margEff_E = margEff_E_bin))
    }
  } else{
    if (firstconti == TRUE){
      return (list(margEff_P = margEff_P_conti,
                   margEff_E = margEff_E_conti))
    } else {
      return ("Error: No marginal effects available.")
    }
  }
}
ksmooth_wrap = function(y, x, kernel){
  res = as.matrix(as.data.frame(ksmooth(x,y, kernel)))
  return(res)
}

library(AER)
library(dplyr)
library(dplyrExtras)
#function to calculate Newey two-step estimates for ivtobit 
#and marginal effects on P(y>0) and on E[y|y>0]

#implementation of Newey two-step estimator based on:
#- Newey, EFFICIENT ESTIMATION OF LIMITED DEPENDENT VARIABLE MODELS WITH ENDOGENOUS EXPLANATORY VARIABLES,
#  Journal of Econometrics 36 (1987), p.231-250
#- http://www.stata.com/manuals13/rivtobit.pdf
#- http://www.stata.com/manuals13/rivprobit.pdf p.11-12

#input (all as strings): ystr - dependent variable, Ystr - endogenous explanatory variable (only one!),
#X1str - exogenous explanatory variables, X2str - excluded instruments
#General remark: Please note that I indicate the maximization condition for the estimates calculated here,
#                but I don't use the function censReg ('censReg' package) calculating maximum likelihood
#                estimates for reasons of computational time but the command tobit ('AER' package) which
#                should deliver the same numerical results
Newey_SB_ivtobit = function(data, ystr, Ystr, X1str, X2str){
  #Extract the observations of the variables of interest from the data
  X1 = subset(data,select=c(X1str))
  X2 = subset(data,select=c(X2str))
  y = subset(data, select=c(ystr))
  Y = subset(data, select=c(Ystr))
  
  #(I) OLS regression: Y_t = X1_t pi1 + X2_t pi2 + V_t = X_t pi + V_t
  x1_formula = paste(X1str, collapse = " + ")
  x2_formula = paste(X2str, collapse = " + ")
  formula_1 = paste(Ystr," ~ ", x1_formula, " + ",x2_formula, sep="")
  first = lm(formula_1, data=data)
  #save estimate of pi in pi_l_hat which equals theta_1_l_hat
  pi_l_hat=coef(first)
  #Use pi_l_hat to calculate D_hat
  D_hat = cbind(pi_l_hat, rbind(diag(x=1, nrow=ncol(X1)+1, ncol=ncol(X1)+1), 
                                matrix(0,nrow=ncol(X2) , ncol(X1)+1)))
  
  #(II) two-step estimator: 
  #max{alpha, lambda, sigma^2, psi} sum{t}l(y_t, X_t alpha + V_t_hat lambda, sigma^2, psi)/n
  #--> take residuals from OLS (I) as V_t_hat
  #    tobit regression: y_t = max(0, X_t alpha + V_t_hat lambda + e_t) with e_t|Y ~ N(0,sigma^2)
  formula_2 = as.formula(paste0(ystr," ~ ",x1_formula," + ",x2_formula," + residuals(first)"))
  second = tobit(formula_2, data=data)
  n_2 = length(coef(second))
  #save estimate of alpha in alpha_l_hat, estimate of lambda in lambda_l_hat
  alpha_l_hat = coef(second)[1:(n_2-1)]
  lambda_l_hat = coef(second)[n_2]
  #calculate estimate of (J^(-1))alpha_alpha
  J_alpha_alpha = vcov(second)[(1:(n_2-1)), (1:(n_2-1))]
  
  #(III) 2SIV estimator:
  #max{delta, rho, sigma^2, psi} sum{t}l(y_t, Z_t delta + V_t_hat rho, sigma^2, psi)
  #with Z_t = [Y_t, X_1_t]
  #--> take residuals from OLS (I) as V_t_hat
  #    tobit regression: y_t = max(0, Z_t delta + V_t_hat rho + e_t)
  #                          = max(0, Y_t beta + X_1_t gamma0 + V_t_hat rho + e_t)
  formula_3 = as.formula(paste0(ystr," ~ ", Ystr, " + ", x1_formula, " + residuals(first)"))
  third = tobit(formula_3, data=data)
  #save estimate of beta in beta_l_hat and name the elements
  beta_l_hat = coef(third)[2]
  names(beta_l_hat)=paste0(Ystr, "_hat")
  
  #(IV) OLS regression of Y_t (lambda_hat - beta_hat) on X_t
  #--> take lambda_l_hat from (II) as lambda_hat and beta_l_hat from (III) as beta_hat
  Y_tilde = as.numeric(as.matrix(Y*(lambda_l_hat - beta_l_hat)))
  formula_4 = as.formula(paste0("Y_tilde", " ~ ", x1_formula, " + ", x2_formula))
  fourth = lm(formula_4, data=data)
  #save estimate of covariance matrix in vcov_tilde
  #this equals Nu__hat^2 (X'X/n)^(-1)
  vcov_tilde = vcov(fourth)
  
  #Calculate Omega_l_hat
  omega_l_hat = J_alpha_alpha + vcov_tilde
  
  #Calculate Newey estimator delta_Al_hat and the estimate of its asymptotic covariance matrix
  delta_Al_hat = solve(t(D_hat)%*%solve(omega_l_hat)%*%D_hat)%*%t(D_hat)%*%solve(omega_l_hat)%*%alpha_l_hat
  vcov_hat_delta = solve(t(D_hat)%*%solve(omega_l_hat)%*%D_hat)
  #add according variable names and column names to vcov_hat_delta
  row.names(vcov_hat_delta)=c(Ystr, "Intercept", colnames(X1))
  colnames(vcov_hat_delta)=c(Ystr, "Intercept", colnames(X1))
  #calculate standard errors of delta_Al_hat
  se_delta_hat = sqrt(diag(vcov_hat_delta))
  
  #Output preparation
  df = data.frame(Estimate = as.vector(delta_Al_hat), SE = as.vector(se_delta_hat))
  
  #Adding t values, p values and significance stars
  #and setting row names of resulting data frame
  df = mutate(df, tvalue = delta_Al_hat/se_delta_hat, 
              pvalue = 2*pt(-abs(tvalue), nrow(data - length(X1str)-length(Ystr) - 1)),
              significance = as.factor(ifelse(pvalue<0.1, ifelse(pvalue<0.05,ifelse(pvalue<0.01, "***" , "**") , "*"), "")))
  row.names(df)=c(Ystr, "Intercept", colnames(X1))
  ###########################################################
  #marginal effects
  row.names(delta_Al_hat)=c(Ystr, "Intercept", colnames(X1))
  names(se_delta_hat)=c(Ystr, "Intercept", colnames(X1))
  
  #estimates for tau1, tau2
  #lntau1 from third regression, tau2sqr from first OLS regression
  n_tobitSM_coef = nrow(summary(third)$coef)
  lntau1 = summary(third)$coef[n_tobitSM_coef,1]
  names(lntau1)="lntau1"
  tau2sqr = var(residuals(first))
  names(tau2sqr)="tau2sqr"
  
  #rho_hat from estimation (III)
  rho_hat_SM = coef(third)[length(coef(third))]
  names(rho_hat_SM)="rho_hat_SM"
  
  #used estimates: all (except tau2sqr) from third regression
  #in order to get a suitable covariance matrix for the estimates
  gamma0_SM_hat = coef(third)[-c(2,length(coef(third)))]
  names(gamma0_SM_hat)=c("Intercept", colnames(X1))
  beta_SM_hat = coef(third)[2]
  names(beta_SM_hat)=c(Ystr)
  X1_mean = colMeans(cbind(1,X1), na.rm=TRUE)
  names(X1_mean)[1]="Intercept"
  names(X1_mean)=paste(names(X1_mean),"mean", sep="")
  #formula input for deltaMethod
  # formula written normally
  # p_marge = pnorm((X1_mean%*%gamma0_SM_hat + beta_SM_hat)/sqrt(var.res)) - 
  #                        pnorm((X1_mean%*%gamma0_SM_hat)/sqrt(var.res))
  formula_p = paste("pnorm((",
                    paste(names(gamma0_SM_hat), names(X1_mean), sep="*", collapse= " + "), 
                    "+ ", names(beta_SM_hat), 
                    ")/sqrt(exp(2*lntau1) + rho_hat_SM^2 * tau2sqr)",
                    ") - pnorm((",
                    paste(names(gamma0_SM_hat), names(X1_mean), sep="*", collapse= " + "),
                    ")/sqrt(exp(2*lntau1) + rho_hat_SM^2 * tau2sqr)", ")")
  
  #variance/covariance input for deltaMethod
  vcov_help = vcov(third)
  vcov_1 = vcov_help[,-2]
  vcov_1 = cbind(vcov_help[,2],vcov_1)
  vcov_2 = vcov_1[-2,]
  vcov_2 = rbind(vcov_1[2,],vcov_2)
  vcov_p = vcov_2
  vcov_p=rbind(cbind(vcov_p,0),0)
  n_p = nrow(vcov_p)
  row.names(vcov_p)[(n_p-2):n_p]=c("rho_hat_SM", "lntau1", "tau2sqr")
  colnames(vcov_p)[(n_p-2):n_p]=c("rho_hat_SM", "lntau1", "tau2sqr")
  
  margeP = deltaMethod(object=c(beta_SM_hat, gamma0_SM_hat, rho_hat_SM, lntau1, tau2sqr),
                       vcov.=vcov_p, g=formula_p, constants = X1_mean)
  row.names(margeP)="marginal effect P(y>0)"
  
  #marginal effect on E[y|y>0]
  #used estimates: all from third regression in order to get a suitable covariance matrix for the estimates
  
  #variance/covariance matrix of input parameters of marginal effects calculation
  vcov_E = vcov_p[-nrow(vcov_p),]
  vcov_E = vcov_E[,-ncol(vcov_E)]
  
  #calculate marginal effect
  n_res=length(residuals(first))
  xbeta_u = X1_mean%*%gamma0_SM_hat + beta_SM_hat + rho_hat_SM *residuals(first)
  xbeta_d = X1_mean%*%gamma0_SM_hat + rho_hat_SM*residuals(first)
  tau1 = exp(lntau1)
  lambda_u = dnorm(xbeta_u/tau1)/pnorm(xbeta_u/tau1)
  lambda_d = dnorm(xbeta_d/tau1)/pnorm(xbeta_d/tau1)
  
  margeE = (1/n_res)*sum(xbeta_u + tau1*lambda_u -
                           (xbeta_d + tau1*lambda_d))
  
  #calculate gradient of margeE
  C_beta = 1-(1/n_res)*sum(lambda_u*((xbeta_u/tau1) + lambda_u))
  C_gamma = (1/n_res)*X1_mean*(sum(lambda_d*((xbeta_d/tau1)+lambda_d)) - 
                                 sum(lambda_u*((xbeta_u/tau1)+lambda_u)))
  C_rho = (1/n_res)*sum(residuals(first)*lambda_d*((xbeta_d/tau1)+lambda_d)) - 
    (1/n_res)*sum(residuals(first)*lambda_u*((xbeta_u/tau1)+lambda_u))
  C_lntau1 = (1/n_res)*tau1*sum(lambda_u-lambda_d) + 
    (1/n_res)*(1/(tau1)^2)*sum((xbeta_u*xbeta_u*lambda_u) - (xbeta_d*xbeta_d*lambda_d)) + 
    (1/n_res)*(1/tau1)*sum((xbeta_u*lambda_u*lambda_u) - (xbeta_d*lambda_d*lambda_d))
  C=as.vector(c(C_beta, C_gamma, C_rho, C_lntau1))
  
  #calculate variance/covariance matrix of margeE according to the delta method
  vcov_margE = t(C)%*%vcov_E%*%C
  se_margE = sqrt(vcov_margE)
  
  
  #prepare RESULTS for output
  df_marge = data.frame(Estimate=c(margeP[1,1],margeE), SE=c(margeP[1,2],se_margE))
  
  #Adding t values, p values, significance stars and add row names
  df_marge = mutate(df_marge, tvalue = Estimate/SE,
                    pvalue = 2*pt(-abs(tvalue), nrow(data - length(X1str)-length(Ystr) - 1)),
                    significance = as.factor(ifelse(pvalue<0.1, ifelse(pvalue<0.05, ifelse(pvalue<0.01, "***" , "**") , "*"), "")))
  rownames(df_marge)=c("marginal effect on P(y>0)", "marginal effect on E[y|y>0]")
  #row.names(df_marge)=c("marginal effect on P(y>0)", "marginal effect on E[y|y>0]")
  
  results = list("coefficients" = df, "marginalEffects" = df_marge)
  return(results)
}


#This function converts a list of data frames containing at least estimates,
#standard errors and significance levels into one single data frame
#which can be used for a clear presentation of all results.

#input: -listOfModels = list of data frames containing regression results
#        for different models; data frames must have the characteristics:
#           * 1 column name contains the term "Estimate"
#           * 1 column name contains the term "SE"
#           * 1 column name is "significance"
#           * each row contains (at least) estimates, std. errors and
#             significance levels for one variable
convert = function(listOfModels){
  #Extract the variable names contained in any of the models
  coefnames = c()
  for (i in 1:length(listOfModels)){
    coefnames = c(coefnames,row.names(listOfModels[[i]]))
  }
  coefnames = unique(coefnames)
  
  #Dummy indicating if currently the first model (first element of
  #the list) is treated
  firstModel = TRUE
  
  #Go through all data frames in the listOfModels...
  for (m in 1:length(listOfModels)){
    #Separate the column containing the significance level
    #from the rest of the data frame
    model = select(listOfModels[[m]], -significance)
    model = round(model, digits=3)
    sig = select(listOfModels[[m]], significance)
    
    #dummy indicating if the first variable in coefnames is currently treated
    first=TRUE
    
    #Go through all variables contained in any of the models' data frames.
    for (i in 1:length(coefnames)) {
      #If it is the first variable: Create a new data frame
      #containing only information on estimate, standard error and significance 
      #level for each variable which is contained in any of the models.
      #If the current model does not use the respective variable,
      #an empty line is added to the data frame. In all other cases,
      #the information is extracted from the row of the current data frame m
      #in which this information for the currently treated variable is contained
      #and written into the new data frame.
      #If it is not the first variable, the respective information/an empty line
      #is added to the data frame just described.
      if (first == TRUE){
        first = FALSE
        if (any(which(row.names(model)==coefnames[i])) == FALSE){
          results = data.frame(estimate = "", se = "", 
                               row.names = coefnames[i])
        } else {
          results = data.frame(estimate = model[which(row.names(model)==coefnames[i]),
                                                grep("Estimate+", colnames(model))],
                               se = paste0("(",as.character(model[which(row.names(model)==coefnames[i]),
                                                                  grep("SE+", colnames(model))]),")", 
                                           as.character(sig[which(row.names(model)==coefnames[i]),1])),
                               row.names = coefnames[i])
        }
      } else {
        if (any(which(row.names(model)==coefnames[i])) == FALSE){
          help = data.frame(estimate = "", se = "",
                            row.names = coefnames[i])
          results = rbind(results, help)
        } else {
          help = data.frame(estimate = model[which(row.names(model)==coefnames[i]),
                                             grep("Estimate+", colnames(model))],
                            se = paste0("(",as.character(model[which(row.names(model)==coefnames[i]),
                                                               grep("SE+", colnames(model))]),")", 
                                        as.character(sig[which(row.names(model)==coefnames[i]),1])),
                            row.names = coefnames[i])
          results = rbind(results, help)
        }
      }
    }
    #The resulting data frame for model m contains information on estimate,
    #standard error and significance level; Each row refers to one variable
    #in coefnames such that all variables appearing in any model are considered
    #in this resulting data frame.
    
    #Name the columns of the resulting data frame
    colnames(results)[(0:((ncol(results)-1)/2))*2+1] = paste0(names(listOfModels)[m], ": ",
                                                              colnames(results)[(0:((ncol(results)-1)/2))*2+1])
    colnames(results)[-((0:((ncol(results)-1)/2))*2+1)] = paste0("(",colnames(results)[-((0:((ncol(results)-1)/2))*2+1)],
                                                                 ".",as.character(m),")")
    #Create the final data frame results_df containing all resulting
    #data frames of the single models.
    if (firstModel == TRUE){
      results_df = results
      firstModel = FALSE
    } else {
      results_df = cbind(results_df, results)
    }
  }
  #Return the final data frame
  return(results_df)
}


#This function estimates the effects of each explanatory variable
#[when changin from its p_lower quantile to its p_upper quantile
#(or from 0 to 1)]
#on P(y>0|x) and E[y|x,y>0] for the Tobit model.

#input: - Tobit model estimated by command "tobit" from package 'AER'
#       - corresponding data frame
#       - probability for lower quantile
#       - probability for upper quantile
effects_tobit_dM_AER = function(tobit,data, p_lower, p_upper){
  #load package 'car' for later use of delta method
  library(car)
  
  #extract estimates from estimated Tobit model for coefficients and
  #for the logarithmized variance of the error term epsilon
  tobit_coef = summary(tobit)$coef[,1]
  
  #define number of estimates
  n = length(tobit_coef)
  
  #extract estimates of the coefficients {beta_j} from Tobit
  tobit_beta = tobit_coef[-n]
  
  #concerning all estimates (tobit_coef):
  #rename "intercept" "beta0" and "log(Sigma)" "lnsigma"
  #save renamed coefficients in renamedcoefs
  renamedcoefs = tobit_coef
  names(renamedcoefs)[c(1,n)] = c("beta0", "lnsigma")
  
  #extract explanatory variables from data
  cols = names(tobit_beta)[-1]
  X = data[,cols]
  #calculate the mean of all explanatory variables and name them appropriately
  if ( length(cols)==1){
    X_mean=mean(X)
  } else{
    X_mean=colMeans(X, na.rm=TRUE)
  }
  names(X_mean) = paste(cols, "mean", sep="")
  
  #determine binary variables and continuous variables
  firstfound = FALSE
  for (i in 1:(n-2)) {
    if (length(cols)==1){
      num_contained = unique(X)
    } else{
      num_contained = unique(X[,i])
    }
    if (length(num_contained)==2){
      if ((num_contained[1]==0 & num_contained[2]==1) | (num_contained[1]==1 & num_contained[2]==0)){
        if (firstfound==FALSE){
          firstfound = TRUE
          bin = c(i)
        } else{
          bin = c(bin, i)
        }
      }
    }
    if (length(num_contained)==1){
      if (num_contained==0 | num_contained==1){
        if (firstfound==FALSE){
          firstfound = TRUE
          bin = c(i)
        } else{
          bin = c(bin, i)
        }
      }
    }
  }
  
  conti = seq(1:(n-2))[-bin]
  
  #for each type of effect:
  #- define the formula as a string
  #- calculate estimate of effect and standard error with delta method
  #- add t value, corresponding p value and significance level
  #  to resulting data frame
  #- name the resulting data frame
  
  #effects for continuous variables
  firstconti = FALSE
  for (i in conti){
    #quantiles --> different if variable is logarithmized
    if(any(grep("ln",cols[i]))==FALSE){
      q_l = quantile(data[,cols[i]],p=p_lower)
      q_u = quantile(data[,cols[i]],p=p_upper)
      
    } else{
      q_l_real = quantile(exp(data[,cols[i]]),p=p_lower)
      q_l = log(q_l_real)
      q_u_real = quantile(exp(data[,cols[i]]),p=p_upper)
      q_u = log(q_u_real)
    }
    #according expl. variables values
    X_mean_conti_upper = X_mean
    X_mean_conti_upper[i] = q_u
    names(X_mean_conti_upper) = paste0(names(X_mean_conti_upper), "upper")
    X_mean_conti_lower = X_mean
    X_mean_conti_lower[i] = q_l
    names(X_mean_conti_lower) = paste0(names(X_mean_conti_lower), "lower")
    
    #effect on P(y>0|x)
    formula_P_conti_i = paste("pnorm((beta0 + ", paste(cols, names(X_mean_conti_upper), sep="*", collapse= " + ") ,
                              ")/exp(lnsigma))", " - ", 
                              "pnorm((beta0 + ", paste(cols, names(X_mean_conti_lower), sep="*", collapse= " + "), 
                              ")/exp(lnsigma))", sep="")
    Eff_P_conti_i = deltaMethod(object = renamedcoefs, vcov. = vcov(tobit),
                                g = formula_P_conti_i, constants = c(X_mean_conti_upper, X_mean_conti_lower))    
    Eff_P_conti_i = mutate(Eff_P_conti_i, tvalue = Estimate/SE, 
                           pvalue = 2*pt(-abs(tvalue), nrow(data) - (n-1)),
                           significance = as.factor(ifelse(pvalue<0.1, 
                                                           ifelse(pvalue<0.05,
                                                                  ifelse(pvalue<0.01, "***" , "**") , "*"), "")))
    rownames(Eff_P_conti_i)[1] = cols[i]
    colnames(Eff_P_conti_i)[1] = "Estimate P(y>0)"
    
    #effect on E[y|x,y>0]
    names(q_l)="q_l"
    names(q_u)="q_u"
    formula_E_conti_i = paste("(q_u - q_l)*", cols[i], " + exp(lnsigma) * (dnorm((beta0 + ", 
                              paste(cols, names(X_mean_conti_upper), sep="*", collapse= " + ") ,
                              ")/exp(lnsigma))/pnorm((beta0 + ",
                              paste(cols, names(X_mean_conti_upper), sep="*", collapse= " + ") ,
                              ")/exp(lnsigma)))", " - ",
                              "exp(lnsigma) * (dnorm((beta0 + ",
                              paste(cols, names(X_mean_conti_lower), sep="*", collapse= " + ") ,
                              ")/exp(lnsigma))/pnorm((beta0 + ",
                              paste(cols, names(X_mean_conti_lower), sep="*", collapse= " + ") ,
                              " - ", cols[i], ")/exp(lnsigma)))", sep="")
    Eff_E_conti_i = deltaMethod(object = renamedcoefs, vcov. = vcov(tobit),
                                g = formula_E_conti_i, 
                                constants = c(X_mean_conti_upper, X_mean_conti_lower, q_l, q_u))
    Eff_E_conti_i = mutate(Eff_E_conti_i, tvalue = Estimate/SE, 
                           pvalue = 2*pt(-abs(tvalue), nrow(data) - (n-1)),
                           significance = as.factor(ifelse(pvalue<0.1, 
                                                           ifelse(pvalue<0.05,
                                                                  ifelse(pvalue<0.01, "***" , "**") , "*"), "")))
    rownames(Eff_E_conti_i)[1] = cols[i]
    colnames(Eff_E_conti_i)[1] = "Estimate E[y|y>0]"
    
    #Save results in data frames for each kind of effect
    #of all continuous variables
    if (firstconti==FALSE){
      firstconti = TRUE
      Eff_P_conti = Eff_P_conti_i
      Eff_E_conti = Eff_E_conti_i
    } else{
      Eff_P_conti = rbind(Eff_P_conti, Eff_P_conti_i)
      Eff_E_conti = rbind(Eff_E_conti, Eff_E_conti_i)
    }
  }
  
  #effect for binary variables equal the according marginal effects
  #marginal effects for binary variables
  firstbin = FALSE
  for (i in bin){
    #adapt the mean of the explanatory variables for the calculations below
    #and save them in X_mean_bin
    X_mean_bin = X_mean
    X_mean_bin[i] = 1
    
    #marginal effect on P(y>0|x)
    formula_P_bin_i = paste("pnorm((beta0 + ", paste(cols, names(X_mean_bin), sep="*", collapse= " + ") ,
                            ")/exp(lnsigma))", " - ", 
                            "pnorm((beta0 + ", paste(cols, names(X_mean_bin), sep="*", collapse= " + "), 
                            " - ", cols[i] , ")/exp(lnsigma))", sep="")
    margEff_P_bin_i = deltaMethod(object = renamedcoefs, vcov. = vcov(tobit),
                                  g = formula_P_bin_i, constants = X_mean_bin)
    margEff_P_bin_i = mutate(margEff_P_bin_i, tvalue = Estimate/SE, 
                             pvalue = 2*pt(-abs(tvalue), nrow(data) - (n-1)),
                             significance = as.factor(ifelse(pvalue<0.1, 
                                                             ifelse(pvalue<0.05,
                                                                    ifelse(pvalue<0.01, "***" , "**") , "*"), "")))
    rownames(margEff_P_bin_i)[1] = cols[i]
    colnames(margEff_P_bin_i)[1] = "Estimate P(y>0)"
    
    #marginal effect on E[y|x,y>0]
    formula_E_bin_i = paste(cols[i], " + exp(lnsigma) * (dnorm((beta0 + ", 
                            paste(cols, names(X_mean_bin), sep="*", collapse= " + ") ,
                            ")/exp(lnsigma))/pnorm((beta0 + ",
                            paste(cols, names(X_mean_bin), sep="*", collapse= " + ") ,
                            ")/exp(lnsigma)))", " - ",
                            "exp(lnsigma) * (dnorm((beta0 + ",
                            paste(cols, names(X_mean_bin), sep="*", collapse= " + ") ,
                            " - ", cols[i], ")/exp(lnsigma))/pnorm((beta0 + ",
                            paste(cols, names(X_mean_bin), sep="*", collapse= " + ") ,
                            " - ", cols[i], ")/exp(lnsigma)))", sep="")
    margEff_E_bin_i = deltaMethod(object = renamedcoefs, vcov. = vcov(tobit),
                                  g = formula_E_bin_i, constants = X_mean_bin)
    margEff_E_bin_i = mutate(margEff_E_bin_i, tvalue = Estimate/SE, 
                             pvalue = 2*pt(-abs(tvalue), nrow(data) - (n-1)),
                             significance = as.factor(ifelse(pvalue<0.1, 
                                                             ifelse(pvalue<0.05,
                                                                    ifelse(pvalue<0.01, "***" , "**") , "*"), "")))
    rownames(margEff_E_bin_i)[1] = cols[i]
    colnames(margEff_E_bin_i)[1] = "Estimate E[y|y>0]"
    
    #Save results in data frames for each kind of marginal effect
    #of all binary variables
    if (firstbin==FALSE){
      firstbin = TRUE
      Eff_P_bin = margEff_P_bin_i
      Eff_E_bin = margEff_E_bin_i
    } else{
      Eff_P_bin = rbind(Eff_P_bin, margEff_P_bin_i)
      Eff_E_bin = rbind(Eff_E_bin, margEff_E_bin_i)
    }
    
  }
  
  #Create the data frame or warning to return
  if (firstbin ==TRUE) {
    if (firstconti == TRUE){
      return (list(Eff_P = rbind(Eff_P_bin, Eff_P_conti),
                   Eff_E = rbind(Eff_E_bin, Eff_E_conti)))
    } else {
      return (list(Eff_P = Eff_P_bin,
                   Eff_E = Eff_E_bin))
    }
  } else{
    if (firstconti == TRUE){
      return (list(Eff_P = Eff_P_conti,
                   Eff_E = Eff_E_conti))
    } else {
      return ("Error: No effects available.")
    }
  }
}