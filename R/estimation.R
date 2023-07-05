#' Main function of the diagacc package
#'
#' `estimation` returns the diagnostic accuracy estimates of the package
#'
#' A typical function call will be of the form
#' estimation(df, ID = "Patient", index = "test", reference = "reference", logit_transform = TRUE)
#' where df is the used dataframe. ID is a reference variable for subjects. index is the test for which
#' accuracy will be estimated. reference is the gold standard against which the index test is compared.
#' logit_transform is an optional transform of the model to ensure estimates being in the unit interval (recommended).
#' cutoff is the diagnostic cutoff value which, if specified, will be used to estimate sensitivity and specificity.
#' reader and method are (not yet functional!) reference variables for factorial designs.

estimation <- function(df, ID, index, reference, cutoff, reader, method, logit_transform) {

#Sets default value of 0.5 for the diagnostic cutoff if the argument is missing
if(missing(cutoff)) {
   cutoff <- 0.5
}

#Adds reader column if it is missing, otherwise renames the variable
if (missing(reader)) {
  df$reader <- 1
} else {
  colnames(df)[colnames(df) == reader] <- 'reader'
}

#Adds method column if it is missing, otherwise renames the variable
if (missing(method)) {
  df$method <- 1
} else {
  colnames(df)[colnames(df) == method] <- 'method'
}

#Creates reader/method data.frame to reference factors for the evaluation
readers_x_methods <- data.frame(expand.grid(unique(df$reader), unique(df$method)))
colnames(readers_x_methods) <- c("reader", "method")

#Renames custom names of index and reference test names
colnames(df)[colnames(df) == ID] <- 'ID'
colnames(df)[colnames(df) == index] <- 'index'
colnames(df)[colnames(df) == reference] <- 'reference'


#Adds identifier column for complete and incomplete cases
df <- df %>% group_by(ID) %>%
  mutate(identifier = case_when(min(reference) == max(reference) ~ 0, min(reference) != max(reference) ~ 1)) %>%
  ungroup()


  ##AUC##
  result_AUC <- AUC_procedure(df)
  result_AUC <- result_AUC$estimates

  ##Sensitivity##
  df_sens <- df
  df_sens$index[df_sens$reference == 0] <- cutoff
  result_sens <- AUC_procedure(df_sens)$estimates



  ##Specificity##
  df_spec <- df
  df_spec$index[df_spec$reference == 1] <- cutoff
  result_spec <- AUC_procedure(df_spec)$estimates

  ##Predictive Values##
  df_spec_alt <- df_spec
  df_spec_alt$reader <- df_spec_alt$reader*10000
  df_spec_alt$method <- df_spec_alt$method*10000
  df_both <- rbind(df_sens, df_spec_alt)

  result_pv <- AUC_procedure(df_both)$estimates
  result_sens_pv <- c(result_pv[1:(length(result_pv)/2)])
  result_spec_pv <- c(result_pv[(length(result_pv)/2 + 1):length(result_pv)])


  prev <- length(df$reference[df$reference == 1])/length(df$reference)

  PPV <- (result_sens_pv*prev)/((result_sens_pv*prev) + (1-result_spec_pv)*(1-prev))
  NPV <- (result_spec_pv*(1-prev))/(result_spec_pv*(1-prev)+(1-result_sens_pv)*prev)



  patient_matrices   <- patient_units_f(df, readers_x_methods)
  groups_status_full <- status_matrices_f(df, readers_x_methods)
  c_s_matrices       <- case_status_matrices_f(df, readers_x_methods)

  auc_cov  <- cov_matrix(df, groups_status_full, c_s_matrices, patient_matrices, readers_x_methods)
  sens_cov <- cov_matrix(df_sens, groups_status_full, c_s_matrices, patient_matrices, readers_x_methods)
  spec_cov <- cov_matrix(df_spec, groups_status_full, c_s_matrices, patient_matrices, readers_x_methods)


  df_both <- rbind(df_sens, df_spec)
  patient_matrices   <- patient_units_f(df_both, readers_x_methods)
  groups_status_full <- status_matrices_f(df_both, readers_x_methods)
  c_s_matrices       <- case_status_matrices_f(df_both, readers_x_methods)
  pv_cov   <- cov_matrix(df_both, groups_status_full, c_s_matrices, patient_matrices, readers_x_methods)


if (logit_transform == FALSE) {

  result  <-  data.frame(result_AUC,
                         result_AUC - sqrt(1/length(unique(df$ID)) * auc_cov) * qt(p=1 - (1-0.95)/2, df=length(unique(df$ID))-1),
                         result_AUC + sqrt(1/length(unique(df$ID)) * auc_cov) * qt(p=1 - (1-0.95)/2, df=length(unique(df$ID))-1))

  result[2,] <- c(result_sens,
                  result_sens - sqrt(1/length(unique(df$ID)) * sens_cov) * qt(p=1 - (1-0.95)/2, df=length(unique(df$ID))-1),
                  result_sens + sqrt(1/length(unique(df$ID)) * sens_cov) * qt(p=1 - (1-0.95)/2, df=length(unique(df$ID))-1))

  result[3,] <- c(result_spec,
                  result_spec - sqrt(1/length(unique(df$ID)) * spec_cov) * qt(p=1 - (1-0.95)/2, df=length(unique(df$ID))-1),
                  result_spec + sqrt(1/length(unique(df$ID)) * spec_cov) * qt(p=1 - (1-0.95)/2, df=length(unique(df$ID))-1))

#  result[4,] <- c(PPV,
#                  PPV - sqrt(1/length(unique(df_both$ID)) * pv_cov) * qt(p=1 - (1-0.95)/2, df=length(unique(df_both$ID))-1),
#                  PPV + sqrt(1/length(unique(df_both$ID)) * pv_cov) * qt(p=1 - (1-0.95)/2, df=length(unique(df_both$ID))-1))

#  result[5,] <- c(NPV,
#                  NPV - sqrt(1/length(unique(df$ID)) * pv_cov) * qt(p=1 - (1-0.95)/2, df=length(unique(df$ID))-1),
#                  NPV + sqrt(1/length(unique(df$ID)) * pv_cov) * qt(p=1 - (1-0.95)/2, df=length(unique(df$ID))-1))

} else {

  auc_cov_logit  <- logit_jacobian(result_AUC)%*%auc_cov%*%logit_jacobian(result_AUC)
  sens_cov_logit <- logit_jacobian(result_sens)%*%sens_cov%*%logit_jacobian(result_sens)
  spec_cov_logit <- logit_jacobian(result_spec)%*%spec_cov%*%logit_jacobian(result_spec)

  result <- data.frame(result_AUC,
                  expit(logit(result_AUC) - (sqrt(auc_cov_logit) * qt(p=1 - (1-0.95)/2, df=length(unique(df$ID))-1))/sqrt(length(unique(df$ID)))),
                  expit(logit(result_AUC) + (sqrt(auc_cov_logit) * qt(p=1 - (1-0.95)/2, df=length(unique(df$ID))-1))/sqrt(length(unique(df$ID)))))

  result[2,] <- c(result_sens,
                  expit(logit(result_sens) - (sqrt(sens_cov_logit) * qt(p=1 - (1-0.95)/2, df=length(unique(df$ID))-1))/sqrt(length(unique(df$ID)))),
                  expit(logit(result_sens) + (sqrt(sens_cov_logit) * qt(p=1 - (1-0.95)/2, df=length(unique(df$ID))-1))/sqrt(length(unique(df$ID)))))

  result[3,] <- c(result_spec,
                  expit(logit(result_spec) - (sqrt(spec_cov_logit) * qt(p=1 - (1-0.95)/2, df=length(unique(df$ID))-1))/sqrt(length(unique(df$ID)))),
                  expit(logit(result_spec) + (sqrt(spec_cov_logit) * qt(p=1 - (1-0.95)/2, df=length(unique(df$ID))-1))/sqrt(length(unique(df$ID)))))

#  result[4,] <- c(PPV,
#                  PPV - sqrt(1/length(unique(df_both$ID)) * pv_cov) * qt(p=1 - (1-0.95)/2, df=length(unique(df_both$ID))-1),
#                  PPV + sqrt(1/length(unique(df_both$ID)) * pv_cov) * qt(p=1 - (1-0.95)/2, df=length(unique(df_both$ID))-1))

#  result[5,] <- c(NPV,
#                  NPV - sqrt(1/length(unique(df$ID)) * pv_cov) * qt(p=1 - (1-0.95)/2, df=length(unique(df$ID))-1),
#                  NPV + sqrt(1/length(unique(df$ID)) * pv_cov) * qt(p=1 - (1-0.95)/2, df=length(unique(df$ID))-1))

}




  colnames(result) <- c("Estimate", "Lower", "Upper")
  rownames(result) <- c("AUC", "Sensitivity", "Specificity")#, "PPV", "NPV")
  return(result)

}

