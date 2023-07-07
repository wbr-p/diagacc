cov_matrix <- function(df, groups_status_full, c_s_matrices, patient_matrices, readers_x_methods) {

  #Global Ranksums
  glob_rank_sums <- df %>% group_by(reader, method) %>% mutate(ranks = rank(index)) %>%
    group_by(ID, identifier, reference, reader, method) %>% summarise(ranksum = sum(ranks)) %>% ungroup()

  #Internal ranksums
  int_rank_sums  <- df %>% group_by(reference, reader, method) %>% mutate(ranks = rank(index)) %>%
    group_by(ID, identifier, reference, reader, method) %>% summarise(ranksum = sum(ranks)) %>% ungroup()

  #Difference of global and internal ranksums by case and summed up differences
  diff_complete     <- cbind(glob_rank_sums, glob_rank_sums$ranksum - int_rank_sums$ranksum)
  colnames(diff_complete)[7] <- "difference"
  diff_complete_sum <- diff_complete %>% group_by(identifier, reference, reader, method) %>% summarise(sums = sum(difference)) %>% ungroup()

  diff_comp <- diff_complete[diff_complete$identifier == 1,]
  diff_incomp <- diff_complete[diff_complete$identifier == 0,]

  #number of cases
  comp_case <- df %>% group_by(identifier, reference) %>% summarise(count = n_distinct(ID))

  #####################################################################################################################

  if (length(diff_comp$ID) > 0) {

    #Weight of Covariance Matrix of complete cases
    covma_comp <- (length(unique(df$ID))*comp_case$count[comp_case$identifier == 1 & comp_case$reference ==1])/
      (comp_case$count[comp_case$identifier == 1 & comp_case$reference ==1]-1)*
      ((solve(groups_status_full$healthy)*
          solve(groups_status_full$diseased))**2)

    #initialising covariance matrix
    mat = matrix(0, length(readers_x_methods[,1]), length(readers_x_methods[,1]))

    i <- 1
    #comp covariance matrix
    for (value in unique(diff_comp$ID)) {

      mat <- mat +
        (diff_comp$difference[diff_comp$ID == value & diff_comp$reference == 1] -
           diff_comp$difference[diff_comp$ID == value & diff_comp$reference == 0] -
           (diag(patient_matrices$comp_dis[[i]]$completed.count, nrow = length(patient_matrices$comp_dis[[i]]$completed.count), ncol = length(patient_matrices$comp_dis[[i]]$completed.count))*
              solve(c_s_matrices$comp_dis)*
              diff_complete_sum$sums[diff_complete_sum$identifier == 1 & diff_complete_sum$reference == 1] -
              diag(patient_matrices$comp_health[[i]]$completed.count, nrow = length(patient_matrices$comp_health[[i]]$completed.count), ncol = length(patient_matrices$comp_health[[i]]$completed.count))*
              solve(c_s_matrices$comp_health)*
              diff_complete_sum$sums[diff_complete_sum$identifier == 1 & diff_complete_sum$reference == 0]))*
        t((diff_comp$difference[diff_comp$ID == value & diff_comp$reference == 1] -
             diff_comp$difference[diff_comp$ID == value & diff_comp$reference == 0] -
             (diag(patient_matrices$comp_dis[[i]]$completed.count, nrow = length(patient_matrices$comp_dis[[i]]$completed.count), ncol = length(patient_matrices$comp_dis[[i]]$completed.count))*
                solve(c_s_matrices$comp_dis)*
                diff_complete_sum$sums[diff_complete_sum$identifier == 1 & diff_complete_sum$reference == 1] -
                diag(patient_matrices$comp_health[[i]]$completed.count, nrow = length(patient_matrices$comp_health[[i]]$completed.count), ncol = length(patient_matrices$comp_health[[i]]$completed.count))*
                solve(c_s_matrices$comp_health)*
                diff_complete_sum$sums[diff_complete_sum$identifier == 1 & diff_complete_sum$reference == 0])))
      i <- i+1
    }

    covma_comp <- covma_comp[1,1]*mat
    covma_comp

  } else {

    covma_comp <- matrix(0, length(readers_x_methods[,1]), length(readers_x_methods[,1]))
  }

  #####################################################################################################################

  if (length(diff_incomp$ID[diff_incomp$reference == 1]) > 0) {

    #Weight of Covariance Matrix of diseased incomplete cases
    covma_incomp_dis <- (length(unique(df$ID))*comp_case$count[comp_case$identifier == 0 & comp_case$reference ==1])/
      (comp_case$count[comp_case$identifier == 0 & comp_case$reference ==1]-1)*
      ((solve(groups_status_full$healthy)*
          solve(groups_status_full$diseased))**2)

    #initialising covariance matrix
    mat = matrix(0, length(readers_x_methods[,1]), length(readers_x_methods[,1]))

    i <- 1
    #incomp diseased covariance matrix
    for (value in unique(diff_incomp$ID[diff_incomp$reference == 1])) {

      mat <- mat +
        ((diff_incomp$difference[diff_incomp$ID == value & diff_incomp$reference == 1] -
           (diag(patient_matrices$incomp_dis[[i]]$completed.count, nrow = length(patient_matrices$incomp_dis[[i]]$completed.count), ncol = length(patient_matrices$incomp_dis[[i]]$completed.count))*
              solve(c_s_matrices$incomp_dis)*
              diff_complete_sum$sums[diff_complete_sum$identifier == 0 & diff_complete_sum$reference == 1]))*
           t(diff_incomp$difference[diff_incomp$ID == value & diff_incomp$reference == 1] -
               (diag(patient_matrices$incomp_dis[[i]]$completed.count, nrow = length(patient_matrices$incomp_dis[[i]]$completed.count), ncol = length(patient_matrices$incomp_dis[[i]]$completed.count))*
                  solve(c_s_matrices$incomp_dis)*
                  diff_complete_sum$sums[diff_complete_sum$identifier == 0 & diff_complete_sum$reference == 1])))
     i <- i+1
    }

    covma_incomp_dis <- covma_incomp_dis[1,1]*mat
    covma_incomp_dis

  } else {

    covma_incomp_dis <- matrix(0, length(readers_x_methods[,1]), length(readers_x_methods[,1]))
  }

  #####################################################################################################################

  if (length(diff_incomp$ID[diff_incomp$reference == 0]) > 0) {

    #Weight of Covariance Matrix of healthy incomplete cases
    covma_incomp_health <- (length(unique(df$ID))*comp_case$count[comp_case$identifier == 0 & comp_case$reference == 0])/
      (comp_case$count[comp_case$identifier == 0 & comp_case$reference == 0]-1)*
      ((solve(groups_status_full$healthy)*
          solve(groups_status_full$diseased))**2)

    #initialising covariance matrix
    mat = matrix(0, length(readers_x_methods[,1]), length(readers_x_methods[,1]))

    i <- 1
    #incomp healthy covariance matrix
    for (value in unique(diff_incomp$ID[diff_incomp$reference == 0])) {

      mat <- mat +
        ((diff_incomp$difference[diff_incomp$ID == value & diff_incomp$reference == 0] -
           (diag(patient_matrices$incomp_health[[i]]$completed.count, nrow = length(patient_matrices$incomp_health[[i]]$completed.count), ncol = length(patient_matrices$incomp_health[[i]]$completed.count))*
              solve(c_s_matrices$incomp_health)*
              diff_complete_sum$sums[diff_complete_sum$identifier == 0 & diff_complete_sum$reference == 0]))*
           t(diff_incomp$difference[diff_incomp$ID == value & diff_incomp$reference == 0] -
               (diag(patient_matrices$incomp_health[[i]]$completed.count, nrow = length(patient_matrices$incomp_health[[i]]$completed.count), ncol = length(patient_matrices$incomp_health[[i]]$completed.count))*
                  solve(c_s_matrices$incomp_health)*
                  diff_complete_sum$sums[diff_complete_sum$identifier == 0 & diff_complete_sum$reference == 0])))
      i<-i+1
    }

    covma_incomp_health <- covma_incomp_health[1,1]*mat
    covma_incomp_health

  } else {

    covma_incomp_health <- matrix(0, length(readers_x_methods[,1]), length(readers_x_methods[,1]))
  }

  covma_comp[is.na(covma_comp)] <- 0
  covma_incomp_dis[is.na(covma_incomp_dis)] <- 0
  covma_incomp_health[is.na(covma_incomp_health)] <- 0



  cov_matrix <- covma_comp + covma_incomp_dis + covma_incomp_health

  return(cov_matrix)

}
