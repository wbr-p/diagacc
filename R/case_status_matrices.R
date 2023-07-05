case_status_matrices_f <- function(df, readers_x_methods) {

#status and case diagonal matrices
groups_case <- df %>% group_by(reference, identifier, reader, method) %>% summarise(count = n()) %>%
  group_by(reference, identifier) %>% group_split()


groups_case_full <- list()

for (i in 1:length(groups_case)) {

  checks <- data.frame("check" = do.call(paste0, readers_x_methods) %in% do.call(paste0, groups_case[[i]][,3:4]))
  completed <- checks %>% group_by(check) %>% mutate(count = case_when(check == TRUE ~ 1, check == FALSE ~ 0))
  completed$count[completed$count == 1] <- groups_case[[i]]$count
  final <- data.frame(groups_case[[i]] %>% slice(rep(row_number(1), length(checks[,1]))),completed$count)
  final[,3:4] <- readers_x_methods
  groups_case_full[[i]] <- final

}

#Naming ((to be optimized!!))
t1 <- map(groups_case_full, ~filter(.x, reference == 0 & identifier == 0))
t1 <- t1[which(lapply(t1, nrow) != 0)]

t2 <- map(groups_case_full, ~filter(.x, reference == 1 & identifier == 0))
t2 <- t2[which(lapply(t2, nrow) != 0)]

t3 <- map(groups_case_full, ~filter(.x, reference == 0 & identifier == 1))
t3 <- t3[which(lapply(t3, nrow) != 0)]

t4 <- map(groups_case_full, ~filter(.x, reference == 1 & identifier == 1))
t4 <- t4[which(lapply(t4, nrow) != 0)]


groups_case_comp <- list("incomp_health" = t1, "incomp_dis" = t2, "comp_health" = t3, "comp_dis" = t4)


if (length(t1) ==1) {
  if (length(t1[[1]]$completed.count) > 1) {
groups_case_comp[[1]] <- diag(t1[[1]]$completed.count)
  } else {
    groups_case_comp[[1]] <- t1[[1]]$completed.count
  }
} else { groups_case_comp[[1]] <- 0}

if (length(t2) ==1) {
  if (length(t2[[1]]$completed.count) > 1) {
groups_case_comp[[2]] <- diag(t2[[1]]$completed.count)
} else {
  groups_case_comp[[2]] <- t2[[1]]$completed.count}
} else { groups_case_comp[[2]] <- 0}


if (length(t3) ==1) {
  if (length(t3[[1]]$completed.count) > 1) {
  groups_case_comp[[3]] <- diag(t3[[1]]$completed.count)
} else {
  groups_case_comp[[3]] <- t3[[1]]$completed.count}
} else { groups_case_comp[[3]] <- 0}


if (length(t4) ==1) {
  if (length(t4[[1]]$completed.count) > 1) {
    groups_case_comp[[4]] <- diag(t4[[1]]$completed.count)
  } else {
    groups_case_comp[[4]] <- t4[[1]]$completed.count}
} else { groups_case_comp[[4]] <- 0}


return(groups_case_comp)
}
