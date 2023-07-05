status_matrices_f <- function(df, readers_x_methods) {

#Total unit counts divided by status
groups_status <- df %>% group_by(reference, reader, method) %>% summarise(count = n()) %>%
  group_by(reference) %>% group_split()

groups_status_full <- list()

for (i in 1:length(groups_status)) {

  checks <- data.frame("check" = do.call(paste0, readers_x_methods) %in% do.call(paste0, groups_status[[i]][,2:3]))
  completed <- checks %>% group_by(check) %>% mutate(count = case_when(check == TRUE ~ 1, check == FALSE ~ 0))
  completed$count[completed$count == 1] <- groups_status[[i]]$count
  final <- data.frame(groups_status[[i]] %>% slice(rep(row_number(1), length(checks[,1]))),completed$count)
  final[,2:3] <- readers_x_methods
  groups_status_full[[i]] <- final

}

if (length(groups_status_full[[1]]$completed.count)> 1) {
groups_status_full[[1]] <- diag(groups_status_full[[1]]$completed.count)
groups_status_full[[2]] <- diag(groups_status_full[[2]]$completed.count)
} else {
  groups_status_full[[1]] <- groups_status_full[[1]]$completed.count
  groups_status_full[[2]] <- groups_status_full[[2]]$completed.count
}

names(groups_status_full) <- c("healthy", "diseased")

return(groups_status_full)
}
