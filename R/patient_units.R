patient_units_f <- function(df, readers_x_methods) {

  #splits df by ID, identifier and reference with unit counts for every factor added
  groups <- df %>% group_by(ID, identifier, reference, reader, method) %>% summarise(count = n()) %>%
    group_by(reference, ID, identifier) %>% group_split()


  groups_adapted <- list()


  #checks for missing factors in patients and adds the row with a unit count of 0
  for (i in 1:length(groups)) {

    checks <- data.frame("check" = do.call(paste0, readers_x_methods) %in% do.call(paste0, groups[[i]][,4:5]))
    completed <- checks %>% group_by(check) %>% mutate(count = case_when(check == TRUE ~ 1, check == FALSE ~ 0))
    completed$count[completed$count == 1] <- groups[[i]]$count
    final <- data.frame(groups[[i]] %>% slice(rep(row_number(1), length(checks[,1]))),completed$count)
    final[,4:5] <- readers_x_methods
    groups_adapted[[i]] <- final

  }

counts <- groups_adapted

  comp_health <- list()
  comp_dis    <- list()
  incomp_health <- list()
  incomp_dis    <- list()

  i <- 1
  j <- 1
  k <- 1
  l <- 1

  for (n in 1:length(counts)) {
    if (counts[[n]]$identifier[1] == 1 & counts[[n]]$reference[1] == 0) {
      comp_health[[i]] <- counts[[n]]
      i <- i+1
    } else if (counts[[n]]$identifier[1] == 1 & counts[[n]]$reference[1] == 1) {
      comp_dis[[j]] <- counts[[n]]
      j <- j+1
    } else if (counts[[n]]$identifier[1] == 0 & counts[[n]]$reference[1] == 0) {
      incomp_health[[k]] <- counts[[n]]
      k <- k+1
    } else if (counts[[n]]$identifier[1] == 0 & counts[[n]]$reference[1] == 1) {
      incomp_dis[[l]] <- counts[[n]]
      l <- l+1
    }
  }

count_list <- list("comp_health" = comp_health, "comp_dis" = comp_dis, "incomp_health" = incomp_health, "incomp_dis" = incomp_dis)
}
