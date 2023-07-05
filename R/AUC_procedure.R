AUC_procedure <- function(df) {
  ##Relevant counting values##
  
  #Overall numbers of units under specific condition
  unit_all <- df %>% group_by(reader,method) %>% summarise(unit_count = n()) %>% ungroup()
  
  #Sum of units in status i
  unit_by_status <- df %>% group_by(reader, method, reference) %>% summarise(unit_count = n()) %>% ungroup()
  
  #Sum of units by case, status and ID
  unit_by_ID <- df %>% group_by(reader, method, ID) %>% summarise(unit_count = n()) %>% ungroup()
  
  #Numbers of complete and incomplete cases
  complete_and_incomplete <- df %>% group_by(ID) %>% summarise(case = min(reference)+max(reference)) %>%
    group_by(case) %>% summarise(count = n()) %>% ungroup()
  
  #Adds identifier for complete and incomplete cases to dataframe
  df <- df %>% group_by(ID) %>% mutate(case = min(reference)+max(reference)) %>% ungroup()
  
  #Sum of units by case, status and ID
  unit_by_ID <- df %>% group_by(reader, method, reference, case, ID) %>% summarise(unit_count = n()) %>% ungroup()
  
  ## Function to calculate AUC and by extension sensitivity and specificity ##
  
  
  #Sums ranks of within a condition, status, case and ID. Then sums those within condition, status and case. Again sums
  #those within condition and status.
  weighted_ranks          <- df %>% group_by(reader,method) %>% mutate(rank = rank(index)) %>%
    group_by(reference, case, ID, .add = TRUE) %>% summarise(unit_sum = sum(rank)) %>%
    group_by(reader, method, reference, case) %>% summarise(id_sum = sum(unit_sum)) %>%
    group_by(reader, method, reference) %>% summarise(case_sum = sum(id_sum)) %>% ungroup()
  
  #Calculates the weighted mean rank.
  weighted_ranks$case_sum <- (1/unit_by_status$unit_count)*weighted_ranks$case_sum
  
  #Calculates the difference between the weighted ranks within status groups.
  
  if (prod(weighted_ranks$reference) == 0) {
  
  weighted_ranks          <- weighted_ranks %>% group_by(reader,method) %>% 
    summarise(weighted_ranks = case_sum[reference == 1]-case_sum[reference == 0]) %>% ungroup()
  } else {
  weighted_ranks <- weighted_ranks[,-3]
  weighted_ranks <- weighted_ranks %>% mutate(weighted_ranks = case_sum)
  }
  
  #Calculates the AUC 
  AUC <- 1/(unit_all$unit_count)*(weighted_ranks$weighted_ranks) + 1/2
  result <- cbind(AUC)
  
  return(list("estimates" = result, "All units counted" = unit_all, "All units by status" = unit_by_status, "complete and incomplete cases" = complete_and_incomplete, "Unit by ID" = unit_by_ID))
}