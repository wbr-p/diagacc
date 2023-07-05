data_generation <- function(prevalence, ID, reader, method, cluster_size) {
  
  #Kartesisches Produkt von {1,..,Anzahl Reader} und {1,...,Anzahl Methods}.
  readers_x_methods <- data.frame(crossing(reader,method))
  
  df       <- data.frame(ID) %>% slice(rep(1:n(), each = cluster_size))
  df[,2]   <- rbinom(length(df[,1]),1,prevalence)
  df       <- df %>% slice(rep(1:n(), each = length(readers_x_methods[,1])))
  df[,3]   <- rep(readers_x_methods[,1], cluster_size*length(ID))
  df[,4]   <- rep(readers_x_methods[,2], cluster_size*length(ID))
  df[,5]   <- rnorm(df$ID, 0, 1)
  colnames(df) <- c("ID", "reference", "reader", "method", "index")
  
  return(df)
  
}