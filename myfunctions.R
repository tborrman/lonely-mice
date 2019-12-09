clean_qvals <- function(d) {
  print(nrow(d))
  # Remove NA cases
  d <- na.omit(d)
  print(nrow(d))
  print(head(d))
  
  # Fix zeros
  d[d$X712.PV.PFC_fdrDonut == 0, "X712.PV.PFC_fdrDonut"] = 1.0e-42
  d[d$X931.PV.PFC_fdrDonut == 0, "X931.PV.PFC_fdrDonut"] = 1.0e-42

  print(head(d))

  # Change q values to log10(qvals)
  d$X712.PV.PFC_fdrDonut <- log10(d$X712.PV.PFC_fdrDonut)
  d$X931.PV.PFC_fdrDonut <- log10(d$X931.PV.PFC_fdrDonut)
  print(head(d))
  
  # Remove duplicates
  d <- d[!duplicated(d[,c('chr1','x1','x2','chr2','y1','y2')]),]
  print(nrow(d))
  return(d)
  
  
}

