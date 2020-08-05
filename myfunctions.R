clean_qvals <- function(d) {
  print(nrow(d))
  # Remove NA cases
  d <- na.omit(d)
  print(nrow(d))
  print(head(d))
  
  # Fix zeros
  d[d$GH_fdrDonut == 0, "GH_fdrDonut"] = 1.0e-42
  d[d$jSI_fdrDonut == 0, "jSI_fdrDonut"] = 1.0e-42

  print(head(d))

  # Change q values to log10(qvals)
  d$GH_fdrDonut <- log10(d$GH_fdrDonut)
  d$jSI_fdrDonut <- log10(d$jSI_fdrDonut)
  print(head(d))
  
  # Remove duplicates
  d <- d[!duplicated(d[,c('chr1','x1','x2','chr2','y1','y2')]),]
  print(nrow(d))
  return(d)
  
  
}

