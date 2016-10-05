
#The most soncervative LOD value is 22 with lowest number of false positives, but trowing away true posistives. The most relaxed value LOD of 27.
LOD <- 27 

Label <- function(ct, outcome, LOD){
  if (ct > LOD ){
    return ("offscale") 
  } else if (ct < LOD & outcome == "Fail"){
    return("missing") 
  } else
    return("inrange")
}

Convert <- function(ct, label, LOD){
  if (label == "offscale") return (LOD) 
  else if (label == "missing") return(NA) 
  else return(ct)
}

#Loop labeling expression data and converiting Missing and Off Scale Data loop
df$Label  <- rep(NA, nrow(df))
df$Conv   <- rep(NA, nrow(df))

for(i in (1:nrow(df))){
  df[i,"Label"] = Label(df[i,"Ct"],
                        df[i,"Ct_Outcome"],
                        LOD)
  df[i,"Conv"] = Convert(df[i,"Ct"],
                         df[i,"Label"],
                         LOD)
}

#Computing means of Converted expressions excluding the missing gene values due to technical failure (Ct < LOD and "Fail")  Conv = NA
library(plyr)

df <- ddply(df, .(nCells, Treatment_MDR, Gene), transform, Means= mean(Conv, na.rm = "TRUE"))

df$Log2Ex <- rep(NA, nrow(df))
LogExp <- function(ct, means, label, LOD){  
 if      (label == "offscale") return (0) 
 else if (label == "missing") return (ifelse(is.na(means),NaN,LOD-means) )
 else                        return (LOD - ct)
}

for(i in (1:nrow(df))){
  df[i,"Log2Ex"] = LogExp(df[i,"Ct"], 
                            df[i,"Means"],
                            df[i,"Label"], 
                            LOD)
}

#Average technical replicates see the SingleCellqPCR.r

