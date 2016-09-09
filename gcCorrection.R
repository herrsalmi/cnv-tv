###############################
## Gc content bias correction
###############################
require(GenomicAlignments)

gc.norm <- function(ref, cvg, method = "window"){
  window.size <- 100
  gc <- gcContent(as.character(ref[[1]]), 100)
  pad.right <- function (rle, value=0, len) {
    if (length(rle) > len) stop("length mismatch: len must be >= length(rle)!")
    if (len > length(rle))
      c(rle, Rle(value, len-length(rle)))
    else rle
  }

  if (length(cvg) > length(gc)){
    cvgP <- cvg[1:length(gc)]
  } else {
    #cvgP <- pad.right(cvg, len=length(gc))
    cvgP <- cvg
    gc <- gc[-length(gc)]
  }
  
  gcU <- unique(gc)
  
  nonZeroIndex <- which(cvgP != 0)
  zeroList <- 1:min(nonZeroIndex)
  
  binsMedian <- numeric(length(gcU))
  names(binsMedian) <- gcU
  for (i in 1:length(gcU)) {
    index <- which(gc == gcU[i])
    index <- index[!index %in% zeroList]
    if (length(index) == 0) {
      binsMedian[i] <- 1
      next()
    }
    binsMedian[i] <- median(cvgP[index])
  }
  
  d <- sqrt(median((cvgP[min(nonZeroIndex):length(cvgP)])^2))
  
  adjustedReadCount <- as.numeric(cvgP) * (d / binsMedian[as.character(gc)])
  return(adjustedReadCount)
}

###############################
## LOWESS correction
###############################

# lowess.gc <- function(jtkx, jtky) {
#   jtklow <- lowess(jtkx, log(jtky), f=0.05)
#   jtkz <- approx(jtklow$x, jtklow$y, jtkx)
#   return(exp(log(jtky) - jtkz$y))
# }
# 
# a <- cvg + 1  # to avoid log(0) when we apply log transformation
# a <- a / mean(a)
# gc.content <- gcContent(as.character(ref[[1]]), 100)
# lowratio <- lowess.gc(gc.content[-49390], windowAll[[2]])
# 
