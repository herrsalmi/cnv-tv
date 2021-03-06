#require(ggplot2)
#library(plotly)
library(GenomicAlignments)
library(genlasso)
source("gcCorrection.R")
Rcpp::sourceCpp('gc.content.cpp')

binSize<-100
windowSize <- 0.5e+6
prefix <- "NA07037"
upper.thr <- numeric()
lower.thr <- numeric()
ref <- DNAStringSet()

#' OS type
#' 
#' Detects OS type to set null device used for messages redirection
#' 
#' @return OS type
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
    if (grepl("mingw32", R.version$os))
      os <- "windows"
  }
  tolower(os)
}

## Set null device
if (get_os() == "windows") {
  NULL.DEV <- "NUL"
} else {
  NULL.DEV <- "/dev/null"
}

#' Mean read depth over a sliding windows
#' 
#' @param windowsize sliding window size
#' @param inputseq vector of read depth per base
#' @return list with positions, mean depth and sd
slidingwindowcoverage <- function(windowsize, inputseq) {
  starts <- seq(1, length(inputseq)-windowsize, by=windowsize)
  n <- length(starts)
  chunkbps <- numeric(n)
  chunkstats<- numeric(n)
  for (i in 1:n) {
    chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)] 
    chunkmean <- mean(chunk)
    chunkstdv <- sd(chunk)
    chunkbps[i] <- chunkmean
    chunkstats[i] <- chunkstdv
  }
  return (list(starts,chunkbps,chunkstats))
}

#' Finds segments with uniform intensities
#' 
#' @param x vector of intensities
#' @return list of (starts,ends) of segments
findSegment <- function(x) {
  starts <- 1
  ends <- c()
  for (i in seq(2, length(x))) {
    if (round(x[i-1], 6) == round(x[i], 6)) {
      next()
    } else {
      starts <- c(starts, i)
      ends <- c(ends, i-1)
    }
  }
  ends <- c(ends, length(x))
  return(cbind(starts, ends))
}

#' Groups intensities with sequential indexes
#' 
#' @param x indexes from list of dup/del
#' @return list of (starts,ends) of segments
group <- function(x){
  starts <- x[1]
  ends <- c()
  for (i in seq(2, length(x))) {
    if (x[i-1] + 1 == x[i]) {
      next()
    } else {
      starts <- c(starts, x[i])
      ends <- c(ends, x[i-1])
    }
  }
  ends <- c(ends, x[length(x)])
  return(cbind(starts, ends))
}

#' Computes coverage for each base from a BAM file
#' 
#' @param file BAM file
#' @return coverage values for each base
extractCoverageFromBAM <- function(file) {
  param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE,
                                         isDuplicate=FALSE),
                        what=c("rname", "pos", "cigar"))
  bam <- scanBam(file, param=param)[[1]]
  ## Note that unmapped reads and reads that are PCR/optical duplicates
  ## have already been filtered out by using the ScanBamParam object above.
  irl <- extractAlignmentRangesOnReference(bam$cigar, pos=bam$pos,
                                           f=bam$rname)
  irl <- irl[elementNROWS(irl) != 0] # drop empty elements
  coverage(irl)
}

#' Mean coverage as data.frame
#' 
#' Calls slidingwindowcoverage and returns a formated data.farme
#' 
#' @param cvg coverage values for each base
#' @return data.frame of mean coverage
computeMeanCoverage <- function(cvg){
  myvector_all <- as.vector(cvg)
  windowAll <- slidingwindowcoverage(binSize, myvector_all)
  df <- data.frame(windowAll[[1]],windowAll[[2]],windowAll[[3]])
  colnames(df) <- c("x","mean","sd")
  return(df)
}

#' Detects CNVs in window
#' 
#' Detects CNVs in window defined by (start,end) using LASSO solver
#' and sets thresholds for CNV calling by using a lognormal distribution
#' 
#' @param data mean coverage
#' @param start start position of the window
#' @param end end position of the window
#' @param chrom.name chromosome name
#' @return data.frame of dup/del positions in window
get.cnv <- function(data, start, end, chrom.name){
  # Skip 500Bp from the begening and end of the chromosom
  # if (start < 500)
  #   start <- 501
  # if (end > tail(data$x, n=1)-500)
  #   end <- tail(data$x, n=1)-500
  y <- data[data$x >= start & data$x <= end,]
  y <- y[is.finite(rowSums(y)),]
  if(length(y$mean) == 0)
    return(NULL)
  r <- rle(y$mean)
  if(max(r$length[r$values == 0]) > windowSize/2)
    return(NULL)
  out <- fusedlasso1d(y$mean, y$x)
  sink(file = NULL.DEV)
  cv <- cv.trendfilter(out)
  sink()
  x <- coef(out, lambda = cv$lambda.1se)
  x <- x$beta
  seg <- findSegment(x)
  names(seg) <- NULL
  ampl <- c()
  for (i in 1:length(seg[,1])) {
    ampl <- c(ampl, sum(y$mean[seg[i,1]:seg[i,2]]) / (seg[i,2] - seg[i,1] + 1) )
  }
  seg.length <- (seg[,2] - seg[,1] + 1)
  x.t <- rep(ampl, seg[,2] - seg[,1] + 1)
  names(x.t) <- NULL
  svg(paste0("Plots_", prefix,"/", chrom.name, "_", start, "_", end, ".svg"),width=14,height=7)
  plot(out, lambda = cv$lambda.1se, pch = ".", cex = 2, col = 'blue')
  lines(y$x, x.t, col = "red")
  abline(a = lower.thr, b = 0, col = "green")
  abline(a = upper.thr, b = 0, col = "green")
  dev.off()
  type <- vector(mode = "character")
  type[x.t >= upper.thr] <- "dup"
  type[x.t <= lower.thr] <- "del"
  x.t[x.t < upper.thr & x.t > lower.thr] <- 0
  return(data.frame(intensity = x.t, pos = y$x, type = type))
}

#' Processes a chromosome for CNVs using sliding windows of pre-defined size
#' 
#' @param chrom coverage over chromosome
#' @param chrom.name chromosome name
extract.cnv <- function(chrom, chrom.name) {
  starts <- seq(1, length(chrom)-windowSize, by=windowSize)
  starts <- c(starts, length(chrom))
  n <- length(starts) - 1
  data <- computeMeanCoverage(chrom)
  # GC correction
  data$mean <- gc.norm(ref, data$mean)
  # compute thresholds
  f <- data$mean[data$mean > 0]
  f <- f[is.finite(f)]
  params <- fitdistr(f, "lognormal")
  upper.thr <<- qlnorm(0.975, params$estimate['meanlog'], params$estimate['sdlog'])
  lower.thr <<- qlnorm(0.012, params$estimate['meanlog'], params$estimate['sdlog'])
  
  print(paste0("Processing chromosome ", chrom.name))
  # fixe the pb: min=1
  pb <- txtProgressBar(min = 1, max = n, style = 3)
  for (i in 1:n) {
    #print(paste0("Processing ", chrom.name, ":", starts[i], "-",starts[i+1]))
    setTxtProgressBar(pb, i)
    cnv.list <- get.cnv(data, starts[i], starts[i+1], chrom.name)
    if(is.null(cnv.list) || length(cnv.list) == 0) {
      next()
    }
    i <- which(!is.na(cnv.list$type))
    if (length(i) < 2) {
      next()
    }
    g <- group(i)
    intensities <- vector(mode = "numeric")
    for (index in 1:(length(g)/2)) {
      intensities <- c(intensities, max(cnv.list$intensity[g[index,1]:g[index,2]]))
    }
    cnv <- data.frame("Chrom" = chrom.name, "Start" = cnv.list$pos[g[,1]], "End" = cnv.list$pos[g[,2]], 
                      "Type" = cnv.list$type[g[,1]], "Length" = (g[,2] - g[,1]) * 100, "Intensity" = intensities)
    ## drop element with length=0
    cnv <- cnv[which(cnv[,5] != 0),]
    write.table(cnv, file = paste0(prefix,"_CNV_TV_result.tmp"), append = TRUE, sep = "\t", 
                row.names = FALSE, col.names = FALSE)
    
  }
  close(pb)
}

cnv.filter <- function(){
  d <- read.table(paste0(prefix,"_CNV_TV_result.tmp"), header = T, sep = "\t")
  del <- d[d$Type=="del",]$Intensity
  dup <- d[d$Type=="dup",]$Intensity
  ## del threshold
  if (length(del) > 1) {
    k <- suppressWarnings(kmeans(del,2))
    grp <- which(k$centers == max(k$centers))
    del <- del[which(k$cluster == grp)]
  }
  params = suppressWarnings(fitdistr(del, "poisson"))
  thr.del <- params$estimate + qnorm(0.95)*params$sd/sqrt(length(del))
  
  ## dup threshold
  if (length(dup) > 1) {
    k <- suppressWarnings(kmeans(dup,2))
    grp <- which(k$centers == min(k$centers))
    dup <- dup[which(k$cluster == grp)]
  }
  params = suppressWarnings(fitdistr(dup, "poisson"))
  thr.dup <- params$estimate - qnorm(0.975)*params$sd/sqrt(length(dup))
  
  d <- d[which((d$Type=="del" & d$Intensity < thr.del) | (d$Type=="del" & d$Length >= 800) | (d$Type=="dup" & d$Intensity > thr.dup) | (d$Type=="dup" & d$Length >= 800)),]
  
  write.table(d, file = paste0(prefix,"_CNV_TV_result.txt"), append = FALSE, sep = "\t", 
              row.names = FALSE, col.names = TRUE)
  #file.remove(paste0(prefix,"_CNV_TV_result.tmp"))
}

run.cnv.tv <- function(file, fasta){
  file <- "Data/NA07037.chrom20.ILLUMINA.bwa.CEU.low_coverage.20101123.bam"
  fasta <- "Data/chr20.fa"
  
  file <- "Data/sim3.bam"
  fasta <- "Data/NC_008253.fa"
  
  ref <<- readDNAStringSet(fasta)
  
  cvg <- extractCoverageFromBAM(file)
  
  # Create directory for graphics output
  dir.create(paste0("Plots_", prefix), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  
  write("Chromosome\tStart\tEnd\tType\tLength\tIntensity", file = paste0(prefix,"_CNV_TV_result.tmp"),
        append = FALSE)
  chrom.names <- names(cvg@listData)
  chrom.index <- 1
  for (i in cvg@listData) {
    #meanCvg <<- mean(i)
    extract.cnv(i, chrom.names[chrom.index])
    chrom.index <- chrom.index + 1
  }
  cnv.filter()
  
}


###########################################################################
d <- data.frame(dp=as.numeric(cvg$`20`[60000:2000000]))

ggplot(d, aes(d$dp)) + 
  geom_histogram(aes(y=..density..),binwidth = 1) +
  scale_y_continuous('frequency') +
  geom_vline(xintercept = 9.7, col = "green") 

hist(d$dp, prob = T, xlab = "Read Depth", main = "Histogram of Read Depth")
abline(v = 9.7, col = "green")
m<-mean(d)
std<-sqrt(var(d))
curve(dnorm(x, mean=m, sd=std), add=TRUE, col = "blue")
d <- remove_outliers(d)
d <- d[!is.na(d)]
params <- fitdistr(d[d != 0], "lognormal")
curve(dlnorm(x, meanlog = params$estimate['meanlog'], sdlog = params$estimate['sdlog'], log = FALSE), add=TRUE, col = "red")

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

