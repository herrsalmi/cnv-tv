#require(ggplot2)
#library(plotly)
library(GenomicAlignments)
library(genlasso)


binSize<-100
windowSize <- 0.2e+6

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

if (get_os() == "windows") {
  NULL.DEV <- "NUL"
} else {
  NULL.DEV <- "/dev/null"
}


slidingwindowplot <- function(windowsize, inputseq) {
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

extractCoverageFromBAM <- function(file) {
  param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE,
                                         isDuplicate=FALSE),
                        what=c("rname", "pos", "cigar"))
  bam <- scanBam(file, param=param)[[1]]
  ## Note that unmapped reads and reads that are PCR/optical duplicates
  ## have already been filtered out by using the ScanBamParam object above.
  irl <- extractAlignmentRangesOnReference(bam$cigar, pos=bam$pos,
                                           f=bam$rname)
  irl <- irl[elementLengths(irl) != 0] # drop empty elements
  coverage(irl)
}

read.data <- function(file){
  #all.data <- extractCoverageFromBAM(file)@listData[[1]] # we deal only with one chromosome
  #myvector_all <- as.vector(all.data)
  myvector_all <- as.vector(file)
  windowAll <- slidingwindowplot(binSize, myvector_all)
  df <- data.frame(windowAll[[1]],windowAll[[2]],windowAll[[3]])
  colname <- c("x","mean","sd")
  colnames(df) <- colname
  return(df)
}

get.cnv <- function(data, start, end, chrom.name){
  y <- data[data$x >= start & data$x <= end,]
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
  svg(paste0("Plots/", chrom.name, "_", start, "_", end, ".svg"),width=14,height=7)
  plot(out, lambda = cv$lambda.1se, pch = ".", cex = 2)
  lines(y$x, x.t, col = "red")
  
  #upper.q <- quantile(x.t, prob = 0.975)
  #lower.q <- quantile(x.t, prob = 0.025)
  #f <- (y$mean[y$mean >= lower.q & y$mean <= upper.q])
  f <- y$mean[y$mean > 0]
  if(length(f) == 0) {
    dev.off()
    return(NULL)
  }
  params <- fitdistr(f, "lognormal")
  
  upper.thr <- qlnorm(0.95, params$estimate['meanlog'], params$estimate['sdlog'])
  lower.thr <- qlnorm(0.05, params$estimate['meanlog'], params$estimate['sdlog'])
  
  abline(a = lower.thr, b = 0, col = "green")
  abline(a = upper.thr, b = 0, col = "green")
  dev.off()
  type <- vector(mode = "character")
  type[x.t >= upper.thr] <- "dup"
  type[x.t <= lower.thr] <- "del"
  x.t[x.t < upper.thr & x.t > lower.thr] <- 0
  return(data.frame(intensity = x.t, pos = y$x, type = type))
}

extract.cnv <- function(chrom, chrom.name) {
  starts <- seq(1, length(chrom)-windowSize, by=windowSize)
  starts <- c(starts, length(chrom))
  n <- length(starts) - 1
  data <- read.data(chrom)
  for (i in 1:n) {
    print(paste0("Processing ", chrom.name, ":", starts[i], "-",starts[i+1]))
    cnv.list <- get.cnv(data, starts[i], starts[i+1], chrom.name)
    if(is.null(cnv.list)) {
      next()
    }
    i <- which(!is.na(cnv.list$type))
    g <- group(i)
    intensities <- vector(mode = "numeric")
    for (index in 1:(length(g)/2)) {
      intensities <- c(intensities, max(cnv.list$intensity[g[index,1]:g[index,2]]))
    }
    cnv <- data.frame("Chrom" = chrom.name, "Start" = cnv.list$pos[g[,1]], "End" = cnv.list$pos[g[,2]], 
                      "Type" = cnv.list$type[g[,1]], "Length" = (g[,2] - g[,1]) * 100, "Intensity" = intensities)
    
    write.table(cnv, file = "CNV_TV_result.txt", append = TRUE, sep = "\t", 
                row.names = FALSE, col.names = FALSE)
    
  }
}

run.cnv.tv <- function(file){
  # file <- "Data/rocker.chr10.bam"
  all.data <- extractCoverageFromBAM(file)
  #windowSize <- 0.3e+6
  
  # Create directory for graphics output
  dir.create("Plots", showWarnings = TRUE, recursive = FALSE, mode = "0777")
  
  write("Chromosome\tStart\tEnd\tType\tLength\tIntensity", file = "CNV_TV_result.txt",
              append = FALSE)
  chrom.names <- names(all.data@listData)
  chrom.index <- 1
  for (i in all.data@listData) {
    extract.cnv(i, chrom.names[chrom.index])
    chrom.index <- chrom.index + 1
  }
  
  ## OLD
  
  #data <- read.data(file)
  #cnv.list <- get.cnv(data, 5.935e+7, 5.965e+7)
  #if(is.null(cnv.list)) {
  #  print("No reads")
  #}
  
}


## clustering ################
## Not Run
dat <- cbind(cnv.list$pos[i], i)

n <- length(i)
n1 <- ceiling(n*2/3)

# percentage of variance explained by clusters
p.exp <- rep(0,n1)

# minimum correlation among all components in each cluster  
min.cor <- matrix(1,n1,n1)  

for (i in 2:n1) {
  fit <- kmeans(dat, centers=i, iter.max=100, nstart=100)
  p.exp[i] <- 1- fit$tot.withinss / fit$totss
}

# minimum number of clusters that explain at least 99.9% of variance
min(which(p.exp > 0.999))

## End(Not Run)



