require(ggplot2)
library(data.table)
library(genlasso)

binSize<-100

slidingwindowplot <- function(windowsize, inputseq)
{
  starts <- seq(1, length(inputseq)-windowsize, by=windowsize)
  n <- length(starts)
  chunkbps <- numeric(n)
  chunkstats<- numeric(n)
  for (i in 1:n) {
    chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)] 
    chunkmean <- mean(chunk)
    chunkstdv<-sd(chunk)
    chunkbps[i] <- chunkmean
    chunkstats[i]<-chunkstdv
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

all.data <- fread("all.rufus.txt", header = FALSE)

myvector_all<-as.vector(all.data$V3)
windowAll<-slidingwindowplot(binSize,myvector_all)
df<-data.frame(windowAll[[1]],windowAll[[2]],windowAll[[3]])
colname<-c("x","mean","sd")
colnames(df)<-colname

remove(all.data)

y <- df[df$x>5.945e+7 & df$x<5.955e+7,]

y <- df[df$x>5.5e+7 & df$x<6e+7,]

out <- fusedlasso1d(y$mean, y$x)
cv <- cv.trendfilter(out)

x <- coef(out, lambda = cv$lambda.1se)

x <- x$beta
seg <- findSegment(x)
names(seg) <- NULL
ampl <- c()
for (i in 1:length(seg[,1])) {
  ampl <- c(ampl, sum(y$mean[seg[i,1]:seg[i,2]]) / (seg[i,2] - seg[i,1] + 1) )
}
seg.length <- (seg[,2] - seg[,1] + 1)
x_t <- rep(ampl, seg[,2] - seg[,1] + 1)
names(x_t) <- NULL

plot(out, lambda = cv$lambda.1se, pch = ".", cex = 2)
lines(y$x, x_t, col = "red")

abline(a = quantile(y$mean, prob = 0.05), b = 0, col = "green")
abline(a = quantile(y$mean, prob = 0.95), b = 0, col = "green")





##########################################
plot(NULL, ylim = c(-1, 300), xlim = c(min(y$x), max(y$x)))

x.rocker[x.rocker < quantile(y$mean, prob = 0.95)] <- 0


lines(y$x, x.rocker, col = "red")
lines(y$x, x.rufus, col = "blue")
lines(y$x, x.rural, col = "green")
##########################################

pos <- which(x_t > quantile(y$mean, prob = 0.95))
dat <- cbind(y$x[pos], pos)

n = length(pos)
n1 = ceiling(n*2/3)

# percentage of variance explained by clusters
p.exp = rep(0,n1)

# minimum correlation among all components in each cluster  
min.cor = matrix(1,n1,n1)  

for (i in 2:n1) {
  fit = kmeans(dat, centers=i, iter.max=100, nstart=100)
  p.exp[i] = 1- fit$tot.withinss / fit$totss
}

# minimum number of clusters that explain at least 99% of variance
min(which(p.exp > 0.99))


# number of clusters based on elbow method
find.maximum.distance.point(p.exp[-1]) + 1



plot(y$x[pos], (pos)^2)
