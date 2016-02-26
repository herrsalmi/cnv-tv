## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      CNV-TV
      
      Arguments:
      --arg1=someValue   - numeric, blah blah
      --arg2=someValue   - character, blah blah
      --arg3=someValue   - logical, blah blah
      --help             - print this text
      
      Example:
      ./test.R --arg1=1 --arg2=text --arg3=TRUE \n\n")
  
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

## Arg1 default
if(is.null(args$arg1)) {
  ## do something
}

## Arg2 default
if(is.null(args$arg2)) {
  ## do something
}

## Arg3 default
if(is.null(args$arg3)) {
  ## do something
}


check.arg <- function(bamfl) {
  if(!file.exists(bamfl)){
    stop(paste("The file:",path(obj),"does not exist."))
  }
  
  if(length(seqinfo(bamfl))==0){
    stop("Your BAM file has no sequence information in its header.")
  }
}