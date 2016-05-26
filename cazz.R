#!/usr/local/bin/Rscript --vanilla --default-packages=utils

if(!require(optparse)){
  install.packages("optparse")
  library(optparse)
}

option_list = list(
  make_option(c("-n", "--network"), type="character", default=NULL, 
              help="network file name", metavar="character"),
  make_option(c("-c", "--community"), type="character", default=NULL, 
              help="community file name", metavar="character"),
  make_option(c("-u", "--usecommunity"), type="logical", default=FALSE, 
              help="whether to use community file [default= %default]", metavar="logical")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)