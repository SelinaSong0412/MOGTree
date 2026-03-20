library(data.tree)
library(dplyr)

ten2bin<-function(x,max_h){
  tmp<- x %%   2^((max_h) :1)
  tmp<-tmp %/% 2^((max_h):1 - 1)
  return(  tmp  )
}

bin2ten<-function(vec){
  n<-length(vec)
  sum( vec * 2^(n:1 - 1)  )
}

vis_tree = function(max_h,treeout,number){
  dat<-matrix(0,2^max_h,max_h)
  treeout$cutoff = round(treeout$cutoff, 3)
  for (k in 1L:(max_h)) { # depth is the most number of split to reach one terminal node
    # The 1st split
    if (k==1L) {
      tmp = rep( c(paste("X",treeout$X[k],"<=",treeout$cutoff[k],sep=""),
                   paste("X",treeout$X[k],">",treeout$cutoff[k],sep="")),
                 each = 2^(max_h-k))
      dat[,1] = tmp
    } else {
      res = vector()
      for(j in (2^(k-1)):(2^k-1)) {
        if(!is.na(treeout$X[j])){
          tmp = rep( c(paste("X",treeout$X[j],"<=",treeout$cutoff[j],sep=""),
                       paste("X",treeout$X[j],">",treeout$cutoff[j],sep="")),
                     each = 2^(max_h-k))
        }else{
          tmp = rep( c(NA,NA),
                     each = 2^(max_h-k))
        }
        res = c(res,tmp)
      }
      dat[,k] = res
    }
  }
  dat = unique(dat)
  pathString <- paste("All Data", dat[,1],sep = "/")
  for(nn in 1:nrow(dat)){
    for(ii in 2:max_h){
      if(!is.na(dat[nn,ii]))
        pathString[nn] <- paste(pathString[nn], dat[nn,ii],sep = "/")
    }  
  }
  
  pathString <- paste(pathString, number,sep = "/")
  dat = as.data.frame(dat)
  dat$pathString <- pathString
  tree0 <- as.Node(dat)  
  return(tree0)  
}


