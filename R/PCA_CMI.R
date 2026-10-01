# PCA_CMI
#
# Input------------------------





#
# Output------------------------



#
# This code is organized by wangchuanyuan1@126.com.

PCA_CMI <- function (data,lamda = 1,order0 = 999,G = NULL)  
{ 
  library(infotheo)  
  n_gene = dim(data)[1]
  if(is.null(G))
  {
    nargin <- 3
    G <- matrix(1, ncol = n_gene, nrow = n_gene)
  }
  G <- G + t(G)
  G[G>1] = 1
  Gval <- G
  order <- 0;
  t <- 0
  # start PCA_CMI
  while(t==0)
  {
    order <- order+1
    if (order0 < 999)
    {
      if (order>order0)
      {
        result1 <- c()
        result1$G <- G    
        result1$Gval <- Gval
        result1$order <- order-1
        
        return(result1)
      }
    }
    result <- edgereduce(G,Gval,order,data,t,lamda)
    G <- result$G
    Gval <- result$Gval
    t <- result$t
    print(t)
    if(t==0)
    {
      print('No edge is reduce! Algorithm  finished!')
      break
    }else
    {t = 0}
  }
  result1 <- c()
  result1$G <- G    
  result1$Gval <- Gval
  result1$order <- order-1
  return(result1)
}  


edgereduce = function (G,Gval,order,data,t,lamda)
{ 
  # edgereduce is pca_cmi

  # lamda: threshold
  dis_data <- discretize(data,disc="equalfreq",(NROW(t(data))^1/2))
  if(order==0)
  {
    for(i in 1:(dim(G)[1]-1))
    {
      for (j in (i+1):dim(G)[1])
      {
        if (G[i,j]!=0)
        {
          mi <- mutinformation(t(dis_data[i,]),t(dis_data[j,]),method="emp")
          Gval[i,j] <- mi 
          Gval[j,i] <- mi
          if (mi<lamda) # set threshold 
          {
            G[i,j] <- 0
            G[j,i] <- 0
          }
        }
      }
    }
    t <- t+1
  }else
  {
    for (i in 1:(dim(G)[1]-1))
    {
      for (j in (i+1):dim(G)[1])
      {
        if (G[i,j]!=0)
        {
          adj <- c()
          for (k in 1:dim(G)[1])
          {
            if (G[i,k]!=0 && G[j,k]!=0)
            {
              adj <- c(adj,k);
            }
          }
          if (length(adj) >= order)
          {
            combntnslist <- t(combn(adj,order));
            combntnsrow <- dim(combntnslist)[1];
            cmiv <- 0
            v1 <- dis_data[i,]
            v2 <- dis_data[j,]
            for (k in (1:combntnsrow))
            {
              vcs <- dis_data[combntnslist[k,],];
              a <- condinformation(t(v1),t(v2),t(vcs),method="emp")
              cmiv <- max(cmiv,a)
            }
            Gval[i,j] <- cmiv
            Gval[j,i] <- cmiv
            if (cmiv<lamda)
            {
              G[i,j] <- 0
              G[j,i] <- 0
            }              
            t <- t+1
          }
        }               
      }
    }             
  }
  result <- c()
  result$G <- G
  result$Gval <- Gval
  result$t <- t
  return(result)
}

tril = function (G, flag = -1) 

{
  if (flag == -1)
  {G <- t(G)}
  m <- dim(G)
  for(i in 1: m[1])
  {
    for(j in i: m[2])
    {G[i,j] <- 0}}
  G <- t(G)
  return(G)
}


