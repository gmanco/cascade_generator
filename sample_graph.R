
library(igraph)


setwd('/Users/manco/Google Drive/SurvivalFactorization/synth/BA')
perc = 0.05


g1 = barabasi.game(50, 
                           m = 2,
                           power = 1.6, 
                           out.pref = TRUE,
                           zero.appeal = 0.5,
                           directed = TRUE)
V(g1)$name = c(1:50)   


g2 = barabasi.game(50, 
                           m = 1,
                           power = 1.2, 
                           out.pref = TRUE,
                           zero.appeal = 0.5,
                           directed = TRUE)    
V(g2)$name = c(51:100)   


n_of_edges = ecount(g1)+ecount(g2)

n = ceiling(perc*n_of_edges)

add = matrix(0,n,2)

for (i in 1:n){
  src = sample(1:100,1) 
  if (src > 50)
    dst = sample(1:50,1)
  else
    dst = sample(51:100,1)
  add[i,1] = src
  add[i,2] = dst
}

edges <- rbind(as.data.frame(get.edgelist(g1)),as.data.frame(get.edgelist(g2)))

edges <- rbind(add,edges)
names(edges) = c('Src','Dst')

g = graph.data.frame(edges,directed = TRUE)
plot(g)

write.table(file = "synth.BA.K=2.m=[2-1].p=[1.2-1.6].n=[50-50].overlap=0.05.txt",x = edges,sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)

communities = t(rbind(1:100,union(rep(1,50),rep(2,50))))

write.table(file = "synth.BA.K=2.m=[2-1].p=[1.2-1.6].n=[50-50].overlap=0.05.communities",x = communities,sep = "\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
