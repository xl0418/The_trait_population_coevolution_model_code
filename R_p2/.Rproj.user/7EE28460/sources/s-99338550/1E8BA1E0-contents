library(ggplot2)
library(ggthemes)
library(viridis)
library(ggridges)
library(tidyverse)
source('C:/Liang/Code/Pro2/R_p2/seprate_vec.R', echo=TRUE)
source('C:/Liang/Code/Pro2/R_p2/theme_henrik.R', echo=TRUE)
source('C:/Liang/Code/Pro2/R_p2/dendf2fredf.R', echo=TRUE)

# setwd("C:/Liang/Googlebox/Research/Project2/modelsele/example1/")
setwd("C:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster/")
# 
# generating = 'TP'
# file1_name = paste0('bestmodel',generating,'.csv')
file1_name = paste0('BWMSt20000_d1_f0.csv')


bmdata = read.csv(file1_name)
bmdata_matrix = as.matrix(t(bmdata))
num.iterations <- dim(bmdata_matrix)[2]
bm_df = as.data.frame(as.table(bmdata_matrix))[,2:3]

colnames(bm_df) = c('Iteration','Samples')
levels(bm_df$Iteration) = c(1:num.iterations)


# 
# 
# # Single mountain plot of the evolution of SMC
# ggplot(bm_df, aes(x = Samples, y = Iteration)) +
#   geom_density_ridges_gradient(aes(fill = factor(as.integer(bm_df$Iteration) %% 2)),
#                                scale = 3,size =0.25,alpha = .4, color = "lightblue") +
#   theme_ridges()+theme_henrik(grid='', legend.position='none')+
#   scale_fill_manual(name = 'Richness',values = c('0' = '#2A7FFF', '1' = '#5599FF'))+
#   labs(title = 'Evolution of SMC')+
#   theme(legend.position='none')+
#   scale_y_discrete(limits = rev(levels(bm_df$Iteration)))
# 


# Single heatmap of SMC
finvec=c()
for(i in c(1:num.iterations)){
  model.count <- c(length(which(bmdata_matrix[,i]==0)),
                   length(which(bmdata_matrix[,i]==1)))
  modelnames <- c(0,1)
  finvec = rbind(finvec,cbind(model.count,modelnames,i))
}

htdf = data.frame(finvec)
colnames(htdf) = c('Frequence','value','Iteration')
# htdf$Iteration = as.character(htdf$Iteration)

singlep <- ggplot(htdf, aes(x = value, y = ordered(Iteration, levels =rev(sort(unique(htdf$Iteration)))),
                            fill = Frequence))+ geom_tile() +
  scale_fill_viridis_c(option = "A") +theme_hc()+
  labs(title = "Top 50% goodness of fit",
       x = "Models", y = "Iteration", fill = "Number of samples") +
  theme(legend.position = "bottom", legend.box.just = "bottom")+
  scale_x_discrete(limit = c(0, 1),labels = c("TP","DR"))

singlep



