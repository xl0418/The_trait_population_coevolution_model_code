library(ggplot2)
library(ggthemes)
library(viridis)
library(ggridges)
library(tidyverse)
source('C:/Liang/Code/Pro2/R_p2/seprate_vec.R', echo=TRUE)
source('C:/Liang/Code/Pro2/R_p2/theme_henrik.R', echo=TRUE)
source('C:/Liang/Code/Pro2/R_p2/dendf2fredf.R', echo=TRUE)

setwd("C:/Liang/Googlebox/Research/Project2/BaleenWhales/Est")

# test = 4
# file1_name = paste0('smcdataa',test,'.csv')
# adata = read.csv(file1_name)
# adata_matrix = as.matrix(t(adata))
# 
# file2_name = paste0('smcdatag',test,'.csv')
# gammadata = read.csv(file2_name)
# gammadata_matrix = as.matrix(t(gammadata))

file1_name = paste0('est_a','.csv')
adata = read.csv(file1_name)
adata_matrix = as.matrix(t(adata))

file2_name = paste0('est_g','.csv')
gammadata = read.csv(file2_name)
gammadata_matrix = as.matrix(t(gammadata))


gamma_df = as.data.frame(as.table(gammadata_matrix))[,2:3]
gamma_df = cbind(gamma_df,'g')
colnames(gamma_df) = c('Iteration','Samples','para')
levels(gamma_df$Iteration) = c(1:30)

a_df = as.data.frame(as.table(adata_matrix))[,2:3]
a_df = cbind(a_df,'a')

colnames(a_df) = c('Iteration','Samples','para')
levels(a_df$Iteration) = c(1:30)

com_df = rbind(gamma_df,a_df)



# Single mountain plot of the evolution of SMC
ggplot(a_df, aes(x = Samples, y = Iteration)) + 
  geom_density_ridges_gradient(aes(fill = factor(as.integer(a_df$Iteration) %% 2)),
                               scale = 3,size =0.25,alpha = .4, color = "lightblue") +
  theme_ridges()+theme_henrik(grid='', legend.position='none')+
  scale_fill_manual(name = 'Richness',values = c('0' = '#2A7FFF', '1' = '#5599FF'))+
  labs(title = 'Evolution of SMC')+
  theme(legend.position='none')+
  scale_y_discrete(limits = rev(levels(gamma_df$Iteration)))



# Combined mountain plot of the evolution of SMC
ggplot(com_df, aes(x = Samples, y = Iteration)) + 
  geom_density_ridges_gradient(aes(fill = factor(as.integer(com_df$Iteration) %% 2)),
                               scale = 3,size =0.25,alpha = .4, color = "lightblue") +
  theme_ridges()+theme_henrik(grid='', legend.position='none')+
  scale_fill_manual(name = 'Richness',values = c('0' = '#2A7FFF', '1' = '#5599FF'))+
  labs(title = 'Evolution of SMC')+
  theme(legend.position='none')+
  scale_y_discrete(limits = rev(levels(gamma_df$Iteration)))+
  facet_wrap(~para) #+scale_x_continuous(limits = c(-0.5, 1.5))


# Single heatmap of SMC
finvec=c()
for(i in c(1:30)){
  startp = (i-1)*10000+1
  endp = i*10000
  focalvec = a_df$Samples[startp:endp]
  finvec1 = seprate_vec(vec = focalvec,leftendpoint = -0.2,rightendpoint = 0.2,binwidth = 0.001)
  finvec = rbind(finvec, cbind(finvec1,i))
}

htdf = data.frame(finvec)
colnames(htdf) = c('Frequence','value','Iteration')
# htdf$Iteration = as.character(htdf$Iteration)

singlep <- ggplot(htdf, aes(x = value, y = ordered(Iteration, levels =rev(sort(unique(htdf$Iteration)))),
                      fill = Frequence))+ geom_tile() +
  scale_fill_viridis_c(option = "A") +theme_hc()+
  labs(title = "Evolution of SMC",
       x = expression(alpha), y = "Iteration", fill = "Number of samples") +
  theme(legend.position = "bottom", legend.box.just = "bottom")

singlep


# Combined heatmap of SMC
resolution = 0.0001
htdfg = data.frame(dendf2fredf(dendf = gamma_df,iteration = 30,population = 20000,
                               leftendpoint = -.02,rightendpoint = 0.02,binwidth = resolution))
htdfg = cbind(htdfg,'g')
colnames(htdfg) = c('Frequence','value','Iteration','para')
# htdf$Iteration = as.character(htdf$Iteration)


htdfa = data.frame(dendf2fredf(dendf = a_df,iteration = 30,population = 20000,
                               leftendpoint = -.02,rightendpoint = 0.02,binwidth = resolution))
htdfa = cbind(htdfa,'a')
colnames(htdfa) = c('Frequence','value','Iteration','para')
# htdf$Iteration = as.character(htdf$Iteration)

combineddf = rbind(htdfg,htdfa)
levels(combineddf$para) = c('gamma','alpha')

comp <- ggplot(combineddf, aes(x = value, y = ordered(Iteration, levels =rev(sort(unique(htdfg$Iteration)))),
                       fill = Frequence))+
  geom_tile() +
  scale_fill_viridis_c(option = "A") +theme_hc()+
  labs(title = "Evolution of SMC",
       x = "", y = "Iteration", fill = "Number of samples") +
  theme(legend.position = "bottom", legend.box.just = "bottom",strip.background = element_blank())+
  facet_wrap(~para,labeller = label_parsed)

comp



