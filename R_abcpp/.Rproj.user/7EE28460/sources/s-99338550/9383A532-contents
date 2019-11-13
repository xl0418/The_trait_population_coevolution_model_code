library(ggplot2)

comp.1 = function(x)y = 0.1*x*exp(-0.1*x^2)
comp.2 = function(x)y = 0.5*x*exp(-0.5*x^2)
comp.3 = function(x)y = 1*x*exp(-1*x^2)
comp.4 = function(x)y = 3*x*exp(-1*x^2)+0.5*x*exp(-0.5*x^2)+0.1*x*exp(-0.1*x^2)



ggplot(data.frame(x=c(-5,5)), aes(x)) +
  stat_function(fun=comp.1,size = 1.5, geom="line", aes(colour="comp.1")) +
  stat_function(fun=comp.2,size = 1.5, geom="line", aes(colour="comp.2")) +
  stat_function(fun=comp.3,size = 1.5, geom="line", aes(colour="comp.3")) +
  stat_function(fun=comp.4,size = 1.5, geom="line", aes(colour="comp.4")) +
  scale_colour_manual("Function", values=c("blue","red","green",'purple'), 
                      breaks=c("comp.1","comp.2","comp.3","comp.4"),
                      labels = c("a = 0.1", "a = 0.5","a = 1",'sum'))+
  xlab("Distance")+ylab("Competitive strength")+
  theme_bw()
