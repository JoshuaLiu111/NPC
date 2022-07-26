library('ggplot2')
library('r2symbols')
library("ggrepel")

beta<-0.95
x<-rnorm(100000, 0, 3)
quantile<-quantile(x,beta)

df1 <- data.frame(
  x = x
)
ggplot(df1, aes(x)) +
  stat_ecdf(geom = "step")+
  geom_vline(xintercept=quantile, colour="darkred",linetype="dashed")+
  xlim(-8, 8)+
  labs(x = "L", y = expression(Phi))+
  theme(text = element_text(size = 20)) 
ggsave("cdf.pdf", width = 8, height = 5)

df2 <- data.frame(
  x = x[x>quantile]
)
ggplot(df2, aes(x)) +
  stat_ecdf(geom = "step")+
  geom_vline(xintercept=quantile, colour="darkred",linetype="dashed")+
  xlim(-8, 8)+
  labs(x = "L", y = expression(Phi[beta]))+
  theme(text = element_text(size = 20)) 
ggsave("tail.pdf", width = 8, height = 5)


