library(ggplot2)

r<-10
c<-5
v<-5
g<-5
U<-r-c+g
W<-r-c
E<-c+v
nvp<-function(order,demand){
  profit<-(order>=demand)*((r+v)*demand-(c+v)*order)+(order<demand)*((r-c+g)*order-g*demand)
  return(profit)
}

demand<-sort(rnorm(10000,2000,150))

Q<-seq(1601,2400,20)

Profit<-matrix(nrow = 10000, ncol = 40)

for (i in 1:10000) {
  for (j in 1:40){
    Profit[i,j]<-nvp(Q[j],demand[i])
  }
}


Mean_Profit<-colMeans(Profit)
Sd_Profit<-apply(Profit,sd,MARGIN=2)

lowerquantile<-function(x){
  return(quantile(x,probs = 0.1))
}
Lower_Profit<-apply(Profit,lowerquantile,MARGIN=2)

upperquantile<-function(x){
  return(quantile(x,probs = 0.9))
}
Upper_Profit<-apply(Profit,upperquantile,MARGIN=2)

qProfit<-data.frame(Order=rep(Q,3), Type=rep(c('Mean','Lower','Upper'),each=40), 
                 Profit=c(Mean_Profit,Lower_Profit,Upper_Profit))

ggplot(data=qProfit, aes(x=Order, y=Profit, group=Type)) +
  geom_line(aes(linetype=Type))+
  geom_point()

sProfit<-data.frame(Order=Q, Sd=Sd_Profit)
ggplot(data=sProfit, aes(x=Order,y=Sd_Profit))+
  geom_point()

cVaR<-function(x){
  quantil<-quantile(x,probs = 0.05)
  return(sum(x*(x<quantil))/sum(x<quantil))
}
cVaR_Profit<-apply(Profit,cVaR,MARGIN=2)


Risk_Profit<-data.frame(Order=rep(Q,3), Type=rep(c('Mean','Risk Averse','cVaR'),each=40), 
                        Profit=c(Mean_Profit,0.7*Mean_Profit+0.3*cVaR_Profit,cVaR_Profit))

ggplot(data=Risk_Profit, aes(x=Order, y=Profit, group=Type)) +
  geom_line(aes(linetype=Type),size=1)+
  geom_point(aes(x = Q[which.max(cVaR_Profit)],y = max(cVaR_Profit)),size=2,color = 'darkred')+
  geom_segment(aes(x = Q[which.max(cVaR_Profit)], y = -600, xend = Q[which.max(cVaR_Profit)], 
                   yend = max(cVaR_Profit)),linetype="dashed", color = 'darkred')+
  geom_point(aes(x = Q[which.max(0.7*Mean_Profit+0.3*cVaR_Profit)],
                 y = max(0.7*Mean_Profit+0.3*cVaR_Profit)),size=2,color = 'darkblue')+
  geom_segment(aes(x = Q[which.max(0.7*Mean_Profit+0.3*cVaR_Profit)], y = -600, 
                   xend = Q[which.max(0.7*Mean_Profit+0.3*cVaR_Profit)], 
                   yend = max(0.7*Mean_Profit+0.3*cVaR_Profit)),linetype="dashed", color = 'darkblue')+
  geom_point(aes(x = Q[which.max(Mean_Profit)],y = max(Mean_Profit)),size=2)+
  geom_segment(aes(x = Q[which.max(Mean_Profit)], y = -600, xend = Q[which.max(Mean_Profit)], 
                   yend = max(Mean_Profit)),linetype="dashed")+
  theme(text = element_text(size = 20))

ggsave("mean-cvar.pdf", width = 8, height = 4)
########################

range_c_u<-seq(0.5,9.5,0.25)
range_c_o<-rev(seq(0.5,9.5,0.25))
quan<-range_c_u/(range_c_u+range_c_o)

newnvp<-function(d,x,c_u,c_o){
  P<-1000-(d>=x)*c_u*(d-x)-(d<x)*c_o*(x-d)
  return(P)
}

Opt_record<-matrix(nrow = 37, ncol = 3)

for (k in 1:37){
  Profit<-matrix(nrow = 10000, ncol = 20)
  for (i in 1:10000) {
    for (j in 1:20){
      Profit[i,j]<-newnvp(demand[i],Q[j],range_c_u[k],range_c_o[k])
    }
  }
  Mean_Profit<-colMeans(Profit)
  Sd_Profit<-apply(Profit,sd,MARGIN=2)
  cVaR_Profit<-apply(Profit,cVaR,MARGIN=2)
  Opt_record[k,1]<-Q[which.max(Mean_Profit)]
  Opt_record[k,2]<-Q[which.max(Mean_Profit-0.8*Sd_Profit)]
  Opt_record[k,3]<-Q[which.max(cVaR_Profit)]
}

Tau_Profit<-data.frame(Tau=rep(quan,3), Type=rep(c('Mean','Risk Averse','cVaR'),each=37), 
                        Order=c(Opt_record[,1],Opt_record[,2],Opt_record[,3]))

ggplot(data=Tau_Profit, aes(x=Tau, y=Order, group=Type)) +
  geom_line(aes(linetype=Type))+
  geom_point()
