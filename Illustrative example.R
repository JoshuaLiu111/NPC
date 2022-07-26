library('ggplot2')

r<-20
c<-8
v<-3
g<-7
U<-r-c+g
W<-r-c
E<-c+v
quan_exp<-U/(E+U)

profit_function_linear<-function(order,demand){
  profit<-(order>=demand)*((r+v)*demand-(c+v)*order)+(order<demand)*((r-c+g)*order-g*demand)
  return(profit-100)
}



beta<-0.95

demand<-rnorm(1000,10,2)

df1<-data.frame(demand=demand, order1=-profit_function_linear(5,demand),
                order2=-profit_function_linear(8,demand))

df1_order1<- df1[order(df1$order1,decreasing = T)[1:(1000*0.05)],]
df1_order2<- df1[order(df1$order2,decreasing = T)[1:(1000*0.05)],]
df1_order21<-df1_order2[df1_order2$demand<10,]
df1_order22<-df1_order2[df1_order2$demand>10,]

ggplot(df1,aes(demand,order1))+
  geom_line()+
  geom_area(data=subset(df1,is.element(demand, df1_order1$demand)),fill="#003300", alpha=1.)+
  geom_vline(xintercept=quantile(demand,beta), colour="darkred",linetype="dashed")+
  geom_vline(xintercept=quantile(demand,(1-beta)), colour="darkred",linetype="dashed")+
  labs(x = "Demand", y = "L")+
  theme(text = element_text(size = 20))
ggsave("linear_order1.pdf", width = 8, height = 5)

ggplot(df1,aes(demand,order2))+
  geom_line()+
  geom_area(data=subset(df1,is.element(demand, df1_order21$demand)),fill="#003300", alpha=1.)+
  geom_area(data=subset(df1,is.element(demand, df1_order22$demand)),fill="#003300", alpha=1.)+
  geom_vline(xintercept=quantile(demand,beta), colour="darkred",linetype="dashed")+
  geom_vline(xintercept=quantile(demand,(1-beta)), colour="darkred",linetype="dashed")+
  labs(x = "Demand", y = "L")+
  theme(text = element_text(size = 20))
ggsave("linear_order2.pdf", width = 8, height = 5)


profit_function_nonlinear<-function(order,demand){
  over<-20*demand-8*order-4*(order-demand)
  under<-12*order-0.01*(demand-order)^2
  profit<-(order>=demand)*over+(order<demand)*under
  return(profit-100)
}

df2<-data.frame(demand=demand, order1=-profit_function_nonlinear(5,demand),
                order2=-profit_function_nonlinear(8,demand))

df2_order1<- df2[order(df2$order1,decreasing = T)[1:(1000*0.05)],]
df2_order11<-df2_order1[df2_order1$demand<10,]
df2_order12<-df2_order1[df2_order1$demand>10,]
df2_order2<- df2[order(df2$order2,decreasing = T)[1:(1000*0.05)],]
df2_order21<-df2_order2[df2_order2$demand<10,]
df2_order22<-df2_order2[df2_order2$demand>10,]

ggplot(df2,aes(demand,order1))+
  geom_line()+
  geom_area(data=subset(df2,is.element(demand, df2_order11$demand)),fill="#003300", alpha=1.)+
  geom_area(data=subset(df2,is.element(demand, df2_order12$demand)),fill="#003300", alpha=1.)+
  geom_vline(xintercept=quantile(demand,beta), colour="darkred",linetype="dashed")+
  geom_vline(xintercept=quantile(demand,(1-beta)), colour="darkred",linetype="dashed")+
  labs(x = "Demand", y = "L")+
  theme(text = element_text(size = 20))
ggsave("nonlinear_order1.pdf", width = 8, height = 5)

ggplot(df2,aes(demand,order2))+
  geom_line()+
  geom_area(data=subset(df2,is.element(demand, df2_order21$demand)),fill="#003300", alpha=1.)+
  geom_area(data=subset(df2,is.element(demand, df2_order22$demand)),fill="#003300", alpha=1.)+
  geom_vline(xintercept=quantile(demand,beta), colour="darkred",linetype="dashed")+
  geom_vline(xintercept=quantile(demand,(1-beta)), colour="darkred",linetype="dashed")+
  labs(x = "Demand", y = "L")+
  theme(text = element_text(size = 20))
ggsave("nonlinear_order2.pdf", width = 8, height = 5)
