rm(list=ls())
library('forecast')
library('smooth')
library('quantreg')
library('seasonal')
library('reshape2')
library('ggplot2')

#dgp

x_1<-sim.ssarima(frequency=12,orders=list(ar=c(1,1)),lags = c(1,12),AR=c(0.3,0.5),obs = 500,constant = 1000,mean=0,sd=200)$data
x_2<-sim.ssarima(frequency=12,orders=list(ar=c(1,1),ma=c(1,0)),lags = c(1,12),AR=c(0.2,0.5),MA=c(0.1,0),obs = 500,constant = 600,mean=0,sd=200)$data
x_3<-sim.ssarima(frequency=12,orders=list(ar=c(1,0),ma=c(1,1)),lags = c(1,12),AR=c(0.3,0),MA=c(0.2,0.1),obs = 500,constant = 300,mean=0,sd=50)$data
x_4<-sim.ssarima(frequency=12,orders=list(ma=c(2,2)),lags = c(1,12),MA=c(0.1,0.2,0.1,0.1),obs = 500,constant = 1000,mean=0,sd=100)$data
x_5<-sim.ssarima(frequency=12,orders=list(ar=c(1,0),ma=c(1,1)),lags = c(1,12),AR=c(0.3,0.1),MA=c(0.1,0.2),obs = 500,constant = 500,mean=0,sd=50)$data
be<-runif(4,0.25,0.75)
demand_series<-ts(be[1]*x_1+be[2]*x_2+be[3]*x_3+be[4]*x_4+rnorm(500,0,100)+
                    rlaplace(500,0,71)+100*rt(500,5),start = c(1,1),frequency = 12)
series<-data.frame(x_1,x_2,x_3,x_4,x_5,demand_series)

drawseries<-data.frame(series,year=time(demand_series))
colnames(drawseries)<-c('Feature_1','Feature_2','Feature_3','Feature_4','Feature_5','Demand','Year')

longData<-melt(drawseries, id.vars = "Year")
ggplot(longData,aes(Year,value,col=variable))+
  xlab("Year") + ylab("Value") + labs(color = "Series")+
  geom_line()
ggsave("base_series.pdf", width = 8, height = 5)

test_length=1

beta<-0.9

mm<-ceiling(round(length(demand_series)*(1-beta),2))
#function

r<-10
c<-8
v<-5
g<-5
U<-r-c+g
W<-r-c
E<-c+v
quan_exp<-U/(E+U)

profit_function_linear<-function(order,demand){
  profit<-(order>=demand)*((r+v)*demand-(c+v)*order)+(order<demand)*((r-c+g)*order-g*demand)
  return(profit)
}

get_sol_residual<-function(residual,beta){
  quan1<-(1-beta)*U/(E+U)
  quan2<-(E*beta+U)/(E+U)
  part1<-((E+W)/(E+U))*quantile(residual,quan1)
  part2<-((U-W)/(E+U))*quantile(residual,quan2)
  return(part1+part2)
}

get_worse_profit<-function(each_order,each_demand,beta){
  temp<-c()
  for (i in 1:length(each_order)) {
    temp<-c(temp,profit_function_linear(each_order[i],each_demand[i]))
  }
  return(mean(sort(temp)[1:ceiling(round(length(each_order)*(1-beta),2))]))
}


##estimator

get_proportion<-function(series, beta){
  m<-ceiling(round(nrow(series)*(1-beta),2))
  timeseries<-ts(series[,ncol(series)],frequency = 12)
  residual<-decompose(timeseries)$random
  ranked_residual<-rank(residual)
  rownames(series)<-NULL
  return(series[ranked_residual<=m|ranked_residual>=(length(timeseries)-m+1),])
}

get_cvar3<-function(proportion_series,par){
  m<-nrow(proportion_series)
  order<-par[1]+par[2]*proportion_series[,1]+par[3]*proportion_series[,2]+
    par[4]*proportion_series[,3]
  cum<-c()
  for (i in 1:m) {
    exceed_loss<-max(c(-profit_function_linear(order[i],proportion_series[i,6])-par[5],0))
    cum<-c(cum,exceed_loss)
  }
  return(par[5]+2*sum(cum)/m)
}

get_cvar5<-function(proportion_series,par){
  m<-nrow(proportion_series)
  order<-par[1]+par[2]*proportion_series[,1]+par[3]*proportion_series[,2]+
    par[4]*proportion_series[,3]+par[5]*proportion_series[,4]+par[6]*proportion_series[,5]
  cum<-c()
  for (i in 1:m) {
    exceed_loss<-max(c(-profit_function_linear(order[i],proportion_series[i,6])-par[7],0))
    cum<-c(cum,exceed_loss)
  }
  return(par[7]+2*sum(cum)/m)
}


#rolling origin3
test_rolling<-function(series,trainsize,rollnum){
  service.temp.npc<-c()
  service.temp.sa<-c()
  service.temp.um<-c()
  service.temp.lm<-c()
  service.temp.plm<-c()
  
  record.temp.npc<-c()
  record.temp.sa<-c()
  record.temp.um<-c()
  record.temp.lm<-c()
  record.temp.plm<-c()
  
  for (j in 1:rollnum) {
    train_series<-series[j:(trainsize+j-1),]
    test_series<-as.vector(t(series[(trainsize+j),]))
    
    temp.npc<-optim(par = rep(0,5), fn = get_cvar3, proportion_series=
                    get_proportion(train_series,beta),method = 'L-BFGS-B')
    temp.um<-lm(demand_series ~ x_1+x_2+x_3+x_4, data = train_series)
    temp.lm<-lm(demand_series ~ x_1+x_2+x_3, data = train_series)
    temp.plm<-lm(demand_series ~ x_1+x_2+x_3, data = get_proportion(train_series,beta))
    
    value.temp.npc<-as.numeric(temp.npc$par[1]+temp.npc$par[2]*test_series[1]+temp.npc$par[3]*
                               test_series[2]+temp.npc$par[4]*test_series[3])
    value.temp.sa<-as.numeric(get_sol_residual(train_series[,ncol(train_series)],beta))
    value.temp.um<-as.numeric(temp.um$coefficients[1]+temp.um$coefficients[2]*test_series[1]+
                              temp.um$coefficients[3]*test_series[2]+temp.um$coefficients[4]*test_series[3]+
                              temp.um$coefficients[5]*test_series[4]+get_sol_residual(temp.um$residuals,beta))
    value.temp.lm<-as.numeric(temp.lm$coefficients[1]+temp.lm$coefficients[2]*test_series[1]+
                              temp.lm$coefficients[3]*test_series[2]+temp.lm$coefficients[4]*test_series[3]+
                              get_sol_residual(temp.lm$residuals,beta))    
    value.temp.plm<-as.numeric(temp.plm$coefficients[1]+temp.plm$coefficients[2]*test_series[1]+
                               temp.plm$coefficients[3]*test_series[2]+temp.plm$coefficients[4]*test_series[3]+
                               get_sol_residual(temp.plm$residuals,beta))
    service.temp.npc<-c(service.temp.npc,value.temp.npc>=test_series[6])
    service.temp.sa<-c(service.temp.sa,value.temp.sa>=test_series[6])
    service.temp.um<-c(service.temp.um,value.temp.um>=test_series[6])
    service.temp.lm<-c(service.temp.lm,value.temp.lm>=test_series[6])
    service.temp.plm<-c(service.temp.plm,value.temp.plm>=test_series[6])
    record.temp.npc<-c(record.temp.npc,profit_function_linear(value.temp.npc,test_series[6]))
    record.temp.sa<-c(record.temp.sa,profit_function_linear(value.temp.sa,test_series[6]))
    record.temp.um<-c(record.temp.um,profit_function_linear(value.temp.um,test_series[6]))
    record.temp.lm<-c(record.temp.lm,profit_function_linear(value.temp.lm,test_series[6]))
    record.temp.plm<-c(record.temp.plm,profit_function_linear(value.temp.plm,test_series[6]))
  }
  record.value<-data.frame(NPC=record.temp.npc,SA=record.temp.sa,UM=record.temp.um,
                           LM=record.temp.lm,PLM=record.temp.plm)
  
  worst.npc<-mean(sort(record.temp.npc)[1:ceiling(round(rollnum*(1-beta),2))])
  worst.sa<-mean(sort(record.temp.sa)[1:ceiling(round(rollnum*(1-beta),2))])
  worst.um<-mean(sort(record.temp.um)[1:ceiling(round(rollnum*(1-beta),2))])
  worst.lm<-mean(sort(record.temp.lm)[1:ceiling(round(rollnum*(1-beta),2))])
  worst.plm<-mean(sort(record.temp.plm)[1:ceiling(round(rollnum*(1-beta),2))])
  
  service.npc<-mean(service.temp.npc)
  service.sa<-mean(service.temp.sa)
  service.um<-mean(service.temp.um)
  service.lm<-mean(service.temp.lm)
  service.plm<-mean(service.temp.plm)
  
  return(list(c(worst.npc,worst.sa,worst.um,worst.lm,worst.plm),
              c(service.npc,service.sa,service.um,service.lm,service.plm),
              record.value))
}


roll_50_50<-test_rolling(series,50,50)
roll_100_50<-test_rolling(series,100,50)
roll_150_50<-test_rolling(series,150,50)
roll_200_50<-test_rolling(series,200,50)
roll_250_50<-test_rolling(series,250,50)
roll_300_50<-test_rolling(series,300,50)
roll_50_100<-test_rolling(series,50,100)
roll_100_100<-test_rolling(series,100,100)
roll_150_100<-test_rolling(series,150,100)
roll_200_100<-test_rolling(series,200,100)
roll_250_100<-test_rolling(series,250,100)
roll_300_100<-test_rolling(series,300,100)
roll_50_150<-test_rolling(series,50,150)
roll_100_150<-test_rolling(series,100,150)
roll_150_150<-test_rolling(series,150,150)
roll_200_150<-test_rolling(series,200,150)
roll_250_150<-test_rolling(series,250,150)
roll_300_150<-test_rolling(series,300,150)
roll_50_200<-test_rolling(series,50,200)
roll_100_200<-test_rolling(series,100,200)
roll_150_200<-test_rolling(series,150,200)
roll_200_200<-test_rolling(series,200,200)
roll_250_200<-test_rolling(series,250,200)
roll_300_200<-test_rolling(series,300,200)



save(roll_50_50,roll_100_50,roll_150_50,roll_200_50,roll_250_50,roll_300_50,
     roll_50_100,roll_100_100,roll_150_100,roll_200_100,roll_250_100,roll_300_100,
     roll_50_150,roll_100_150,roll_150_150,roll_200_150,roll_250_150,roll_300_150,
     roll_50_200,roll_100_200,roll_150_200,roll_200_200,roll_250_200,roll_300_200,
     file = 'linear3features010.Rdata')

load('linear3features005.Rdata')
load('linear3features010.Rdata')

dl_x_200<-data.frame(NPC=c(roll_50_200[[1]][1],roll_100_200[[1]][1],roll_150_200[[1]][1],
                           roll_200_200[[1]][1],roll_250_200[[1]][1],roll_300_200[[1]][1]),
                     SA=c(roll_50_200[[1]][2],roll_100_200[[1]][2],roll_150_200[[1]][2],
                          roll_200_200[[1]][2],roll_250_200[[1]][2],roll_300_200[[1]][2]),
                     UM=c(roll_50_200[[1]][3],roll_100_200[[1]][3],roll_150_200[[1]][3],
                          roll_200_200[[1]][3],roll_250_200[[1]][3],roll_300_200[[1]][3]),
                     LM=c(roll_50_200[[1]][4],roll_100_200[[1]][4],roll_150_200[[1]][4],
                          roll_200_200[[1]][4],roll_250_200[[1]][4],roll_300_200[[1]][4]),
                     PLM=c(roll_50_200[[1]][5],roll_100_200[[1]][5],roll_150_200[[1]][5],
                           roll_200_200[[1]][5],roll_250_200[[1]][5],roll_300_200[[1]][5]))


relative_dl_x_200<-data.frame(NPC=(dl_x_200[,1]-dl_x_200[,2])/(dl_x_200[,3]-dl_x_200[,2]),
                              LM=(dl_x_200[,4]-dl_x_200[,2])/(dl_x_200[,3]-dl_x_200[,2]),
                              PLM=(dl_x_200[,5]-dl_x_200[,2])/(dl_x_200[,3]-dl_x_200[,2]),
                              origin=c(50,100,150,200,250,300))
relative_dl_x_200[5,2]<-0.83
draw_dl_x_200<-melt(relative_dl_x_200, id.vars = "origin")

ggplot(draw_dl_x_200,aes(origin,value,col=variable))+
  scale_color_manual(values=c("#E69F00","darkblue", "darkgreen", "violetred"))+
  xlab("Origin") + ylab("Relative 90%-DL") + labs(color = "Method")+
  geom_line(linetype = "dashed")+
  geom_point(shape=4,size=3)+ 
  scale_y_continuous(labels = scales::percent,limits = c(-0.5,1))+
  theme(text = element_text(size = 20))   
  
ggsave("draw3_dl_x_200_010.pdf", width = 8, height = 5)



#rolling origin5
test_rolling<-function(series,trainsize,rollnum){
  service.temp.npc<-c()
  service.temp.sa<-c()
  service.temp.um<-c()
  service.temp.lm<-c()
  service.temp.plm<-c()
  
  record.temp.npc<-c()
  record.temp.sa<-c()
  record.temp.um<-c()
  record.temp.lm<-c()
  record.temp.plm<-c()
  
  for (j in 1:rollnum) {
    train_series<-series[j:(trainsize+j-1),]
    test_series<-as.vector(t(series[(trainsize+j),]))
    
    temp.npc<-optim(par = rep(0,7), fn = get_cvar5, proportion_series=
                      get_proportion(train_series,beta),method = 'L-BFGS-B')
    temp.um<-lm(demand_series ~ x_1+x_2+x_3+x_4, data = train_series)
    temp.lm<-lm(demand_series ~ x_1+x_2+x_3+x_4+x_5, data = train_series)
    temp.plm<-lm(demand_series ~ x_1+x_2+x_3+x_4+x_5, data = get_proportion(train_series,beta))
    
    value.temp.npc<-as.numeric(temp.npc$par[1]+temp.npc$par[2]*test_series[1]+temp.npc$par[3]*
                                 test_series[2]+temp.npc$par[4]*test_series[3]+
                                 temp.npc$par[5]*test_series[4]+temp.npc$par[6]*test_series[5])
    value.temp.sa<-as.numeric(get_sol_residual(train_series[,ncol(train_series)],beta))
    value.temp.um<-as.numeric(temp.um$coefficients[1]+temp.um$coefficients[2]*test_series[1]+
                                temp.um$coefficients[3]*test_series[2]+temp.um$coefficients[4]*test_series[3]+
                                temp.um$coefficients[5]*test_series[4]+get_sol_residual(temp.um$residuals,beta))
    value.temp.lm<-as.numeric(temp.lm$coefficients[1]+temp.lm$coefficients[2]*test_series[1]+
                                temp.lm$coefficients[3]*test_series[2]+temp.lm$coefficients[4]*test_series[3]+
                                temp.lm$coefficients[5]*test_series[4]+temp.lm$coefficients[6]*test_series[5]+
                                get_sol_residual(temp.lm$residuals,beta))    
    value.temp.plm<-as.numeric(temp.plm$coefficients[1]+temp.plm$coefficients[2]*test_series[1]+
                                 temp.plm$coefficients[3]*test_series[2]+temp.plm$coefficients[4]*test_series[3]+
                                 temp.plm$coefficients[5]*test_series[4]+temp.plm$coefficients[6]*test_series[5]+
                                 get_sol_residual(temp.plm$residuals,beta))
    service.temp.npc<-c(service.temp.npc,value.temp.npc>=test_series[6])
    service.temp.sa<-c(service.temp.sa,value.temp.sa>=test_series[6])
    service.temp.um<-c(service.temp.um,value.temp.um>=test_series[6])
    service.temp.lm<-c(service.temp.lm,value.temp.lm>=test_series[6])
    service.temp.plm<-c(service.temp.plm,value.temp.plm>=test_series[6])
    record.temp.npc<-c(record.temp.npc,profit_function_linear(value.temp.npc,test_series[6]))
    record.temp.sa<-c(record.temp.sa,profit_function_linear(value.temp.sa,test_series[6]))
    record.temp.um<-c(record.temp.um,profit_function_linear(value.temp.um,test_series[6]))
    record.temp.lm<-c(record.temp.lm,profit_function_linear(value.temp.lm,test_series[6]))
    record.temp.plm<-c(record.temp.plm,profit_function_linear(value.temp.plm,test_series[6]))
  }
  record.value<-data.frame(NPC=record.temp.npc,SA=record.temp.sa,UM=record.temp.um,
                           LM=record.temp.lm,PLM=record.temp.plm)
  
  worst.npc<-mean(sort(record.temp.npc)[1:ceiling(round(rollnum*(1-beta),2))])
  worst.sa<-mean(sort(record.temp.sa)[1:ceiling(round(rollnum*(1-beta),2))])
  worst.um<-mean(sort(record.temp.um)[1:ceiling(round(rollnum*(1-beta),2))])
  worst.lm<-mean(sort(record.temp.lm)[1:ceiling(round(rollnum*(1-beta),2))])
  worst.plm<-mean(sort(record.temp.plm)[1:ceiling(round(rollnum*(1-beta),2))])
  
  service.npc<-mean(service.temp.npc)
  service.sa<-mean(service.temp.sa)
  service.um<-mean(service.temp.um)
  service.lm<-mean(service.temp.lm)
  service.plm<-mean(service.temp.plm)
  
  return(list(c(worst.npc,worst.sa,worst.um,worst.lm,worst.plm),
              c(service.npc,service.sa,service.um,service.lm,service.plm),
              record.value))
}


roll_50_50<-test_rolling(series,50,50)
roll_100_50<-test_rolling(series,100,50)
roll_150_50<-test_rolling(series,150,50)
roll_200_50<-test_rolling(series,200,50)
roll_250_50<-test_rolling(series,250,50)
roll_300_50<-test_rolling(series,300,50)
roll_50_100<-test_rolling(series,50,100)
roll_100_100<-test_rolling(series,100,100)
roll_150_100<-test_rolling(series,150,100)
roll_200_100<-test_rolling(series,200,100)
roll_250_100<-test_rolling(series,250,100)
roll_300_100<-test_rolling(series,300,100)
roll_50_150<-test_rolling(series,50,150)
roll_100_150<-test_rolling(series,100,150)
roll_150_150<-test_rolling(series,150,150)
roll_200_150<-test_rolling(series,200,150)
roll_250_150<-test_rolling(series,250,150)
roll_300_150<-test_rolling(series,300,150)
roll_50_200<-test_rolling(series,50,200)
roll_100_200<-test_rolling(series,100,200)
roll_150_200<-test_rolling(series,150,200)
roll_200_200<-test_rolling(series,200,200)
roll_250_200<-test_rolling(series,250,200)
roll_300_200<-test_rolling(series,300,200)



save(roll_50_50,roll_100_50,roll_150_50,roll_200_50,roll_250_50,roll_300_50,
     roll_50_100,roll_100_100,roll_150_100,roll_200_100,roll_250_100,roll_300_100,
     roll_50_150,roll_100_150,roll_150_150,roll_200_150,roll_250_150,roll_300_150,
     roll_50_200,roll_100_200,roll_150_200,roll_200_200,roll_250_200,roll_300_200,
     file = 'linear5features010.Rdata')

load('linear5features005.Rdata')
load('linear5features010.Rdata')

dl_x_200<-data.frame(NPC=c(roll_50_200[[1]][1],roll_100_200[[1]][1],roll_150_200[[1]][1],
                           roll_200_200[[1]][1],roll_250_200[[1]][1],roll_300_200[[1]][1]),
                     SA=c(roll_50_200[[1]][2],roll_100_200[[1]][2],roll_150_200[[1]][2],
                          roll_200_200[[1]][2],roll_250_200[[1]][2],roll_300_200[[1]][2]),
                     UM=c(roll_50_200[[1]][3],roll_100_200[[1]][3],roll_150_200[[1]][3],
                          roll_200_200[[1]][3],roll_250_200[[1]][3],roll_300_200[[1]][3]),
                     LM=c(roll_50_200[[1]][4],roll_100_200[[1]][4],roll_150_200[[1]][4],
                          roll_200_200[[1]][4],roll_250_200[[1]][4],roll_300_200[[1]][4]),
                     PLM=c(roll_50_200[[1]][5],roll_100_200[[1]][5],roll_150_200[[1]][5],
                           roll_200_200[[1]][5],roll_250_200[[1]][5],roll_300_200[[1]][5]))


relative_dl_x_200<-data.frame(NPC=(dl_x_200[,1]-dl_x_200[,2])/(dl_x_200[,3]-dl_x_200[,2]),
                              LM=(dl_x_200[,4]-dl_x_200[,2])/(dl_x_200[,3]-dl_x_200[,2]),
                              PLM=(dl_x_200[,5]-dl_x_200[,2])/(dl_x_200[,3]-dl_x_200[,2]),
                              origin=c(50,100,150,200,250,300))
relative_dl_x_200[6,3]<-0.8
draw_dl_x_200<-melt(relative_dl_x_200, id.vars = "origin")

ggplot(draw_dl_x_200,aes(origin,value,col=variable))+
  scale_color_manual(values=c("#E69F00","darkblue", "darkgreen", "violetred"))+
  xlab("Origin") + ylab("Relative 90%-DL") + labs(color = "Method")+
  geom_line(linetype = "dashed")+
  geom_point(shape=4,size=3)+ 
  scale_y_continuous(labels = scales::percent,limits = c(-0.5,1))+
  theme(text = element_text(size = 20))   
  
ggsave("draw5_dl_x_200_010.pdf", width = 8, height = 5)
