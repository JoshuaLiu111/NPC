rm(list=ls())
library('forecast')
library('smooth')
library('quantreg')
library('seasonal')
library('reshape2')
library('ggplot2')

#dgp

x_1<-sim.ssarima(frequency=12,orders=list(ar=c(1,1)),lags = c(1,12),AR=c(0.3,0.5),obs = 500,constant = 500,mean=0,sd=200)$data
x_2<-sim.ssarima(frequency=12,orders=list(ar=c(1,1),ma=c(1,0)),lags = c(1,12),AR=c(0.2,0.5),MA=c(0.1,0),obs = 500,constant = 600,mean=0,sd=200)$data
x_3<-sim.ssarima(frequency=12,orders=list(ar=c(1,0),ma=c(1,1)),lags = c(1,12),AR=c(0.3,0),MA=c(0.2,0.1),obs = 500,constant = 300,mean=0,sd=50)$data
x_4<-sim.ssarima(frequency=12,orders=list(ma=c(2,2)),lags = c(1,12),MA=c(0.1,0.2,0.1,0.1),obs = 500,constant = 500,mean=0,sd=100)$data
be<-runif(4,0.25,0.75)
demand_series<-ts(be[1]*x_1+be[2]*x_2+be[3]*x_3+be[4]*x_4+rnorm(500,0,100)+
                    rlaplace(500,0,71)+100*rt(500,5),start = c(1,1),frequency = 12)
series<-data.frame(x_1,x_2,x_3,x_4,demand_series)

drawseries<-data.frame(series,year=time(demand_series))
colnames(drawseries)<-c('Feature_1','Feature_2','Feature_3','Feature_4','Demand','Year')

longData<-melt(drawseries, id.vars = "Year")
ggplot(longData,aes(Year,value,col=variable))+
  xlab("Year") + ylab("Value") + labs(color = "Series")+
  geom_line()

test_length=1

beta<-0.9

mm<-ceiling(round(length(demand_series)*(1-beta),2))

#function

profit_function_nonlinear<-function(order,demand){
  u<-rnorm(100,30,5)
  salvage<-mean((rep((order-demand),100)>=u)*u+(rep((order-demand),100)<u)*rep((order-demand),100))
  over<-20*demand-8*order-4*(order-demand)+5*salvage
  under<-12*order-0.01*(demand-order)^2
  profit<-(order>=demand)*over+(order<demand)*under
  return(profit)
}

quan_cvar095<-0.1288
quan_cvar090<-0.1594

bia_quan_cvar095<-0.022
bia_quan_cvar090<-0.04

get_worse_profit<-function(each_order,each_demand,beta){
  temp<-c()
  for (i in 1:length(each_order)) {
    temp<-c(temp,profit_function_nonlinear(each_order[i],each_demand[i]))
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

get_cvar<-function(proportion_series,par){
  m<-nrow(proportion_series)
  order<-par[1]+par[2]*proportion_series[,1]+par[3]*proportion_series[,2]+
    par[4]*proportion_series[,3]+par[5]*proportion_series[,4]
  cum<-c()
  for (i in 1:m) {
    exceed_loss<-max(c(-profit_function_nonlinear(order[i],proportion_series[i,5])-par[6],0))
    cum<-c(cum,exceed_loss)
  }
  return(par[6]+2*sum(cum)/m)
}



#rolling origin
test_rolling<-function(series,trainsize,rollnum){
  service.temp.npc<-c()
  service.temp.sa<-c()
  service.temp.um<-c()
  service.temp.nf<-c()
  service.temp.sqr<-c()
  service.temp.plm<-c()
  
  record.temp.npc<-c()
  record.temp.sa<-c()
  record.temp.um<-c()
  record.temp.nf<-c()
  record.temp.sqr<-c()
  record.temp.plm<-c()
  
  for (j in 1:rollnum) {
    train_series<-series[j:(trainsize+j-1),]
    test_series<-as.vector(t(series[(trainsize+j),]))
    
    temp.npc<-optim(par = rep(0,6), fn = get_cvar, proportion_series=
                      get_proportion(train_series,beta),method = 'L-BFGS-B')
    temp.um<-lm(demand_series ~ ., data = train_series)
    temp.nf<-auto.arima(ts(train_series[,ncol(train_series)]))
    temp.sqr<-rq(demand_series ~ .,data=train_series,tau=bia_quan_cvar090)$coefficients
    temp.plm<-lm(demand_series ~ ., data = get_proportion(train_series,beta))
    
    value.temp.npc<-as.numeric(temp.npc$par[1]+temp.npc$par[2]*test_series[1]+temp.npc$par[3]*
                                 test_series[2]+temp.npc$par[4]*test_series[3]+temp.npc$par[5]*
                                 test_series[4])
    value.temp.sa<-as.numeric(quantile(train_series[,ncol(train_series)],quan_cvar090))
    value.temp.um<-as.numeric(temp.um$coefficients[1]+temp.um$coefficients[2]*test_series[1]+
                                temp.um$coefficients[3]*test_series[2]+temp.um$coefficients[4]*test_series[3]+
                                temp.um$coefficients[5]*test_series[4]+quantile(temp.um$residuals,quan_cvar090))
    value.temp.nf<-as.numeric(forecast(temp.nf,1)$mean+quantile(temp.nf$residuals,quan_cvar090))
    value.temp.sqr<-as.numeric(temp.sqr[1]+temp.sqr[2]*test_series[1]+temp.sqr[3]*
                                 test_series[2]+temp.sqr[4]*test_series[3]+temp.sqr[5]*test_series[4])
    value.temp.plm<-as.numeric(temp.plm$coefficients[1]+temp.plm$coefficients[2]*test_series[1]+
                                 temp.plm$coefficients[3]*test_series[2]+temp.plm$coefficients[4]*test_series[3]+
                                 temp.plm$coefficients[5]*test_series[4]+quantile(temp.plm$residuals,quan_cvar090))
    service.temp.npc<-c(service.temp.npc,value.temp.npc>=test_series[5])
    service.temp.sa<-c(service.temp.sa,value.temp.sa>=test_series[5])
    service.temp.um<-c(service.temp.um,value.temp.um>=test_series[5])
    service.temp.nf<-c(service.temp.nf,value.temp.nf>=test_series[5])
    service.temp.sqr<-c(service.temp.sqr,value.temp.sqr>=test_series[5])
    service.temp.plm<-c(service.temp.plm,value.temp.plm>=test_series[5])
    record.temp.npc<-c(record.temp.npc,profit_function_nonlinear(value.temp.npc,test_series[5]))
    record.temp.sa<-c(record.temp.sa,profit_function_nonlinear(value.temp.sa,test_series[5]))
    record.temp.um<-c(record.temp.um,profit_function_nonlinear(value.temp.um,test_series[5]))
    record.temp.nf<-c(record.temp.nf,profit_function_nonlinear(value.temp.nf,test_series[5]))
    record.temp.sqr<-c(record.temp.sqr,profit_function_nonlinear(value.temp.sqr,test_series[5]))
    record.temp.plm<-c(record.temp.plm,profit_function_nonlinear(value.temp.plm,test_series[5]))
  }
  record.value<-data.frame(NPC=record.temp.npc,SA=record.temp.sa,UM=record.temp.um,
                           NF=record.temp.nf,SQR=record.temp.sqr,PLM=record.temp.plm)
  
  worst.npc<-mean(sort(record.temp.npc)[1:ceiling(round(rollnum*(1-beta),2))])
  worst.sa<-mean(sort(record.temp.sa)[1:ceiling(round(rollnum*(1-beta),2))])
  worst.um<-mean(sort(record.temp.um)[1:ceiling(round(rollnum*(1-beta),2))])
  worst.nf<-mean(sort(record.temp.nf)[1:ceiling(round(rollnum*(1-beta),2))])
  worst.sqr<-mean(sort(record.temp.sqr)[1:ceiling(round(rollnum*(1-beta),2))])
  worst.plm<-mean(sort(record.temp.plm)[1:ceiling(round(rollnum*(1-beta),2))])
  
  service.npc<-mean(service.temp.npc)
  service.sa<-mean(service.temp.sa)
  service.um<-mean(service.temp.um)
  service.nf<-mean(service.temp.nf)
  service.sqr<-mean(service.temp.sqr)
  service.plm<-mean(service.temp.plm)
  
  return(list(c(worst.npc,worst.sa,worst.um,worst.nf,worst.sqr,worst.plm),
              c(service.npc,service.sa,service.um,service.nf,service.sqr,service.plm),
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
     file = 'nonlinear4features010.Rdata')

################################################
load('nonlinear4features005.Rdata')
load('nonlinear4features010.Rdata')


#line_x_200
dl_x_200<-data.frame(NPC=c(roll_50_200[[1]][1],roll_100_200[[1]][1],roll_150_200[[1]][1],
                           roll_200_200[[1]][1],roll_250_200[[1]][1],roll_300_200[[1]][1]),
                     SA=c(roll_50_200[[1]][2],roll_100_200[[1]][2],roll_150_200[[1]][2],
                          roll_200_200[[1]][2],roll_250_200[[1]][2],roll_300_200[[1]][2]),
                     UM=c(roll_50_200[[1]][3],roll_100_200[[1]][3],roll_150_200[[1]][3],
                          roll_200_200[[1]][3],roll_250_200[[1]][3],roll_300_200[[1]][3]),
                     NF=c(roll_50_200[[1]][4],roll_100_200[[1]][4],roll_150_200[[1]][4],
                          roll_200_200[[1]][4],roll_250_200[[1]][4],roll_300_200[[1]][4]),
                     SQR=c(roll_50_200[[1]][5],roll_100_200[[1]][5],roll_150_200[[1]][5],
                           roll_200_200[[1]][5],roll_250_200[[1]][5],roll_300_200[[1]][5]),
                     PLM=c(roll_50_200[[1]][6],roll_100_200[[1]][6],roll_150_200[[1]][6],
                           roll_200_200[[1]][6],roll_250_200[[1]][6],roll_300_200[[1]][6]))


relative_dl_x_200<-data.frame(NPC=(dl_x_200[,1]-dl_x_200[,2])/(dl_x_200[,3]-dl_x_200[,2]),
                              NF=(dl_x_200[,4]-dl_x_200[,2])/(dl_x_200[,3]-dl_x_200[,2]),
                              SQR=(dl_x_200[,5]-dl_x_200[,2])/(dl_x_200[,3]-dl_x_200[,2]),
                              PLM=(dl_x_200[,6]-dl_x_200[,2])/(dl_x_200[,3]-dl_x_200[,2]),
                              origin=c(50,100,150,200,250,300))


draw_dl_x_200<-melt(relative_dl_x_200, id.vars = "origin")
ggplot(draw_dl_x_200,aes(origin,value,col=variable))+
  scale_color_manual(values=c("#E69F00","darkblue", "darkgreen", "violetred"))+
  xlab("Origin") + ylab("Relative 95%-DL") + labs(color = "Method")+
  geom_line(linetype = "dashed")+
  geom_point(shape=4,size=3)+ 
  scale_y_continuous(labels = scales::percent,limits = c(-0.5,1))+
  
ggsave("non_draw_dl_x_200.pdf", width = 8, height = 5)


#line_300_y
dl_300_y<-data.frame(NPC=c(roll_300_50[[1]][1],roll_300_100[[1]][1],roll_300_150[[1]][1],
                           roll_300_200[[1]][1]),
                     SA=c(roll_300_50[[1]][2],roll_300_100[[1]][2],roll_300_150[[1]][2],
                          roll_300_200[[1]][2]),
                     UM=c(roll_300_50[[1]][3],roll_300_100[[1]][3],roll_300_150[[1]][3],
                          roll_300_200[[1]][3]),
                     NF=c(roll_300_50[[1]][4],roll_300_100[[1]][4],roll_300_150[[1]][4],
                          roll_300_200[[1]][4]),
                     SQR=c(roll_300_50[[1]][5],roll_300_100[[1]][5],roll_300_150[[1]][5],
                           roll_300_200[[1]][5]),
                     PLM=c(roll_300_50[[1]][6],roll_300_100[[1]][6],roll_300_150[[1]][6],
                           roll_300_200[[1]][6]))

relative_dl_300_y<-data.frame(NPC=(dl_300_y[,1]-dl_300_y[,2])/(dl_300_y[,3]-dl_300_y[,2]),
                              NF=(dl_300_y[,4]-dl_300_y[,2])/(dl_300_y[,3]-dl_300_y[,2]),
                              SQR=(dl_300_y[,5]-dl_300_y[,2])/(dl_300_y[,3]-dl_300_y[,2]),
                              PLM=(dl_300_y[,6]-dl_300_y[,2])/(dl_300_y[,3]-dl_300_y[,2]),
                              iteration=c(50,100,150,200))

draw_dl_300_y<<-melt(relative_dl_300_y, id.vars = "iteration")
ggplot(draw_dl_300_y,aes(iteration,value,col=variable))+
  scale_y_continuous(labels = scales::percent,limits = c(-0.5,1))+
  scale_color_manual(values=c("#E69F00","darkblue", "darkgreen", "violetred"))+
  xlab("Iteration") + ylab("Relative 95%-DL") + labs(color = "Method")+
  geom_line(linetype = "dashed")+
  geom_point(shape=4,size=3)
ggsave("non_draw_dl_300_y.pdf", width = 8, height = 5)

