rm(list=ls())
library('forecast')
library('smooth')
library('quantreg')
library('seasonal')
library('reshape2')
library('ggplot2')
library('readxl')

#data
source<-as.data.frame(read_excel('Foodbank.xlsx',col_names = T))


ggplot(source, aes(x=Week, y=Visit)) +
  geom_line(size = 0.7) + 
  labs(x = "Week", y =  "Number of Visit") +
  theme(text = element_text(size = 20)) 
ggsave("real_ts.pdf", width = 8, height = 5)

source<-source[,-14]

test_length=1

beta<-0.95

mm<-ceiling(round(nrow(source)*(1-beta),2))

#function

profit_function_nonlinear<-function(order,demand){
  over<-15*(order-demand)
  under<-0.01*(demand-order)^2
  profit<-(order>=demand)*over+(order<demand)*under
  return(profit)
}

quan_cvar095<-0.15
quan_cvar090<-0.2

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
  timeseries<-ts(series[,ncol(series)],frequency = 4)
  residual<-decompose(timeseries)$random
  ranked_residual<-rank(residual)
  rownames(series)<-NULL
  return(series[ranked_residual<=m|ranked_residual>=(length(timeseries)-m+1),])
}

get_cvar0<-function(proportion_series,par){
  m<-nrow(proportion_series)
  order<-par[1]+par[2]*proportion_series[,9]+par[3]*proportion_series[,10]+
    par[4]*proportion_series[,11]+par[5]*proportion_series[,12]
  cum<-c()
  for (i in 1:m) {
    exceed_loss<-max(c(profit_function_nonlinear(order[i],proportion_series[i,13])-par[6],0))
    cum<-c(cum,exceed_loss)
  }
  return(par[6]+2*sum(cum)/m)
}



get_cvar1<-function(proportion_series,par){
  m<-nrow(proportion_series)
  order<-par[1]+par[2]*proportion_series[,5]+par[3]*proportion_series[,6]+
         par[4]*proportion_series[,7]+par[5]*proportion_series[,8]+
         par[6]*proportion_series[,9]+par[7]*proportion_series[,10]+
         par[8]*proportion_series[,11]+par[9]*proportion_series[,12]
  cum<-c()
  for (i in 1:m) {
    exceed_loss<-max(c(profit_function_nonlinear(order[i],proportion_series[i,13])-par[6],0))
    cum<-c(cum,exceed_loss)
  }
  return(par[10]+2*sum(cum)/m)
}



get_cvar2<-function(proportion_series,par){
  m<-nrow(proportion_series)
  order<-par[1]+par[2]*proportion_series[,1]+par[3]*proportion_series[,2]+
         par[4]*proportion_series[,3]+par[5]*proportion_series[,4]+
         par[6]*proportion_series[,5]+par[7]*proportion_series[,6]+
         par[8]*proportion_series[,7]+par[9]*proportion_series[,8]+
         par[10]*proportion_series[,9]+par[11]*proportion_series[,10]+
         par[12]*proportion_series[,11]+par[13]*proportion_series[,12]
  cum<-c()
  for (i in 1:m) {
    exceed_loss<-max(c(profit_function_nonlinear(order[i],proportion_series[i,13])-par[6],0))
    cum<-c(cum,exceed_loss)
  }
  return(par[14]+2*sum(cum)/m)
}



#rolling origin
test_rolling<-function(series,trainsize,rollnum){
  value.sa<-c()
  value.um<-c()
  value.nf<-c()
  value.plm0<-c()
  value.plm1<-c()
  value.plm2<-c()
  value.lm0<-c()
  value.lm1<-c()
  value.lm2<-c()
  value.npc0<-c()
  value.npc1<-c()
  value.npc2<-c()
  
  record.temp.sa<-c()
  record.temp.um<-c()
  record.temp.nf<-c()
  record.temp.plm0<-c()
  record.temp.plm1<-c()
  record.temp.plm2<-c()
  record.temp.lm0<-c()
  record.temp.lm1<-c()
  record.temp.lm2<-c()
  record.temp.npc0<-c()
  record.temp.npc1<-c()
  record.temp.npc2<-c()
  
  
  for (j in 1:rollnum) {
    train_series<-as.data.frame(series[j:(trainsize+j-1),])
    test_series<-as.vector(t(series[(trainsize+j),]))
    
    temp.um<-lm(Visit ~ ., data = train_series)
    temp.nf<-auto.arima(ts(train_series[,13]))
    temp.plm0<-lm(Visit ~ Holiday + Season_1 + Season_2 + Season_3, data = get_proportion(train_series,beta),na.action=na.omit)
    temp.plm1<-lm(Visit ~ Birth + Death + Covid + Crime + 
                    Holiday + Season_1 + Season_2 + Season_3, data = get_proportion(train_series,beta))
    temp.plm2<-lm(Visit ~ ., data = get_proportion(train_series,beta))
    temp.lm0<-lm(Visit ~ Holiday + Season_1 + Season_2 + Season_3, data = train_series)
    temp.lm1<-lm(Visit ~ Birth + Death + Covid + Crime + 
                   Holiday + Season_1 + Season_2 + Season_3, data = train_series)
    temp.lm2<-lm(Visit ~ ., data = train_series)
    temp.npc0<-optim(par = rep(0,6), fn = get_cvar0, proportion_series=
                      get_proportion(train_series,beta),method = 'L-BFGS-B')
    temp.npc1<-optim(par = rep(0,10), fn = get_cvar1, proportion_series=
                       get_proportion(train_series,beta),method = 'L-BFGS-B')
    temp.npc2<-optim(par = rep(0,14), fn = get_cvar2, proportion_series=
                       get_proportion(train_series,beta),method = 'L-BFGS-B')
    
    value.temp.sa<-as.numeric(quantile(train_series[,ncol(train_series)],quan_cvar095))
    value.temp.um<-as.numeric(temp.um$coefficients[1]+temp.um$coefficients[2]*test_series[1]+
                                temp.um$coefficients[3]*test_series[2]+temp.um$coefficients[4]*test_series[3]+
                                temp.um$coefficients[5]*test_series[4]+temp.um$coefficients[6]*test_series[5]+
                                temp.um$coefficients[7]*test_series[6]+temp.um$coefficients[8]*test_series[7]+
                                temp.um$coefficients[9]*test_series[8]+temp.um$coefficients[10]*test_series[9]+
                                temp.um$coefficients[11]*test_series[10]+temp.um$coefficients[12]*test_series[11]+
                                temp.um$coefficients[13]*test_series[12]++quantile(temp.um$residuals,quan_cvar095))
    value.temp.nf<-as.numeric(forecast(temp.nf,1)$mean+quantile(temp.nf$residuals,quan_cvar095))
    value.temp.plm0<-as.numeric(temp.plm0$coefficients[1]+temp.plm0$coefficients[2]*test_series[9]+
                                 temp.plm0$coefficients[3]*test_series[10]+temp.plm0$coefficients[4]*test_series[11]+
                                 temp.plm0$coefficients[5]*test_series[12]+quantile(temp.plm0$residuals,quan_cvar095))
    value.temp.plm1<-as.numeric(temp.plm1$coefficients[1]+temp.plm1$coefficients[2]*test_series[5]+
                                  temp.plm1$coefficients[3]*test_series[6]+temp.plm1$coefficients[4]*test_series[7]+
                                  temp.plm1$coefficients[5]*test_series[8]+temp.plm1$coefficients[6]*test_series[9]+
                                  temp.plm1$coefficients[7]*test_series[10]+temp.plm1$coefficients[8]*test_series[11]+
                                  temp.plm1$coefficients[9]*test_series[12]+quantile(temp.plm1$residuals,quan_cvar095))
    value.temp.plm2<-as.numeric(temp.plm2$coefficients[1]+temp.plm2$coefficients[2]*test_series[1]+
                                  temp.plm2$coefficients[3]*test_series[2]+temp.plm2$coefficients[4]*test_series[3]+
                                  temp.plm2$coefficients[5]*test_series[4]+temp.plm2$coefficients[6]*test_series[5]+
                                  temp.plm2$coefficients[7]*test_series[6]+temp.plm2$coefficients[8]*test_series[7]+
                                  temp.plm2$coefficients[9]*test_series[8]+temp.plm2$coefficients[10]*test_series[9]+
                                  temp.plm2$coefficients[11]*test_series[10]+temp.plm2$coefficients[12]*test_series[11]+
                                  temp.plm2$coefficients[13]*test_series[12]+quantile(temp.plm2$residuals,quan_cvar095))
    value.temp.lm0<-as.numeric(temp.lm0$coefficients[1]+temp.lm0$coefficients[2]*test_series[9]+
                                  temp.lm0$coefficients[3]*test_series[10]+temp.lm0$coefficients[4]*test_series[11]+
                                  temp.lm0$coefficients[5]*test_series[12]+quantile(temp.lm0$residuals,quan_cvar095))
    value.temp.lm1<-as.numeric(temp.lm1$coefficients[1]+temp.lm1$coefficients[2]*test_series[5]+
                                  temp.lm1$coefficients[3]*test_series[6]+temp.lm1$coefficients[4]*test_series[7]+
                                  temp.lm1$coefficients[5]*test_series[8]+temp.lm1$coefficients[6]*test_series[9]+
                                  temp.lm1$coefficients[7]*test_series[10]+temp.lm1$coefficients[8]*test_series[11]+
                                  temp.lm1$coefficients[9]*test_series[12]+quantile(temp.lm1$residuals,quan_cvar095))
    value.temp.lm2<-as.numeric(temp.lm2$coefficients[1]+temp.lm2$coefficients[2]*test_series[1]+
                                  temp.lm2$coefficients[3]*test_series[2]+temp.lm2$coefficients[4]*test_series[3]+
                                  temp.lm2$coefficients[5]*test_series[4]+temp.lm2$coefficients[6]*test_series[5]+
                                  temp.lm2$coefficients[7]*test_series[6]+temp.lm2$coefficients[8]*test_series[7]+
                                  temp.lm2$coefficients[9]*test_series[8]+temp.lm2$coefficients[10]*test_series[9]+
                                  temp.lm2$coefficients[11]*test_series[10]+temp.lm2$coefficients[12]*test_series[11]+
                                  temp.lm2$coefficients[13]*test_series[12]+quantile(temp.lm2$residuals,quan_cvar095))
    value.temp.npc0<-as.numeric(temp.npc0$par[1]+temp.npc0$par[2]*test_series[9]+
                                  temp.npc0$par[3]*test_series[10]+temp.npc0$par[4]*test_series[11]+
                                  temp.npc0$par[5]*test_series[12])
    value.temp.npc1<-as.numeric(temp.npc1$par[1]+temp.npc1$par[2]*test_series[5]+
                                  temp.npc1$par[3]*test_series[6]+temp.npc1$par[4]*test_series[7]+
                                  temp.npc1$par[5]*test_series[8]+temp.npc1$par[6]*test_series[9]+
                                  temp.npc1$par[7]*test_series[10]+temp.npc1$par[8]*test_series[11]+
                                  temp.npc1$par[9]*test_series[12])
    value.temp.npc2<-as.numeric(temp.npc2$par[1]+temp.npc2$par[2]*test_series[1]+
                                  temp.npc2$par[3]*test_series[2]+temp.npc2$par[4]*test_series[3]+
                                  temp.npc2$par[5]*test_series[4]+temp.npc2$par[6]*test_series[5]+
                                  temp.npc2$par[7]*test_series[6]+temp.npc2$par[8]*test_series[7]+
                                  temp.npc2$par[9]*test_series[8]+temp.npc2$par[10]*test_series[9]+
                                  temp.npc2$par[11]*test_series[10]+temp.npc2$par[12]*test_series[11]+
                                  temp.npc2$par[13]*test_series[12])
    
    value.sa<-c(value.sa,value.temp.sa)
    value.um<-c(value.um,value.temp.um)
    value.nf<-c(value.nf,value.temp.nf)
    value.plm0<-c(value.plm0,value.temp.plm0)
    value.plm1<-c(value.plm1,value.temp.plm1)
    value.plm2<-c(value.plm2,value.temp.plm2)
    value.lm0<-c(value.lm0,value.temp.lm0)
    value.lm1<-c(value.lm1,value.temp.lm1)
    value.lm2<-c(value.lm2,value.temp.lm2)
    value.npc0<-c(value.npc0,value.temp.npc0)
    value.npc1<-c(value.npc1,value.temp.npc1)
    value.npc2<-c(value.npc2,value.temp.npc2)
    
    record.temp.sa<-c(record.temp.sa,profit_function_nonlinear(value.temp.sa,test_series[13]))
    record.temp.um<-c(record.temp.um,profit_function_nonlinear(value.temp.um,test_series[13]))
    record.temp.nf<-c(record.temp.nf,profit_function_nonlinear(value.temp.nf,test_series[13]))
    record.temp.plm0<-c(record.temp.plm0,profit_function_nonlinear(value.temp.plm0,test_series[13]))
    record.temp.plm1<-c(record.temp.plm1,profit_function_nonlinear(value.temp.plm1,test_series[13]))
    record.temp.plm2<-c(record.temp.plm2,profit_function_nonlinear(value.temp.plm2,test_series[13]))
    record.temp.lm0<-c(record.temp.lm0,profit_function_nonlinear(value.temp.lm0,test_series[13]))
    record.temp.lm1<-c(record.temp.lm1,profit_function_nonlinear(value.temp.lm1,test_series[13]))
    record.temp.lm2<-c(record.temp.lm2,profit_function_nonlinear(value.temp.lm2,test_series[13]))
    record.temp.npc0<-c(record.temp.npc0,profit_function_nonlinear(value.temp.npc0,test_series[13]))
    record.temp.npc1<-c(record.temp.npc1,profit_function_nonlinear(value.temp.npc1,test_series[13]))
    record.temp.npc2<-c(record.temp.npc2,profit_function_nonlinear(value.temp.npc2,test_series[13]))
    
  }
  record.value<-data.frame(SA=value.sa,UM=value.um,NF=value.nf,
                           PLM0=value.plm0,PLM1=value.plm1,PLM2=value.plm2,
                           LM0=value.lm0,LM1=value.lm1,LM2=value.lm2,
                           NPC0=value.npc0,NPC1=value.npc1,NPC2=value.npc2)
  
  worst.sa<-sort(record.temp.sa)[1:ceiling(round(rollnum*(1-beta),2))]
  worst.um<-sort(record.temp.um)[1:ceiling(round(rollnum*(1-beta),2))]
  worst.nf<-sort(record.temp.nf)[1:ceiling(round(rollnum*(1-beta),2))]
  worst.plm0<-sort(record.temp.plm0)[1:ceiling(round(rollnum*(1-beta),2))]
  worst.plm1<-sort(record.temp.plm1)[1:ceiling(round(rollnum*(1-beta),2))]
  worst.plm2<-sort(record.temp.plm2)[1:ceiling(round(rollnum*(1-beta),2))]
  worst.lm0<-sort(record.temp.lm0)[1:ceiling(round(rollnum*(1-beta),2))]
  worst.lm1<-sort(record.temp.lm1)[1:ceiling(round(rollnum*(1-beta),2))]
  worst.lm2<-sort(record.temp.lm2)[1:ceiling(round(rollnum*(1-beta),2))]
  worst.npc0<-sort(record.temp.npc0)[1:ceiling(round(rollnum*(1-beta),2))]
  worst.npc1<-sort(record.temp.npc1)[1:ceiling(round(rollnum*(1-beta),2))]
  worst.npc2<-sort(record.temp.npc2)[1:ceiling(round(rollnum*(1-beta),2))]
  
  record.worst<-data.frame(SA=worst.sa,UM=worst.um,NF=worst.nf,
                           PLM0=worst.plm0,PLM1=worst.plm1,PLM2=worst.plm2,
                           LM0=worst.lm0,LM1=worst.lm1,LM2=worst.lm2,
                           NPC0=worst.npc0,NPC1=worst.npc1,NPC2=worst.npc2)
  
  
  return(list(record.worst,record.value))
}

roll<-test_rolling(source,60,44)

save(roll, file = 'real_result.Rdata')

###################################


extreme<-c(4,18,23,28,35,39)
real_demand<-source[61:104,13][extreme]

estimate_demand<-roll[[2]][extreme,]
data_demand<-cbind(estimate_demand,real_demand)

cost<-roll[[1]]
cost_vec<-colMeans(cost)
cost_vec[1]-cost_vec[2]

MAE<-matrix(nrow = 6,ncol = 12)
for (i in 1:6) {
  for (j in 1:12) {
    MAE[i,j]<-abs(data_demand[i,j]-data_demand[i,13])
  }
}
MRAE<-matrix(nrow = 6,ncol = 10)
for (i in 1:6) {
  for (j in 1:10) {
    MRAE[i,j]<-MAE[i,(j+2)]/MAE[i,1]
  }
}
MRAE[is.infinite(MRAE)] <- 0
colMeans(MRAE)
