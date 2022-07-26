rm(list=ls())
library('forecast')
library('smooth')
library('quantreg')
library('seasonal')
library('reshape2')
library('ggplot2')
library('readxl')

#data
source<-as.data.frame(read_excel('Inventory.xlsx',col_names = T))
source$Day<-1:365
ggplot(source, aes(x=Day, y=Demand)) +
  geom_line(size = 0.7) + 
  labs(x = "Day", y =  "Demand") +
  theme(text = element_text(size = 20)) 
ggsave("inventory_ts.pdf", width = 8, height = 5)

source<-source[,-13]

test_length=1

beta<-0.95

mm<-ceiling(round(nrow(source)*(1-beta),2))

#function

r<-10
c<-6
v<-0
g<-1
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
  timeseries<-ts(series[,ncol(series)],frequency = 7)
  residual<-decompose(timeseries)$random
  ranked_residual<-rank(residual)
  rownames(series)<-NULL
  return(series[ranked_residual<=m|ranked_residual>=(length(timeseries)-m+1),])
}

get_cvar<-function(proportion_series,par){
  m<-nrow(proportion_series)
  order<-par[1]+par[2]*proportion_series[,1]+par[3]*proportion_series[,2]+
    par[4]*proportion_series[,3]+par[5]*proportion_series[,4]+
    par[6]*proportion_series[,5]+par[7]*proportion_series[,6]+
    par[8]*proportion_series[,7]+par[9]*proportion_series[,8]+
    par[10]*proportion_series[,9]+par[11]*proportion_series[,10]+
    par[12]*proportion_series[,11]
  cum<-c()
  for (i in 1:m) {
    exceed_loss<-max(c(-profit_function_linear(order[i],proportion_series[i,12])-par[13],0))
    cum<-c(cum,exceed_loss)
  }
  return(par[13]+2*sum(cum)/m)
}

#rolling origin
test_rolling<-function(series,trainsize,rollnum){
  value.sa<-c()
  value.um<-c()
  value.npc<-c()

  record.temp.sa<-c()
  record.temp.um<-c()
  record.temp.npc<-c()

  
  for (j in 1:rollnum) {
    train_series<-as.data.frame(series[j:(trainsize+j-1),])
    test_series<-as.vector(t(series[(trainsize+j),]))
    
    temp.um<-lm(Demand ~ ., data = train_series)
    temp.npc<-optim(par = rep(0,13), fn = get_cvar, proportion_series=
                       get_proportion(train_series,beta),method = 'L-BFGS-B')

    value.temp.sa<-as.numeric(get_sol_residual(train_series[,ncol(train_series)],beta))
    value.temp.um<-as.numeric(temp.um$coefficients[1]+temp.um$coefficients[2]*test_series[1]+
                                temp.um$coefficients[3]*test_series[2]+temp.um$coefficients[4]*test_series[3]+
                                temp.um$coefficients[5]*test_series[4]+temp.um$coefficients[6]*test_series[5]+
                                temp.um$coefficients[7]*test_series[6]+temp.um$coefficients[8]*test_series[7]+
                                temp.um$coefficients[9]*test_series[8]+temp.um$coefficients[10]*test_series[9]+
                                temp.um$coefficients[11]*test_series[10]+temp.um$coefficients[12]*test_series[11]+
                                get_sol_residual(temp.um$residuals,beta))
    value.temp.npc<-as.numeric(temp.npc$par[1]+temp.npc$par[2]*test_series[1]+
                                  temp.npc$par[3]*test_series[2]+temp.npc$par[4]*test_series[3]+
                                  temp.npc$par[5]*test_series[4]+temp.npc$par[6]*test_series[5]+
                                  temp.npc$par[7]*test_series[6]+temp.npc$par[8]*test_series[7]+
                                  temp.npc$par[9]*test_series[8]+temp.npc$par[10]*test_series[9]+
                                  temp.npc$par[11]*test_series[10]+temp.npc$par[12]*test_series[11])
    
    value.sa<-c(value.sa,value.temp.sa)
    value.um<-c(value.um,value.temp.um)
    value.npc<-c(value.npc,value.temp.npc)
    
    record.temp.sa<-c(record.temp.sa,profit_function_linear(value.temp.sa,test_series[12]))
    record.temp.um<-c(record.temp.um,profit_function_linear(value.temp.um,test_series[12]))
    record.temp.npc<-c(record.temp.npc,profit_function_linear(value.temp.npc,test_series[12]))
    
  }
  record.value<-data.frame(SA=value.sa,UM=value.um,NPC=value.npc)
  
  worst.sa<-sort(record.temp.sa)[1:ceiling(round(rollnum*(1-beta),2))]
  worst.um<-sort(record.temp.um)[1:ceiling(round(rollnum*(1-beta),2))]
  worst.npc<-sort(record.temp.npc)[1:ceiling(round(rollnum*(1-beta),2))]
  
  record.worst<-data.frame(SA=worst.sa,UM=worst.um,NPC=worst.npc)
  
  
  return(list(record.worst,record.value))
}

roll<-test_rolling(source,200,165)
##########################################
npc_order<-roll[[2]][,3]
plot_df<-data.frame(Day=201:365,Demand=source$Demand[201:365],Order=npc_order)
longData<-melt(plot_df, id.vars = "Day")
ggplot(longData,aes(Day,value,col=variable))+
  xlab("Day") + ylab("Demand/Order") + labs(color = "Series")+
  scale_color_manual(values=c("black","#E69F00"))+
  geom_line()+
  theme(text = element_text(size = 20)) 
ggsave("inventory_out.pdf", width = 8, height = 5)

profit<-roll[[1]]
longprofit<-melt(profit)
ggplot(longprofit, aes(y=value, x=variable,color=variable)) +
  geom_boxplot() + 
  scale_color_manual(values=c("#999999","darkred","#E69F00"))+
  coord_cartesian(ylim = c(600, 850)) +
  theme(text = element_text(size = 20)) +
  labs(x = "Method", y =  "Realised Profit")  
ggsave("inventory_box_0.95.pdf", width = 6, height = 5)
  
  