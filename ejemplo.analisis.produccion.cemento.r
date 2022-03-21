

#### andres
rm(list=ls())       
graphics.off()

library(lmtest)
install.packages("fpp2")


#------------cargar datos 
require(fpp2)	
beca<-data("qcement")
beca
View(beca)
# Total quarterly production of Portland cement in Australia
# (in millions of tonnes) from 1956:Q1 to 2014:Q1.

y = qcement   
# y = ts(y,frequency=4,start=c(1956,1),end=c(1994,3))
frequency(y)
colnames(y)
#[1] 4

#----- validacion cruzada
yi = window(y,start=c(1956,1),end=c(1992,3))
yf = window(y,start=c(1992,4),end=c(1994,3))
m = 8

# alternativa:
# m = 8
# n = length(y)
# yi = ts(y[1:(n-m)],frequency=4)
# yf = ts(y[(n-m+1):n],frequency=4)


lyi = log(yi)
 

#----modelo con tendencia lineal y estacionalidad
#----con variables indicadoras estacionales
require(forecast)

ti = seq(1,length(yi))

It = seasonaldummy(yi)

mod1 = lm(yi ~ ti + It)
summary(mod1)


r1 = mod1$residuals

# analisis de los residuos  

par(mfrow=c(2,2))
plot(ti,r1,type='o',ylab='residuo')
abline(h=0,lty=2)
plot(density(r1),xlab='x',main= '')
qqnorm(r1)
qqline(r1,col=2)


# valores ajustados

yhat1 = mod1$fitted.values


plot(ti,yi,type='o',col='darkgray')
lines(ti,yhat1,col='red')

#----modelo con tendencia lineal y estacionalidad
#----con variables trigonometricas

It.trig = fourier(yi,2)
mod2 = lm(yi ~ ti + It.trig)
summary(mod2)
r2 = mod1$residuals

M = summary(mod2)
print(xtable(M))

# analisis de los residuos  

par(mfrow=c(2,2))
plot(ti,r2,type='o',ylab='residuo')
abline(h=0,lty=2)
plot(density(r2),xlab='x',main= '')
qqnorm(r2)
qqline(r2,col=2)

yhat2 = mod2$fitted.values
 
B=head(cbind(yhat1,yhat2),n=12)
colnames(B)=c("lineal+indicadoras","lineal+trigonom") 
print(xtable(B))


#---------modelo log cuadratico con estacionalidad
ti2 = ti*ti
mod3 = lm(lyi ~ ti + ti2 + It)
summary(mod3)
 
# El modelo exponencial-cuadratico-estacional
T = length(yi)
Xt = cbind(rep(1,T),ti,ti2,It)
Ds = data.frame(yi,Xt)
theta.0 = mod3$coefficient

mod4 = nls(yi~exp(Xt%*%theta),
data=Ds, start= list(theta=theta.0))

summary(mod4)

yhat4 = fitted(mod4)

plot(ti,yi,type='o',col='darkgray')
lines(ti,yhat1,col='red')
lines(ti,yhat4,col='blue')

#-------- graficas

plot(ti,yi,type='b')
lines(ti,yhat1,col = 'blue')
lines(ti,yhat4,col = 'red')

source("medidas.r")

rbind(medidas(mod1,yi,5),
medidas(mod4,yi,5))

# -------pronosticos 

T = length(yi)
Itf = seasonaldummy(yi,m)
tf = seq(T+1,T+m,1)
 
pron1 = predict(mod1,data.frame(ti = tf, It=I(Itf)))

# pronosticos con el exp cuad 

tf2 = tf*tf
Xtf = cbind(rep(1,m),tf,tf2,Itf)

pron4 = predict(mod4,
data.frame(Xt = I(Xtf)))

par(mfrow=c(1,1))
plot(tf,yf, type = 'o',ylim=c(1.3,2.0))
lines(tf,pron1, type = 'b', pch = 3,col='red' )
lines(tf,pron4, type = 'b', pch = 5,col='blue' )
legend("topleft", 
c("Obs","Indic", "Exp-cuad-estac"), 
pch = c(1, 3, 5),
col = c("black","red","blue"))

#-------medidas de calidad de pronosticos

pron4 = ts(pron4,frequency=4,start=c(1992,4),end=c(1994,3))
pron1 = ts(pron1,frequency=4,start=c(1992,4),end=c(1994,3))


A=rbind(
accuracy(pron1,yf), 
accuracy(pron4,yf))

rownames(A) = c("Lin-Ind","Exp-Cuad-Ind")

library(xtable)
print(xtable(A))


