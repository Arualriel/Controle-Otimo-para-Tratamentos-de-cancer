#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 01:41:28 2020

@author: laura
"""

import numpy as np
import matplotlib.pyplot as plt

global p0,q0,r0,s,m,gama,alfa,beta,delta,xi,b,d,mi
#####parametros######
p0,q0,r0,s,m,gama=.20,.30,0.0001,10.,1.,0.01 ##s=10. ou 0.01
alfa,beta,delta,teta=0.0529,0.00291,0.3743,1.
xi,b,d,mi=0.0347,5,0.0667,0.0


def F(X):
    p=X[0,0]
    q=X[1,0]
    r=X[2,0]
    f1=-xi*p*(np.log(p/q))-teta*p*r
    f2=b*p-(mi+d*p**(2./3.))*q
    f3=alfa*(p-beta*p**2)*r+gama-delta*r
    saida=np.matrix([[f1],[f2],[f3]])
    return saida

X0=np.zeros((3,1))
X0[0,0],X0[1,0],X0[2,0]=p0,q0,r0


num=100 #numero de iteracoes

X1=np.zeros((3,1)) #vetor do passo seguinte

dt=1. #taxa de variacao do tempo
t=np.arange(0,num,dt)
Rx=np.zeros((3,num)) #matriz de resultados
Rx[:,0]=X0[:,0]




for k in range(1,num):
    K1=F(X0)
    K2=F(X0+dt*K1*0.5)
    K3=F(X0+dt*K2*0.5)
    K4=F(X0+dt*K3)
    X1=X0+(dt/6)*(K1+2*K2+2*K3+K4)
    
    Rx[0,k]=X1[0,0]
    Rx[1,k]=X1[1,0]
    Rx[2,k]=X1[2,0]
    X0=X1

######grafico 3d
xyr=Rx[0,:]
yyr=Rx[1,:]
zyr=Rx[2,:]
fig0y=plt.figure()
ax0y=fig0y.add_subplot(111,projection='3d')
ax0y.scatter(xyr,yyr,zyr)
plt.show()
####grafico 2d
#fig1=plt.figure()
#ax1=fig1.add_subplot(3,1,1)
#ax2=fig1.add_subplot(3,1,2)
#ax3=fig1.add_subplot(3,1,3)
linha1,=plt.plot( t, xyr, label="Volume tumoral", color='red',ls='-.')

plt.xlabel("t")
plt.ylabel("P")
plt.legend(handles=[linha1])
plt.show()

linha2,=plt.plot( t, yyr, label="Capacidade de carga tumoral", color='navy',ls='-.')

plt.xlabel("t")
plt.ylabel("Q")
plt.legend(handles=[linha2])
plt.show()

linha3,=plt.plot( t, zyr, label="Densidade de células imunocompetentes", color='green',ls='-.')

plt.xlabel("t")
plt.ylabel("R")
plt.legend(handles=[linha3])

iss=0
maiors=0
for i in range(num):
    if Rx[0,i]>=maiors:
        maiors=Rx[0,i]
        iss=i
        

maiorr=0
ir=0
for i in range (num):
    if Rx[1,i]>=maiorr:
        maiorr=Rx[1,i]
        ir=i

ic=0
maiorc=0
for i in range (num):
    if Rx[2,i]>=maiorc:
        maiorc=Rx[2,i]
        ic=i

print ("pico p",maiors,"indice p",iss)

print ("pico q",maiorr,"indice q",ir)

print ("pico r",maiorc,"indice r",ic)

#ax1.plot(t,xr,c='b')
#ax2.plot(t,yr,c='m')
#ax3.plot(t,zr,c='r')
print(Rx[:,num-1])

plt.show()
#plt.scatter(xr,yr)


