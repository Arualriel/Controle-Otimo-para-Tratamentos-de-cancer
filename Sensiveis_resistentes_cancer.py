#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 13:57:24 2019

@author: alunos
"""


#bibliotecas

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#funcao f

global s1,s2

global g1,g2

global alfa,beta

global eta
#####parametros######
s1,s2,g1,g2,alfa,eta=0.192,0.084,0.15,0.02,1./3.,0.02
beta=2./3. - alfa
b,d=5.85,0.00873

def F1(X):
    S=X[0,0]
    R=X[1,0]
    C=X[2,0]
    f1=-s1*S*np.log((R+S)/C)-g1*S+g2*R
    f2=-s2*R*np.log((R+S)/C)+g1*S-g2*R
    f3=b*((R+S)**((2./3.)-alfa))*(C**(1.-beta))-d*((R+S)**(2./3.))*C-eta*C
    saida=np.matrix([[f1],[f2],[f3]])
    return saida



#P 0 = 0.1031 and Q 0 = 0.8869

#condicao inicial

S,R,C= 10.*10**(-6),2.*10**(-6),15.*10**(-6) #metros cúbicos
#S,R,q=1.0,.2,1.5
#X0=np.matrix([[x0],[y0]])
Y0=np.zeros((3,1))
Y0[0,0],Y0[1,0],Y0[2,0]=S,R,C


#variaveis auxiliares

num=70 #numero de iteracoes

Y1=np.zeros((3,1)) #vetor do passo seguinte

dt=1. #taxa de variacao do tempo
t=np.arange(0,num,dt)
Ry=np.zeros((3,num)) #matriz de resultados

Ry[:,0]=Y0[:,0]

#aplicandoo o metodo

for k in range(1,num):
    K1=F1(Y0)
    K2=F1(Y0+dt*K1*0.5)
    K3=F1(Y0+dt*K2*0.5)
    K4=F1(Y0+dt*K3)
    Y1=Y0+(dt/6)*(K1+2*K2+2*K3+K4)
    
    Ry[0,k]=Y1[0,0]
    Ry[1,k]=Y1[1,0]
    Ry[2,k]=Y1[2,0]
    Y0=Y1

######grafico 3d
xyr=Ry[0,:]
yyr=Ry[1,:]
zyr=Ry[2,:]
fig0y=plt.figure()
ax0y=fig0y.add_subplot(111,projection='3d')
ax0y.scatter(xyr,yyr,zyr)
plt.show()
####grafico 2d
#fig1=plt.figure()
#ax1=fig1.add_subplot(3,1,1)
#ax2=fig1.add_subplot(3,1,2)
#ax3=fig1.add_subplot(3,1,3)
linha1,=plt.plot( t, xyr, label="Células sensíveis", color='green',ls='-.')

plt.xlabel("t")
plt.ylabel("S")
plt.legend(handles=[linha1])
plt.show()

linha2,=plt.plot( t, yyr, label="Células resistentes", color='blue',ls='-.')

plt.xlabel("t")
plt.ylabel("R")
plt.legend(handles=[linha2])
plt.show()

linha3,=plt.plot( t, zyr, label="Capacidade de carga", color='red',ls='-.')

plt.xlabel("t")
plt.ylabel("C")
plt.legend(handles=[linha3])

iss=0
maiors=0
for i in range(num):
    if Ry[0,i]>=maiors:
        maiors=Ry[0,i]
        iss=i
        

maiorr=0
ir=0
for i in range (num):
    if Ry[1,i]>=maiorr:
        maiorr=Ry[1,i]
        ir=i

ic=0
maiorc=0
for i in range (num):
    if Ry[2,i]>=maiorc:
        maiorc=Ry[2,i]
        ic=i

print ("pico s",maiors,"indice s",iss)

print ("pico r",maiorr,"indice r",ir)

print ("pico c",maiorc,"indice c",ic)

#ax1.plot(t,xr,c='b')
#ax2.plot(t,yr,c='m')
#ax3.plot(t,zr,c='r')
print(Ry[:,num-1])

plt.show()
#plt.scatter(xr,yr)
