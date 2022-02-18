#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 15:01:20 2019

@author: laura
"""

#bibliotecas

import numpy as np

#funcao f

def F(X):
    x=X[0,0]
    y=X[1,0]
    xi=0.084
    #b=5.85
    #mu=0.02
    f2=-xi*x*np.log(x/y)
    f1=0
    saida=np.matrix([[f1],[f2]])
    return saida
#def f(y):
#    f=-y
#    return f
#condicao inicial

import matplotlib.pyplot as plt


x02,y02=0.4,0.3
x04,y04=0.4,0.3
x0e,y0e=0.4,0.3
#X0=np.matrix([[x0],[y0]])


X02=np.zeros((2,1))
X02[0,0],X02[1,0]=x02,y02


X04=np.zeros((2,1))
X04[0,0],X04[1,0]=x04,y04

X0e=np.zeros((2,1))
X0e[0,0],X0e[1,0]=x0e,y0e   
    
#variaveis auxiliares
   
num=20 #numero de iteracoes
    
X1=np.zeros((2,1)) #vetor do passo seguinte
    
dt=0.1 #taxa de variacao do tempo
    
R2=np.zeros((2,num)) #matriz de resultados
    
R2[:,0]=X02[:,0]
    
R4=np.zeros((2,num)) #matriz de resultados
    
R4[:,0]=X04[:,0]

Re=np.zeros((2,num)) #matriz de resultados
    
Re[:,0]=X0e[:,0]


#aplicandoo o metodo
X=np.zeros(num)
for k in range(1,num):

    K1=F(X04)
    K2=F(X04+dt*K1*0.5)
    K3=F(X04+dt*K2*0.5)
    K4=F(X04+dt*K3)
    
    X12=X04+(dt/6)*(K1+2.0*K2+2.0*K3+1.0*K4)
    
    X14=X02+0.5*dt*(F(X02)+F(X02+dt*F(X02))) 
    
    X1e=X0e+dt*F(X0e)
    
    R2[0,k]=X12[0,0]
    R2[1,k]=X12[1,0]
    X02=X12
    R4[0,k]=X14[0,0]
    R4[1,k]=X14[1,0]
    X04=X14
    Re[0,k]=X1e[0,0]
    Re[1,k]=X1e[1,0]
    X0e=X1e
    X[k]=k

    
xr2=R2[0,:]
yr2=R2[1,:]

xr4=R4[0,:]
yr4=R4[1,:]

xre=Re[0,:]
yre=Re[1,:]


linha1,=plt.plot( X, yr4, label="Runge-Kutta 4", color='orange')

linha2,=plt.plot( X, yr2, label="Runge-Kutta 2", color='blue')

linha3,=plt.plot( X, yre, label="Euler", color='green')

plt.xlabel("x")
plt.ylabel("y")
plt.legend(handles=[linha1, linha2,linha3])


    
plt.plot(yr2)
plt.plot(yr4)
plt.plot(yre)
#plt.scatter(xr,yr)
plt.plot()
