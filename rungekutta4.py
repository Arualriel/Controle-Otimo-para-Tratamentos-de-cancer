#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 14:11:25 2019

@author: laura
"""


#bibliotecas

import numpy as np

#funcao f

def F(X):
    x=X[0,0]
    y=X[1,0]
    xi=0.084
    b=5.85
    mu=0.02
    f1 = x*(1.-y)
    f2 = y*(x-1.)
    #f1=-xi*x*np.log(x/y)
    #f2=b*x-(mu-b*x**(2./3.))*y
    saida=np.matrix([[f1],[f2]])
    return saida

#condicao inicial

x0,y0=0.5,0.5
#X0=np.matrix([[x0],[y0]])
X0=np.zeros((2,1))
X0[0,0],X0[1,0]=x0,y0


#variaveis auxiliares

num=100 #numero de iteracoes

X1=np.zeros((2,1)) #vetor do passo seguinte

dt=0.1 #taxa de variacao do tempo

R=np.zeros((2,num-1)) #matriz de resultados

R[:,0]=X0[:,0]
#print(R[:,0],X0[:,0])


#aplicandoo o metodo

for k in range(1,num-1):
    K1=F(X0)
    K2=F(X0+dt*K1*0.5)
    K3=F(X0+dt*K2*0.5)
    K4=F(X0+dt*K3)
    X1=X0+(dt/6)*(K1+2*K2+2*K3+K4)
    
    R[0,k]=X1[0,0]
    R[1,k]=X1[1,0]
    X0=X1

import matplotlib.pyplot as plt

xr=R[0,:]
yr=R[1,:]
t=np.zeros(num)
for i in range(num-1):
    t[i]=dt*i


fig1 = plt.figure()
fig2 = plt.figure()

ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(111)

ax1.plot( t, yr, label="p", color='orange')
ax1.plot( t, xr, label="q", color='blue')


#ax1.xlabel("t")
#ax1.ylabel("x,y")



#ax1.ylim(0,0.10)


numt = len(t)-2
print(num, numt)

print(numt)

print('***',np.shape(t),np.shape(xr))

ax1.plot(t[0:numt],xr[0:numt])
ax1.plot(t[0:numt],yr[0:numt])
ax1.text(t[numt],xr[numt],'aqui')


ax2.plot(xr[:numt],yr[:numt])
#plt.scatter(xr,yr)

plt.show()
    
    
    