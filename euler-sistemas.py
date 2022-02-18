#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 13:56:12 2019

@author: laura
"""
#bibliotecas

import numpy as np

#funcao f

def F(X):
    x=X[0,0]
    y=X[1,0]
    f1=(x*(1-y))
    f2=(y*(x-1))
    saida=np.matrix([[f1],[f2]])
    return saida

#condicao inicial

x0,y0=0.5,0.5
#X0=np.matrix([[x0],[y0]])
X0=np.zeros((2,1))
X0[0,0],X0[1,0]=x0,y0


#variaveis auxiliares

num=1000 #numero de iteracoes

X1=np.zeros((2,1)) #vetor do passo seguinte

dt=0.01 #taxa de variacao do tempo

R=np.zeros((2,num)) #matriz de resultados

R[:,0]=X0[:,0]

#aplicandoo o metodo

for k in range(1,num):
    X1=X0+dt*F(X0)
    
    R[0,k]=X1[0,0]
    R[1,k]=X1[1,0]
    X0=X1

import matplotlib.pyplot as plt

xr=R[0,:]
yr=R[1,:]

plt.plot(xr,yr)
#plt.scatter(xr,yr)

    
    
    