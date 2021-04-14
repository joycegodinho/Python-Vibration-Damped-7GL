# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 16:27:31 2021

@author: Joyce
"""

from __future__ import print_function
import numpy as np
import sympy 

from scipy.linalg import eigh
from matplotlib import pyplot as plt


### Modelagem do sistema 

## Gravidade
g = 9.81

## Dimensões
L0 = 2.385 #[m]
L1 = (3/5)*1.7
L2 = 1.7-L1
L3 = 1.4

## Massas e Inércias
m_roda = 7 # kg
m_carro = 230 #kg

Ix_carro = 47
Iz_carro = 20

## Constantes Elásticas
k2a = 150*9.81/0.01
k3a = 150*9.81/0.01
k4a = 150*9.81/0.01
k5a = 150*9.81/0.01
k2b = k2a/8
k3b = k3a/8
k4b = 1.5*k2b
k5b = 1.5*k3b

## Constantes de Amortecimento
c1 = 0
c2 = 0.5
c3 = 0.5
c4 = 0.5
c5 = 0.5
c6 = 0
c7 = 0

#### Declarando variaveis simbólicas

y1_dot, y2_dot, y3_dot, y4_dot, y5_dot, gama_dot, teta_dot, y1, y2, y3, y4, y5, gama, teta = sympy.symbols('y1_dot, y2_dot, y3_dot, y4_dot, y5_dot, gama_dot, teta_dot, y1, y2, y3, y4, y5, gama, teta')      
var = [y2, y3, y4, y5, y1, gama, teta]
var_dot = [y2_dot, y3_dot, y4_dot, y5_dot, y1_dot, gama_dot, teta_dot]

#### energia cinetica, potencial

Ec = m_carro*(y1_dot**2)/2 + Ix_carro*(teta_dot**2)/2 + Iz_carro*(gama_dot**2)/2 + m_roda*((y2_dot**2)/2 + (y3_dot**2)/2 + (y4_dot**2)/2 + (y5_dot**2)/2) 
Ep = ((y1 + (gama*L1) + teta*(L3/2) - y3)**2)*k3b/2 + ((y1 - (gama*L2) + teta*(L3/2)- y2)**2)*k2b/2  + ((y1 + (gama*L2) - teta*(L3/2) - y4)**2)*k4b/2 + ((y1 - (gama*L2) - teta*(L3/2) - y5)**2)*k5b/2  + (y2**2)*k2a/2 + (y3**2)*k3a/2 + (y4**2)*k4a/2 + (y5**2)*k5a/2 + m_roda*g*(y1 + y2 + y3 + y4 + y5)

## Força dissipativa
Q = [- c2*y2_dot, - c3*y3_dot, - c4*y4_dot, - c5*y5_dot, - c1*y1_dot, - c6*gama_dot, - c7*teta_dot]

##### Lagrangeano 
L = Ec - Ep
L = sympy.simplify(L)

#### Derivadas em relacao aos deslocamentos
#### Derivadas em relacao as velocidades
#### Equecoes do movimento para cada grau de liberdade 
dL_dy = []
dL_dy_dot = []

graus = []

for i in range(len(var)):
    dL_dy.append(0)
    dL_dy_dot.append(0)
    graus.append(0)
    dL_dy[i] = sympy.diff(L, var[i])
    dL_dy_dot[i] = sympy.diff(L, var_dot[i])
    graus[i] = dL_dy_dot[i] + dL_dy[i]
    #print(graus[i],'\n')
    
  
    
### Montando Matriz de Rigidez
K = []

for i in range(len(var)):
    k = []
    for j in range(len(var)):
        k.append(0)
    K.append(k)
    

for i in range(len(var)):
    for j in range(len(var)):
       
       K[i][j] = graus[i].coeff(var[j],1)

for i in range(len(var)):
    for j in range(len(var)):
        if K[i][j] == K[j][i]:
           pass
        else:
           print('Erro: Matriz de rigidez assimétrica')
           break
print('A matriz é simétrica') 

### Montando Matriz de Amortecimento
C = []

## montando matriz generica n x n
for i in range(len(var)):
    c = []
    for j in range(len(var)):
        c.append(0)
    C.append(c)
    
## Preenchendo a matriz genérica
for i in range(len(var)):
    for j in range(len(var)):
       
       C[i][j] = - Q[i].coeff(var_dot[j],1)

for i in range(len(var)):
    for j in range(len(var)):
        if C[i][j] == C[j][i]:
           pass
        else:
           print('Erro: Matriz de Amortecimento assimétrica')
           break
print('A matriz é simétrica') 
print(C)

### Montando Matriz de Massa

M = []

for i in range(len(var)):
    m = []
    for j in range(len(var)):
        m.append(0)
    M.append(m)
    

for i in range(len(var)):
    for j in range(len(var)):
       
       M[i][j] = graus[i].coeff(var_dot[j],1)
#print(M)
for i in range(len(var)):
    for j in range(len(var)):
        if M[i][j] == M[j][i]:
           pass
        else:
           print('Erro: Matriz de rigidez assimétrica')
           break
print('A matriz é simétrica')


##### Analise Modal 

#### Matriz A

A = []
for i in range(2):
    a = []
    for j in range(2):
        a.append(0)
    A.append(a)
    
N = []
for i in range(len(var)):
    n = []
    for j in range(len(var)):
        n.append(0)
    N.append(n)
     
for i in range(2):
    for j in range(2):
        if ( i != j ):
            A[i][j] = M
        else:
            if ( i == j == 0 ):
                A[i][j] = C
            else:
                A[i][j] = N
                
for i in range(2):
    for j in range(2):
        if A[i][j] == A[j][i]:
           pass
        else:
           print('Erro: Matriz A assimétrica')
           break
print('A matriz A é simétrica')

#### Matriz Ac

Ac = []
for i in range(2*len(var)):
    ac = []
    for j in range(2*len(var)):
        ac.append(0)
    Ac.append(ac)
    

for i in range(2):
    count_i = 0
    for j in range(2):
        count_j = 0
        
        for k in range(len(var)):
            if (count_j == 13):
                count_j = 0
            for l in range(len(var)):
           
                Ac[count_i][count_j] = A[i][j][k][l]
                count_j = count_j + 1              
        count_i = count_i + 1
        



        
print (Ac)           
 
 
#### Matriz D

# Mneg = []

# for i in range(len(var)):
#     mneg = []
#     for j in range(len(var)):
#         Mneg.append(0)
#     Mneg.append(mneg)
    
# Neg = []


# for i in range(len(var)):
#     neg = []
#     for j in range(len(var)):
#         Neg.append(-1)
#     Neg.append(neg)
    

# for i in range(len(var)):
#     for j in range(len(var)):
       
#        Mneg[i][j] = int(Neg[i][j]) * int(M[i][j])

D = []
for i in range(2):
    d = []
    for j in range(2):
        d.append(0)
    D.append(d)
      
for i in range(2):
    for j in range(2):
        if ( i != j ):
            D[i][j] = N
        else:
            if ( i == j == 0 ):
                D[i][j] = K
            else:
                D[i][j] = M
                
for i in range(2):
    for j in range(2):
        if D[i][j] == D[j][i]:
           pass
        else:
           print('Erro: Matriz D assimétrica')
           break
print('A matriz D é simétrica')
#print (D)  

##### autovalores e autovetores
# for i in range(len(var)):
#     for j in range(len(var)):
#         K[i][j] = float(K[i][j])
#         M[i][j] = float(M[i][j])
        
# A = np.array(A)
# D = np.array(D)        
# [wn, B] = eigh(D,A)
# wn = np.sqrt(wn**2)
# #print(wn)
# #print(B) 

# freq_n = np.sqrt(((wn)**2)**(1/2))
# #print(freq_n)

# meff = []
# for i in range(len(var)):
#     meff.append(0)
#     meff[i] = np.dot(np.dot(B[:,i].T,A),B[:,i])
    
# meff = np.diag(meff)
# for i in range(len(var)):
#     for j in range(len(var)):
#         if meff[i][j] == meff[j][i]:
#            pass
#         else:
#            print('Erro: Matriz meff assimétrica')
#            break
# print('A matriz meff é simétrica')

# #### Condições iniciais
# x0 = [0,0.3,0,0,0,0,0]
# v0 = [0,0,0,0,0,0,0]
# y0 = np.dot(np.dot(B.T,A),x0)
# y0dot = np.dot(np.dot(B.T,A),v0)

# t_inc0 = 0.005 #segundos
# t_end0 = 10 #segundos    
# y = []
# t = np.arange(0, t_end0, t_inc0)
# for r in range(len(var)):
#     y.append(0)
#     y[r] = y0[r]*np.cos(freq_n[r]*t) + y0dot[r]/(freq_n[r]*np.sin(freq_n[r]*t))
    
# x = np.dot(B,y)
# x1 = x[0][:]
# x2 = x[1][:]
# x3 = x[2][:]
# x4 = x[3][:]
# x5 = x[4][:]
# x6 = x[5][:]
# x7 = x[6][:]
# #for i in range(len(var)):

# plt.plot(t,x[5][:])

# np.savetxt('7GDL_RespostaAnimBlue', np.transpose([x1,x2,x3,x4,x5,x6,x7]), fmt='%1.5f')


