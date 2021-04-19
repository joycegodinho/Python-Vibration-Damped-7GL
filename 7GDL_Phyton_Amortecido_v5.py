# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 16:27:31 2021

@author: Joyce & Lucas
"""

from __future__ import print_function
import numpy as np
import sympy 

from scipy.linalg import eigh
from scipy.sparse.linalg import lobpcg
from matplotlib import pyplot as plt


### Modelagem do sistema 

## Gravidade
g = 9.81

## Dimensões
L0 = 2.385 #[m] #entre eixos
L1 = (3/5)*1.7 #dist roda frente ao CG
L2 = 1.7-L1 #dist roda traseira ao CG
L3 = 1.4 #Bitola

## Massas e Inércias
m_roda = 7 # kg
m_carro = 230 #kg

Ix_carro = 47 #valor qualqr
Iz_carro = 20

## Constantes Elásticas Estimadas
k2a = 150*9.81/0.01
k3a = 150*9.81/0.01
k4a = 150*9.81/0.01
k5a = 150*9.81/0.01
k2b = k2a/8
k3b = k3a/8
k4b = 1.5*k2b
k5b = 1.5*k3b

## Constantes de Amortecimento (chutes)
c1 = 0
c2 = 200
c3 = 200
c4 = 200
c5 = 200
c6 = 0
c7 = 0
c_ = [c1, c2, c3, c4, c5, c6, c7]
#### Declarando variaveis simbólicas

y1, y2, y3, y4, y5, gama, teta, y1_dot, y2_dot, y3_dot, y4_dot, y5_dot, gama_dot, teta_dot, y1_dot2, y2_dot2, y3_dot2, y4_dot2, y5_dot2, gama_dot2, teta_dot2, t_ = sympy.symbols('y1, y2, y3, y4, y5, gama, teta, y1_dot, y2_dot, y3_dot, y4_dot, y5_dot, gama_dot, teta_dot, y1_dot2, y2_dot2, y3_dot2, y4_dot2, y5_dot2, gama_dot2, teta_dot2, t_')      
var_dot2 = [y2_dot2, y3_dot2, y4_dot2, y5_dot2, y1_dot2, gama_dot2, teta_dot2] #Variável de "velocidade" para derivada da energia cinética
var = [y2, y3, y4, y5, y1, gama, teta] #Variável de deslocamento para derivada da energia potencial
var_dot = [y2_dot, y3_dot, y4_dot, y5_dot, y1_dot, gama_dot, teta_dot] #Variável de velocidade pra derivada da energia dissipada em amortecimento viscoso

#### energia cinetica, potencial e dissipativa

Ec = m_carro*(y1_dot2**2)/2 + Ix_carro*(teta_dot2**2)/2 + Iz_carro*(gama_dot2**2)/2 + m_roda*((y2_dot2**2)/2 + (y3_dot2**2)/2 + (y4_dot2**2)/2 + (y5_dot2**2)/2) 
Ep = ((y1 + (gama*L1) + teta*(L3/2) - y3)**2)*k3b/2 + ((y1 - (gama*L2) + teta*(L3/2)- y2)**2)*k2b/2  + ((y1 + (gama*L2) - teta*(L3/2) - y4)**2)*k4b/2 + ((y1 - (gama*L2) - teta*(L3/2) - y5)**2)*k5b/2  + (y2**2)*k2a/2 + (y3**2)*k3a/2 + (y4**2)*k4a/2 + (y5**2)*k5a/2 + m_roda*g*(y1 + y2 + y3 + y4 + y5)
Ed = c2*(1/2)*(y1_dot - y2_dot - teta_dot*L3/2 + gama_dot*L1)**2 + c3*(1/2)*(y1_dot - y3_dot + teta_dot*L3/2 + gama_dot*L1)**2 + c4*(1/2)*(y1_dot - y4_dot + teta_dot*L3/2 - gama_dot*L2)**2 + c5*(1/2)*(y1_dot - y5_dot - teta_dot*L3/2 - gama_dot*L2)**2

## Força dissipativa
#Q = [- c2*y2_dot, - c3*y3_dot, - c4*y4_dot, - c5*y5_dot, - c1*y1_dot, - c6*gama_dot, - c7*teta_dot]

##### Lagrangeano 
L = Ec - Ep
L = sympy.simplify(L)

#### Derivadas em relacao aos deslocamentos
#### Derivadas em relacao as velocidades
#### Equecoes do movimento para cada grau de liberdade 
dL_dy = []
dL_dy_dot = []
#dtdL = []
dEd_dy_dot = []
graus = []

for i in range(len(var)):
    dL_dy.append(0)
    dL_dy_dot.append(0)
    dEd_dy_dot.append(0)
    #dtdL.append(0)
    graus.append(0)
    dL_dy[i] = sympy.diff(L, var[i])
    dL_dy_dot[i] = sympy.diff(L, var_dot2[i])
    #dtdL[i] = sympy.diff(dL_dy_dot[i], t_)
    dEd_dy_dot[i] = sympy.diff(Ed,var_dot[i])
    
    graus[i] = dL_dy_dot[i] - dL_dy[i] + dEd_dy_dot[i] #Equaçao do movimento da qual se extrai as matrizes M, K e C
    
    print( dEd_dy_dot[i],'\n')

  
    
### Montando Matriz de Rigidez
K = []

for i in range(len(var)):
    k = []
    for j in range(len(var)):
        k.append(0)
    K.append(k)
    

for i in range(len(var)):
    for j in range(len(var)):
       
        K[i][j] = float(graus[i].coeff(var[j],1))

for i in range(len(var)): #Verificando Simetria
    for j in range(len(var)):
        if K[i][j] == K[j][i]:
            pass
        else:
            print('Erro: Matriz K assimétrica')
            break
print('A matriz K é simétrica') 



K = np.reshape(K, (len(var),len(var)))


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
       
        C[i][j] = float(graus[i].coeff(var_dot[j],1))

for i in range(len(var)): #Verificando Simetria
    for j in range(len(var)):
        if C[i][j] == C[j][i]:
            pass
        else:
            print('Erro: Matriz de Amortecimento assimétrica')
            break
print('A matriz C é simétrica') 

        

C = np.reshape(C, (len(var),len(var)))
#print(C)

### Montando Matriz de Massa

M = []

for i in range(len(var)):
    m = []
    for j in range(len(var)):
        m.append(0)
    M.append(m)
    

for i in range(len(var)):
    for j in range(len(var)):
       
        M[i][j] = float(graus[i].coeff(var_dot2[j],1))
#print(M)
for i in range(len(var)): #Verificando Simetria
    for j in range(len(var)):
        if M[i][j] == M[j][i]:
            pass
        else:
            print('Erro: Matriz de Massa assimétrica')
            break
print('A matriz M é simétrica')

M = np.reshape(M, (len(var),len(var)))


##### Resolvendo problema de autovalores e autovetores oara calcular frequencias naturais e modos de vibrar
       
[wn, B] = eigh(K,M)
wn = np.sqrt(wn**2) #Obtendo o módulo das frequências naturais ao quadrado 
# #print(wn)
# #print(B) 
 

freq_n = np.sqrt(((wn)**2)**(1/2)) #Frequencia natural



meff = []
for i in range(len(var)): #Normalizando Matriz de massa pelo matriz de modal (autovetores)
    meff.append(0)
    meff[i] = np.dot(np.dot(B[:,i].T,M),B[:,i])

meffdiag = np.diag(1/np.sqrt(meff)) #Matriz diagonal com as massas efetivas meff para cada grau de liberdade

Ceff = []
for i in range(len(var)): #Normalizando Matriz de Amortecimento pela matriz modal
    Ceff.append(0)
    Ceff[i] = np.dot(np.dot(B[:,i].T,C),B[:,i]) #vetor de amortecimento efetivo [2.qsi(i).wn(i)]

Ceffdiag = np.diag(1/np.sqrt(Ceff))
#B = B*np.diag(1/np.sqrt(meff));

qsi = []
for i in range(len(var)): #Calculando taxa de amortecimento para cada grau de liberdade
    qsi.append(0)
    qsi[i] = Ceff[i]/(2*freq_n[i])      

freq_d = []
for i in range(len(var)): #Calculando frequencia natural amortecida
    freq_d.append(0)
    if qsi[i] <= 1:
        freq_d[i] = freq_n[i]*np.sqrt(1-qsi[i])
    else:
        freq_d[i] = freq_n[i]*np.sqrt(qsi[i]-1)
        
meff = np.diag(meff) 
for i in range(len(var)): #Verificando Simetria
    for j in range(len(var)):
        if meff[i][j] == meff[j][i]:
            pass
        else:
            print('Erro: Matriz meff assimétrica')
            break
print('A matriz meff é simétrica')
#### Constante Ksi de amortecimento


#### Condições iniciais
x0 = [0,0,0,0,0,0,0]
v0 = [0,100,100,0,0,0,0]

y0 = np.dot(np.dot(B.T,M),x0) #Condições iniciais em coordenadas modais
y0dot = np.dot(np.dot(B.T,M),v0)

t_inc0 = 0.05 #segundos
t_end0 = 2 #segundos    
t = np.arange(0, t_end0, t_inc0)

y = []
exp_qsi = []
for i in range(len(var)): #Pré-calculando primeiro termo da função resposta, "mágica"
    exp_qsi.append(0)
    exp_qsi[i] = np.exp(-qsi[i]*freq_n[i]*t)
    
for r in range(len(var)): #Obtendo a resposta no tempo em coordenadas modais para cada grau de liberdade
    y.append(0)
    y[r] = exp_qsi[r]*(y0[r]*np.cos(freq_d[r]*t) + ((y0dot[r] + qsi[r]*freq_n[r]*y0[r])/freq_d[r])*np.sin(freq_d[r]*t))
    
x = np.dot(B,y) #Convertendo resposta de coordenadas modais para físicas
x1 = x[0][:]
x2 = x[1][:]
x3 = x[2][:]
x4 = x[3][:]
x5 = x[4][:]
x6 = x[5][:]
x7 = x[6][:]

for i in range(len(var)):
    
      plt.plot(t,x[i][:])

# plt.plot(t,x6)
# plt.plot(t,x7)
np.savetxt('7GDL_RespostaAMORTECIDA_RodasLaterais', np.transpose([x1,x2,x3,x4,x5,x6,x7]), fmt='%1.5f')



