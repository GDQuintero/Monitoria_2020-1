#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 15:25:27 2020

@author: gustavo
"""

import numpy as np
import random as rd

teste = str(input("Digite o item (a, b, c ou d): "))

if (teste == "a"):
    N = 128; nf = 1
elif (teste == "b"):
    N = 128; nf = 4
else:
    nf = 10
    N = int(input("Digite o valor de N (128, 256, 512, 1024 ou 2048): "))

aux = np.zeros(N+1); T = 1.0; dx = 1.0 / N; dt = dx; lam = N; M = N

t = np.linspace(0.0, 1.0, M+1); barra = np.linspace(0.0, 1.0, N+1)

uold = np.zeros(N+1); unew = np.zeros(N+1); p = np.zeros(nf)

Adiag = np.zeros(N-1); Asub = np.zeros(N-2); b = np.zeros(N-1)

z = np.zeros(N-1); y = np.zeros(N-1); x = np.zeros(N-1); B = np.zeros(nf)

u = np.zeros((N+1,nf)); ut = np.zeros(N+1); A = np.zeros((nf,nf))

def r(t):
    return 10.0 * (1.0 + np.cos(5.0 * t))
    
def gh(x,aux,h,pi):
    aux[:] = 0.0
    for i in range(len(x)):
        if ((pi - (0.5 * h)) <= x[i]) and (x[i] <= (pi + (0.5 * h))):
            aux[i] = 1.0 / h
        else:
            aux[i] = 0.0
    return aux

def f(t,x,pi):
    return r(t) * gh(x, aux, dx, pi)

def u0(x):
    return 0.0

def g1(t): 
    return 0.0

def g2(t):
    return 0.0

def ldl_denso(A,n):
    
    D = np.zeros(n); L = np.zeros((n,n))
   
    #Primeira iteracao
    D[0] = A[0,0]
    for i in range(n):
        L[i,i] = 1.0
        
    #1's na diagonal
    for i in range(1,n):
        L[i,0] = A[i,0] / D[0]
        
    #Iteracoes intermedias
    for j in range(1,n):
        soma = 0.0
        
        for s in range(j):
            soma = soma + D[s]*(L[j,s]**2)
        D[j] = A[j,j] - soma
        
        for i in range(j+1,n):
            soma = 0.0
            for s in range(j):
                soma = soma + D[s]*L[i,s]*L[j,s]
            L[i,j] = (1.0 / D[j]) * (A[i,j] - soma)  
    
    return D, L
            
def solver_denso(D,L,B,n):
    
    X = np.zeros(nf); Y = np.zeros(nf); Z = np.zeros(nf)
    
    Z[0] = B[0]
    for i in range(1,n):
        soma = 0.0
        for s in range(i):
            soma = soma + L[i,s]*Z[s]
        Z[i] = B[i] - soma
    
    for i in range(n):
        Y[i] = Z[i] / D[i]
    
    X[n-1] = Y[n-1]
    for i in range(n-2,-1,-1):
        soma = 0.0
        for s in range(n-i-1):
            soma = soma + L[n-s-1,i]*X[n-s-1]
        X[i] = Y[i] - soma
        
    return X

def ldl_sparso(Adiag,Asub,n):
	d = np.zeros(n); l = np.zeros(n-1)
	d[0] = Adiag[0]; l[0] = Asub[0] / d[0]

	for i in range(1,n-1):
		d[i] = Adiag[i] - d[i-1] * (l[i-1]**2)
		l[i] = Asub[i] / d[i]

	d[n-1] = Adiag[n-1] - d[n-2] * (l[n-2]**2)

	return d, l

def solver_sparso(l,d,b,x,y,z,n):
	
	z[0] = b[0]

	for i in range(1,n):
		z[i] = b[i] - l[i-1] * z[i-1]

	for i in range(n):
		y[i] = z[i] / d[i]

	x[n-1] = y[n-1]

	for i in range(n-2,-1,-1):
		x[i] = y[i] - l[i] * x[i+1]
		
	return x


def erro_quadratico(ut,u,a,n):
    soma1 = 0.0; soma2 = 0.0
    for i in range(1,n-1):
        for j in range(nf):
            soma2 = soma2 + a[j]*u[i,j]
        soma1 = soma1 + (ut[i] - soma2)**2
        soma2 = 0.0
    E2 = np.sqrt(dx*soma1)
    return E2
            
Adiag[:] = 1.0 + lam 
Asub[:] = -0.5 * lam

d, l = ldl_sparso(Adiag,Asub,N-1)

if (teste == "a"):
    p[0] = 0.35
elif (teste == "b"):
    p = [0.15, 0.3, 0.7, 0.8]
else:
    p = [0.15, 0.2, 0.3, 0.35, 0.5, 0.6, 0.7, 0.73, 0.85, 0.9]
  
    
for i in range(nf): 
    uold[:] = u0(barra); porcentaje = int(100 / nf)
    
    for k in range(M):
        uold[0] = g1(t[k]); uold[N] = g2(t[k])
        fk = f(t[k], barra,p[i]); fkk = f(t[k+1], barra, p[i])
        b[0] = (0.5 * lam * (g1(t[k]) + g1(t[k+1]) + uold[2] - 2.0*uold[1])) + uold[1] + ((0.5 * dt) * (fk[1] + fkk[1]))
        b[1:-1] = (0.5 * lam * (uold[1:-3] -2.0*uold[2:-2] + uold[3:-1])) + uold[2:-2] + (0.5 * dt * (fk[2:-2] + fkk[2:-2]))
        b[N-2] = (0.5 * lam * (g2(t[k]) + g2(t[k+1]) + uold[-3] - 2.0*uold[-2])) + uold[-2] + ((0.5 * dt) * (fk[-2] + fkk[-2]))
        unew = solver_sparso(l,d,b,x,y,z,N-1)
        uold[1:-1] = unew.copy(); x[:] = 0.0; y[:] = 0.0; z[:] = 0.0
        
    u[:,i] = uold.copy(); uold[:] = 0.0
    
    print(str((i+1)*porcentaje)+"%","concluido")
    
if (teste == "a"):
    ut[:] = 7.0 * u[:,0]
    a1 = np.inner(ut, u[:,0]) / np.inner(u[:,0], u[:,0])  
    print()
    print("a1 = ", a1)
    
elif (teste == "b"):
    ut[:] = 2.3*u[:,0] + 3.7*u[:,1] + 0.3*u[:,2] + 4.2*u[:,3]
    for i in range(nf):
        B[i] = np.inner(ut, u[:,i])
        for j in range(nf):
            A[i,j] = np.inner(u[:,i], u[:,j])
            
    D, L = ldl_denso(A,nf)
    a = solver_denso(D,L,B,nf)
    
    print()
    for i in range(nf):
        print("a" + str(i+1)," = ", a[i])
    print()
    print("Erro quadratico: ", erro_quadratico(ut,u,a,N))
        
else:
    file = open('malha.txt','r')
    dados = file.readlines()
    k = int(2048 / N)
    
    for i in range(N):
            ut[i] = float(dados[int(i*k)])
            
    if (teste == "d"):
        epsilon = 0.01
        R = 2.0*(rd.random() - 0.5)
        ut[:] = (1.0 + R*epsilon)*ut[:]        
        
    for i in range(nf):
        B[i] = np.inner(ut, u[:,i])
        for j in range(nf):
            A[i,j] = np.inner(u[:,i], u[:,j])
            
    D, L = ldl_denso(A,nf)
    a = solver_denso(D,L,B,nf)
    print()
    for i in range(nf):
        print("a" + str(i+1)," = ", round(a[i],4))
    print()
    print("Erro quadratico: ", erro_quadratico(ut,u,a,N))
    
    file.close()

print(A)
AA = np.array([[4.,-1.,1.],[-1.,4.25,2.75],[1.,2.75,3.5]])
bb = np.ones(3)


    
    

