from math import *
import matplotlib.pyplot as plt
import numpy as np
import sys

# implementacao do crank nicolson retirada do livro 
# Numerical Analysis, 9th ed., pgs. 735-736
def crank_nicolson(l, T, alpha, m, N):
    # TODO: Inicializar um caralho de variaveis

    # step 1
    h = 1/m
    k = T/N
    gamma = (alpha**2 * k)/h**2
    w[m] = 0
    # step 2
    for i in range (1, m-1):
        # set wi = f(ih) [Initial Values]
    ### steps 3-11 according to Algorithm 6.7
    # step 3
    l[1] = 1 + gamma
    u[1] = -gamma/(2*l1)
    # step 4
    for i in range (2, m-2):
        l[i] = 1 + gamma + ((gamma * u[i-1])/2)
        u[i] = -gamma/(2*l[i])
    # step 5
    l[m-1] = 1 + gamma + (gamma*u[m-2])/2
    # step 6, for j do 7~11
    for j in range(1, n):
        # step 7
        t[j] = j * k
        z[1] = (((1 - gamma) * w[1]) + (gamma/2) * w[2])/l[1]
        # step 8
        for i in range(2, m-1):
            z[i] = ()(((1 - gamma) * w[1]) + (gamma/2) * (w[i+1] + w[i-1] + z[i-1]))/l[i]
        # step 9
        w[m-1] = z[m-1]
        # step 10, a complexidade disso é uma caquinha
        for i in range(m-2, 1, -1):
            w[i] = z[i] - (u[i]*w[i+1])
        # step 11
            # OUTPUT T