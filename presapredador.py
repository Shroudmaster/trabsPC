from math import *
import matplotlib.pyplot as plt
import numpy as np
import sys

def coef(b, passos, a = 0):
	global natPrey

	h = (b - a) / passos
	for i in range(1, passos):
		t = h * i
		aux = (1.5 + np.sin(t)) * natPrey
		vTempo.append(t)
		vAlpha.append(aux)
	return

'''natPrey = natPrey * (maxPrey - v[-1])
'''
def fx(x, y, i):
	global natPrey, mortPrey, maxPresa, vAlpha
	if(x > maxPresa):
		x = maxPresa
	res = (vAlpha[i] * (maxPresa - x)) - (mortPrey * x * y)
	return res

def fy(x, y):
	global natPred, mortPred
	res = -(mortPred * y) + (natPred * x * y)
	return res

def BDF2(b, xVec, yVec, passos, a = 0):
	# atual, anterior, proximo
	h = (b - a) / passos
	
	X = xVec[-1]
	Y = yVec[-1]

	### Contas antigas ###
	#atualx = X + (h/2) * fx(X + h, Y + (h * funx), 0)
	#atualy = Y + (h/2) * fy(X + h, Y + (h * funy))
	#print(str(atualx) + ", " + str(atualy))

	#a1 = 1/2
	#a2 = 1/2
	#p1 = 1
	#q1 = 1

	# inicialização 
	funx = fx(X, Y, 0)
	funy = fy(X + h, Y + funx * h)

	atualx = X + h * ((funx/2) + (funy/2))
	atualy = Y + h * ((funx/2) + (funy/2))

	#print(str(atualx) + ", " + str(atualy))

	anteriorx = X
	anteriory = Y

	vPreyBDF2.append(atualx)
	vPredBDF2.append(atualy)

	for i in range (1, passos-1):

		#predicao RK2
		funx = fx(atualx, atualy, i)
		funy = fy(atualx + h, atualy + funx * h)
		proximox = atualx + h * ((funx/2) + (funy/2))
		proximoy = atualy + h * ((funx/2) + (funy/2))
		
		# tem muitos metodos
		#proximox = atualx + (h * fx(atualx + h/2, atualy + h/2))
		#proximoy = atualy + (h * fy(atualx + h/2, atualy + h/2))
		

		# correcao
		proximox = (1/3) * ((4 * atualx) - anteriorx + (2 * h * fx(proximox, proximoy, i)))
		proximoy = (1/3) * ((4 * atualy) - anteriory + (2 * h * fy(proximox, proximoy)))

		anteriorx = atualx
		anteriory = atualy

		atualx = proximox
		atualy = proximoy

		vPreyBDF2.append(atualx)
		vPredBDF2.append(atualy)

	return

def RK4(b, xVec, yVec, passos, a = 0): #f(t, x(t)) = ax(t) - bx(t)y(t)
	h = (b - a) / passos
	for i in range (passos - 1):
		X = xVec[-1]
		Y = yVec[-1]
		
		KX1 = h * fx(X, Y, i)
		KY1 = h * fy(X, Y)
		
		KX2 = h * fx(X + (KX1/2.0), Y + (KY1/2.0), i)
		KY2 = h * fy(X + (KX1/2.0), Y + (KY1/2.0))
		
		KX3 = h * fx(X + (KX2/2.0), Y + (KY2/2.0), i)
		KY3 = h * fy(X + (KX2/2.0), Y + (KY2/2.0))
		
		KX4 = h * fx(X + KX3, Y + KY3, i)
		KY4 = h * fy(X + KX3, Y + KY3)
		
		X = X + (KX1 + (2 * KX2) + 2 * KX3 + KX4)/6.0
		Y = Y + (KY1 + (2 * KY2) + 2 * KY3 + KY4)/6.0
		
		vPreyRK4.append(X)
		vPredRK4.append(Y)
	return

deltat = 200 #deve ser suficiente
passo = 800

#valores iniciais, onde presa deve ser razoavelmente maior q predador MAS NAO MUITO
presa = 40
maxPresa = 100
predador = 10

# taxas de mortalidade e natalidade dos bichanos
natPrey = 0.50 # lembrando que de acord com a especificação isso é cíclico
mortPrey = 0.05
mortPred= 0.45
natPred = 0.0225

#aqui a gente vai armazenar os resultados pra dps plotar
vPreyRK4 = [presa]
vPredRK4 = [predador]

vPreyBDF2 = [presa]
vPredBDF2 = [predador]

vAlpha = [natPrey]
vTempo = [0]

coef(deltat, passo)
RK4(deltat, vPreyRK4, vPredRK4, passo)
BDF2(deltat, vPreyBDF2, vPredBDF2, passo)

'''
print(vPreyRK4)
print(vPredRK4)
print(vPreyBDF2)
print(vPredBDF2)
print(vTempo)
'''

plt.figure("Modelo Lotka-Volterra", figsize=(8,5))
plt.title("Análise de Resultados")
plt.plot(vTempo, vPreyBDF2, label='presaBDF2')
plt.plot(vTempo, vPredBDF2, label='predadorBDF2')
plt.xlabel('tempo')
plt.ylabel('populacao')
plt.legend()

plt.figure("Modelo Lotka-Volterra", figsize=(8,5))
plt.title("Análise de Resultados")
plt.plot(vTempo, vPreyRK4, label='presaRK4')
plt.plot(vTempo, vPredRK4, label='predadorRK4')
plt.xlabel('tempo')
plt.ylabel('populacao')
plt.legend()
plt.show()