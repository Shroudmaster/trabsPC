from math import *
from matplotlib import *
import numpy as np
import sys

def coef(b, passos, a = 0):
	global natPrey

	h = (b - a) / passos
	for i in range(1, passos - 1):
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

	# inicialização 
	funx = fx(X, Y, 0)
	funy = fy(X, Y)
	atualx = X + (h/2) * fx(X + h, Y + (h * funx), 0)
	atualy = Y + (h/2) * fy(X + h, Y + (h * funy))

	print(str(atualx) + ", " + str(atualy))
	# muitas inicializações possiveis AAAAA
	#atualx = X + h * fx(X + h/2, Y + h/2, 0)
	#atualy = Y + h * fx(X + h/2, Y + h/2, 0)
	print(str(atualx) + ", " + str(atualy))

	anteriorx = X
	anteriory = Y

	vPreyBDF2.append(atualx)
	vPredBDF2.append(atualy)

	for i in range (1, passos-1):

		#predicao RK2
		funx = fx(atualx, atualy, i)
		funy = fy(atualx, atualy)
		proximox = atualx + ((h/2) * fx(atualx + h, atualy + (h * funx), i))
		proximoy = atualy + ((h/2) * fy(atualx + h, atualy + (h * funy)))
		
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

deltat = 2 #deve ser suficiente
passo = 4

#valores iniciais, onde presa deve ser razoavelmente maior q predador MAS NAO MUITO, senao o numero de arrombadinhos cresce exponencialmente (embora ainda falte implementar um sistema de numero maximo de arrombadinhos caso ele exploda)
presa = 45
maxPresa = 200
predador = 10

'''as taxas de natalidade/mortalidade dos bichanos
to sem um base decente de como chutar esses numeros. As regras gerais que tenho visto sao: natPrey e mortPred devem ser os dois maiores (apesar de q eu n vejo nenhuma razao pros predadores terem uma taxa de mortalidade tao grande, mas enfim) e os outros dois devem ser tbm beeeeeeeeeee menores. Inclusive acho q o valor atual pode estar cagado'''
natPrey = 1.0 #lembrando que essa porra eh ciclica, ta la no texto grande
mortPrey = 0.1
mortPred= 0.9
natPred = 0.045

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

print(vPreyRK4)
print(vPredRK4)
print(vPreyBDF2)
print(vPredBDF2)
print(vTempo)

