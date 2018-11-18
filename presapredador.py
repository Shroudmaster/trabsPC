from math import *
from matplotlib import *
import numpy as np
import sys
'''
	Ok, esse eh o cao. A gente vai montar um sistema presa-predador baseado numa funcao probabilistica que deve ser meio zoada de implementar, mas vambora.
	A gente precisa setar alpha beta gama delta (go kappa kappa!), delta t e pah.
	Pegar pseudo codigo do runge kutta ssj4 pg288(388?) do livro no google drive
	Implementar BDF2 (inicializacao e predicao via runge kutta lvl 2) mas isso n vai matar ngm
	Nao sei se o python tem isso ou se sequer ajuda mas sugiro fazer essas porras todas serem macros, pq a gente vai iterar esses metodos ate o cu fazer bico e tal
	A taxa de natalidade (alpha) eh ciclica (1.5 + sin(t)) * alpha). Entao tipo, if (primavera) ++(coelhos.transa);
	Existe um limite no numero maximo de coelhos mesmo sem predadores no sistema. Logo o aumento de presas na falta de predadores fica alpha*(n_max_coelhos - n_atual_coelhos)
	nada vai ser odiosamente copiado via scipy, infelizmente
	Plotar esta merda toda, hence matplotlib. Fazer uns graficos bonitoes pra ficar maneiro no relatorio
	Nunca mais esquecer fone de ouvido quando eu for trabalhar nessa merda
	Mais alguma coisa?
'''
def fx(x, y):
	global natPrey, mortPrey
	res = (natPrey * x) - (mortPrey * x * y)
	return res

def fy(x, y):
	global natPred, mortPred
	res = -(mortPred * y) + (natPred * x * y)
	return res

def RK4(b, xVec, yVec, passos, a = 0): #f(t, x(t)) = ax(t) - bx(t)y(t)
	h = (b - a) / passos;
	for i in range (passos):
		t = h * i
		X = xVec[-1]
		Y = yVec[-1]
		
		KX1 = h * fx(X, Y)
		KY1 = h * fy(X, Y)
		
		KX2 = h * fx(X + (KX1/2.0), Y + (KY1/2.0))
		KY2 = h * fy(X + (KX1/2.0), Y + (KY1/2.0))
		
		KX3 = h * fx(X + (KX2/2.0), Y + (KY2/2.0))
		KY3 = h * fy(X + (KX2/2.0), Y + (KY2/2.0))
		
		KX4 = h * fx(X + KX3, Y + KY3)
		KY4 = h * fy(X + KX3, Y + KY3)
		
		X = X + (KX1 + (2 * KX2) + 2 * KX3 + KX4)/6.0
		Y = Y + (KY1 + (2 * KY2) + 2 * KY3 + KY4)/6.0
		
		vPrey.append(X)
		vPred.append(Y)
		vTempo.append(t)
	return

deltat = 2 #deve ser suficiente
passo = 4

#valores iniciais, onde presa deve ser razoavelmente maior q predador MAS NAO MUITO, senao o numero de arrombadinhos cresce exponencialmente (embora ainda falte implementar um sistema de numero maximo de arrombadinhos caso ele exploda)
presa = 50
predador = 10

#aqui a gente vai armazenar os resultados pra dps plotar
vPrey = [presa]
vPred = [predador]
vTempo = []

'''as taxas de natalidade/mortalidade dos bichanos
to sem um base decente de como chutar esses numeros. As regras gerais que tenho visto sao: natPrey e mortPred devem ser os dois maiores (apesar de q eu n vejo nenhuma razao pros predadores terem uma taxa de mortalidade tao grande, mas enfim) e os outros dois devem ser tbm beeeeeeeeeee menores. Inclusive acho q o valor atual pode estar cagado'''
natPrey = 1.0 #lembrando que essa porra eh ciclica, ta la no texto grande
mortPrey = 0.02
mortPred= 0.9
natPred = 0.03

RK4(deltat, vPrey, vPred, passo)

print(vPrey)
print(vPred)
print(vTempo)

