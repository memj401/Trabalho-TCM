import numpy as np
import time
import math
import copy
import matplotlib.pyplot as plt
import matplotlib.animation as animation


n_Fourier = 0.25 #Número de Fourier  definido a fim de manter a estabilidade da simulação

geracao_interna = 10 #Valor estipulado no trabalho

k = 81  #condutividade termica do aço

tamanho = 3 #comprimento do objeto simulado, em metros


#Funcão utilizada para desenhar o i-ésimo frame da animação
def animar(i,animacao):
  quadro = plt.imshow(animacao[i], cmap='hot', interpolation='nearest', extent=[0,tamanho,0,tamanho])
  return quadro
	
#Função que realiza o método das diferenças finitas	
def diferencas_finitas(numero_nos):
	inicio = time.process_time()

	delta_x = tamanho/(numero_nos)
	x_gap  = math.floor(numero_nos/3) #Coordenada X  da matriz da malha obde se encontra o gap no caso 7
	y_gap  = math.floor(numero_nos/3) #Coordenada Y da matriz da malha obde se encontra o gap no caso 7
	intervalo_gap_x = round(numero_nos/3) #Intervalo de Coordenadas em X no qual se encontra o gap do caso 7
	intervalo_gap_y = round(numero_nos/3) #Intervalo de Coordenadas em X no qual se encontra o gap do caso 7

	matriz = np.zeros([numero_nos,numero_nos], dtype = float)

	matriz[0, 1:numero_nos] = 100 #Condição de Contorno do Topo
	matriz[numero_nos - 1, 1:numero_nos] = 100 #Condição de Contorno do Fundo
	matriz[:numero_nos, 0] = 100 #Condição de Contorno da Esquerda
	matriz[:numero_nos, numero_nos - 1] = 300 #Condição de Contorno da Direita

	matriz[y_gap, x_gap :(x_gap + intervalo_gap_x)] = 20 #Condição de Contorno do Topo do Gap
	matriz[y_gap + intervalo_gap_y , x_gap: (x_gap + intervalo_gap_x)] = 15 #Condição de Contorno do Fundo do Gap
	matriz[y_gap : (y_gap + intervalo_gap_y + 1), x_gap] = 10 #Condição de Contorno da Esquerda do Gap
	matriz[y_gap : (y_gap + intervalo_gap_y + 1), (x_gap + intervalo_gap_x)] = 30 #Condição de Contorno da Direita do Gap
	matriz[y_gap + 1 : (y_gap + intervalo_gap_y), x_gap + 1 : (x_gap + intervalo_gap_x)] = 0 #Condição de Contorno do interior do Gap


	matriz_atual = copy.deepcopy(matriz)
	matriz_anterior = copy.deepcopy(matriz)
	parada = False
	quadros = 0
	animacao = []

	while parada != True:
		quadros += 1
		animacao.append(copy.deepcopy(matriz_atual)) 
		if(quadros%100 == 0):
			print(quadros)
		for x in range(numero_nos):
			for y in range(numero_nos):
				if (x == 0 or y == 0) or (x == (numero_nos - 1) or y == (numero_nos - 1)):
					continue
				elif (x >= x_gap and x <= (x_gap + intervalo_gap_x)) and (y >= y_gap and y <= (y_gap + intervalo_gap_y)):
					continue
				matriz_anterior[y,x] = copy.deepcopy(matriz_atual[y,x])	
				matriz_atual[y,x] = n_Fourier*(matriz_atual[y,x+1] + matriz_atual[y,x-1] + matriz_atual[y+1,x] + matriz_atual[y-1,x]) 
				+ (1 - 4*n_Fourier)*matriz_atual[y,x] 
				+ (geracao_interna*(delta_x**2)/k)
				
				if ((abs(matriz_atual[y,x] - matriz_anterior[y,x])/matriz_atual[y,x]) <= 0.0000001): 
				#Condição de Parada (Diferença entre dois pontos do estado atual e do anterior menor que 0.00001%)
					parada = True
	tempo_consumido = time.process_time() - inicio
	print('Total de Frames: ' + (str)(quadros)) 
	return animacao,quadros,tempo_consumido,matriz_atual


#Função que gera o vídeo da simulação e mostra a imagem do estado permanente
def gerar_midias(numero_nos):
	resultado = diferencas_finitas(numero_nos)
	animacao = resultado[0]
	quadros = resultado[1]
	figura, eixos = plt.subplots()
	eixos.set(title = 'Distribuição de Temperatura - Caso 7 (Malha de ' + (str)(numero_nos**2) + ' nós)', xlabel = 'Eixo X', ylabel = 'Eixo Y' )
	imagem = eixos.imshow(animacao[quadros-1], cmap='hot', interpolation='nearest', extent=[0,tamanho,0,tamanho])
	figura.colorbar(imagem)
	video = animation.FuncAnimation(figura,animar,frames=quadros,fargs=(animacao,))
	writermp4 = animation.FFMpegFileWriter(fps=30, bitrate=1800)
	video.save('sim.mp4', writer=writermp4)
	plt.show()
	
#Função que define o custo operacional, em s, para um série de valores de nós
def custo_operacional(valores_nos):
	tempos = []
	malhas = []
	
	for i in valores_nos:
		tempos.append(diferencas_finitas(i)[2])
	for j in valores_nos:
		malhas.append(j**2)
	
	plt.plot(malhas,tempos)
	plt.title('Custo Operacional por Número de Nós na Malha')
	plt.xlabel('Numeros de Nós')
	plt.ylabel('Tempo de Execução (s)')
	plt.savefig('CustoOperacional.png', bbox_inches = 'tight')

#Função que obtém os perfis de temperatura para X = 0.5,1.5 e 2.5 e Y = 1.5 e 2.5
def obter_perfis(numero_nos):
	y_observado = [math.floor(numero_nos/(tamanho/1.5)), math.floor(numero_nos/(tamanho/2.5))] #Valores de x defindos para análise
	x_observado = [(math.floor(numero_nos/(tamanho/0.5))), (math.floor(numero_nos/(tamanho/1.5))),math.floor(numero_nos/(tamanho/2.5))] #Valores de y definidos para análise
	valores_y = (1.5,2.5)
	valores_x = (0.5,1.5,2.5)
	intervalo = np.linspace(0,tamanho, num = numero_nos)
	estado_permanente = diferencas_finitas(numero_nos)[3]
	
	for i in range(len(x_observado)):
		temperaturas = np.array(estado_permanente[:,x_observado[i]])
		plt.plot(intervalo,temperaturas)
		plt.title('Perfil de Temperatura para X = ' + (str)(valores_x[i]))
		plt.xlabel('Distância em Y (m)')
		plt.ylabel('Temperatura (K)')
		plt.savefig('PerfilTempX'+ (str)(valores_x[i]) + '.png', bbox_inches = 'tight')
		plt.close('all')

	for i in range(len(y_observado)):
		temperaturas = np.array(estado_permanente[y_observado[i],:])
		plt.plot(intervalo,temperaturas)
		plt.title('Perfil de Temperatura para Y = ' + (str)(valores_y[i]))
		plt.xlabel('Distância em X (m)')
		plt.ylabel('Temperatura (K)')
		plt.savefig('PerfilTempY'+ (str)(valores_y[i]) + '.png', bbox_inches = 'tight')
		plt.close('all')


