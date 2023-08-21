# -*- coding: utf-8 -*-
"""
Created on Sun Oct  3 08:17:02 2021

@author: Anderson
"""

"""Experimento V de LAB IV"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, root_scalar





def gaussiana(x, mi, A, sigma):

    return A*np.exp(-((x-mi)**2)/(2*(sigma**2)))


def rodando_gaussiana(arquivo,a=20):
    #Dados para Gaussiana:
    Y_max = arquivo.max()
    i = np.where(arquivo == Y_max)[0][0]
    minimo = max(i-a,0)
    maximo = min(i+a,len(arquivo))
    x_gauss, y_gauss = np.linspace(minimo, maximo, len(arquivo[minimo:maximo])),arquivo[minimo:maximo] 
    media = np.mean(x_gauss)
    sigma  = np.std(x_gauss)
    n = len(x_gauss)
    
    chute = [media, arquivo.max(), sigma]
    
    fit, cov = curve_fit(gaussiana, x_gauss, y_gauss, sigma = [1/np.sqrt(n)]*n, p0 = chute, absolute_sigma = True)
    mean, amplitude, sigma = fit
    return x_gauss, mean, amplitude, sigma, cov
    
    
    
 
def fit(): #CERTO
    #Figura:
    fig, ax = plt.subplots()
    diretorio = './Novo/'
    arquivos = os.listdir(diretorio)
    n = 1
    f = open('arquivo_dados.dat', 'w')
    
    for fname in arquivos:
      arquivo = np.loadtxt(diretorio + fname)
      Y_max = arquivo[20:].max()
      i = np.where(arquivo == Y_max)[0][0]
      if i>600 or n == 23:
          X, Y = np.linspace(20,len(arquivo[20:]), len(arquivo[20:])), arquivo[20:]
          ax.plot(X,Y)
         
          #Gaussiana
          x_gauss, mean, amplitude, sigma, cov = rodando_gaussiana(arquivo[20:])
      else:
          X, Y = np.linspace(0,len(arquivo), len(arquivo)), arquivo
          ax.plot(X,Y)
         
          #Gaussiana
          x_gauss, mean, amplitude, sigma, cov = rodando_gaussiana(arquivo,50)
      ax.plot(x_gauss,gaussiana(x_gauss,mean,amplitude,sigma),c='k')
      ax.set_xlim(1800,2000)
      #ax.set_ylim(0,500)
      num = fname.split('_')[1].replace('m','')
      f.write(num + ' ' + str(mean) + ' ' + str(np.sqrt(cov[0][0])) + '\n')
      n+=1
      print('---')
      ax.set_title("Ajustes Gaussianos dos picos dos sinais")
      ax.set_xlabel("Canais")
      ax.set_ylabel("Contagens")
    f.close()  
    return fig, ax


def reta_calibracao():
    """Parte para a calibracao"""
    arquivo2 = np.loadtxt("DadosCalibracao.txt")
    fig1, ax1 = plt.subplots()
    ax1.scatter(arquivo2[:,0], arquivo2[:,1]*5.486/146)
    coef = np.polyfit(arquivo2[:,0],arquivo2[:,1]*5.486/146,1, cov = True)
    f = np.poly1d(coef[0])
    ax1.plot(arquivo2[:,0], f(arquivo2[:,0]), c = 'k', label = f' a = ({coef[0][0]:.5f} +/- {np.sqrt(coef[1][0][0]):.5f}) MeV/canal')
    ax1.plot(arquivo2[:,0], f(arquivo2[:,0]), c = 'k', label = f' b = ({coef[0][1]:.2f} +/- {np.sqrt(coef[1][1][1]):.2f}) MeV')
    
    
    print("Os coeficientes são: {} e {} " .format(coef[0][0], coef[0][1]))
    print("E sua incerteza é : {} e {} " .format(np.sqrt(coef[1][0][0]), np.sqrt(coef[1][1][1])))
    
    ax1.set_title("Reta de calibração")
    ax1.set_ylabel("Energia (MeV)")
    ax1.set_xlabel("Canal")
    ax1.legend()
    return coef[0][0], np.sqrt(coef[1][0][0]), coef[0][1], coef[1][1][1]
    

def origem():
    file = np.loadtxt("arquivo_dados.dat")
    a, inc_a, b, inc_b = reta_calibracao()
    inn = np.sqrt(file[:,2]**2 + inc_a**2)  #Incerteza energia
    incertezas = np.sqrt(inn**2 + inc_b**2)
    
    #Plotando:
    fig3, ax3 = plt.subplots()
    ax3.scatter(file[:,1]*a + b, file[:,0])
    
    coef = np.polyfit(file[:,1]*a + b,file[:,0],3, cov = True)
    f = np.poly1d(coef[0])
    file[:,1].sort()
    
    
    ax3.errorbar(file[:,1]*a + b, f(file[:,1]*a + b), incertezas, c = 'k')
    ax3.set_title("Gráfico reetrância")
    ax3.set_xlabel("Energia (MeV)")
    ax3.set_ylabel("Distância (mm)")
    print("A reetrancia é: ", f(5.486))
    print(np.sqrt(coef[1][2][2]))
    return file, f(5.486),np.sqrt(coef[1][2][2]), file[:,1]*a + b, incertezas #a raiz e a incerteza do novo valor de d0


def calibracao_distancias(): 
    file, f_alpha, inc_f_alpha, energias, inc_energia= origem()
    
    energias = np.flip(energias)
    energias = np.concatenate((energias, [0, 0]))
    inc_energia = np.concatenate((inc_energia, [0, 0]))
    novas_distancias = file[:,0] - f_alpha
    novas_distancias.sort()
    novas_distancias = np.concatenate((novas_distancias, [43.169, 43.170]))
    inc_nova_dist = inc_f_alpha
    
    #Plotando:
    fig4, ax4 = plt.subplots()
    ax4.scatter(novas_distancias, energias)
    
    ax4.set_title("Energia X distância")
    ax4.set_xlabel("Distância (mm)")
    ax4.set_ylabel("Energia (MeV)")
    
    #ajustando uma curva:
    coef = np.polyfit(novas_distancias, energias, 3, cov = True)
    f = np.poly1d(coef[0])
    
    #achando a raiz:
    a, b, c = coef[0][0], coef[0][1], coef[0][2]
    raiz = root_scalar(f, x0 = 40, method = "newton", fprime = lambda x: a*3*x**2 + 2*b*x + c )
    print()
    print("A raiz é:", raiz.root)
    print()
    
    ax4.errorbar(novas_distancias, f(novas_distancias), inc_energia, inc_nova_dist,  c = 'k')
 
    
    
    #Queremos agora a curva de Brag, ou seja, a derivada da energia com relacao a distancia e, para tal, precisamos da
    #dos coeficientes da curva que ele ajustou:
    
    #Coeficientes:
    print("Os coeficientes para esse caso são: ", coef[0])
    print("As incertezas dos coeficientes sao: {}, {}, {} e {}" .format(np.sqrt(coef[1][0][0]), np.sqrt(coef[1][1][1]), np.sqrt(coef[1][2][2]), np.sqrt(coef[1][3][3])) )
    return coef[0][0], coef[0][1], coef[0][2], np.sqrt(coef[1][0][0]), np.sqrt(coef[1][1][1]), np.sqrt(coef[1][2][2]), novas_distancias, energias
    


def curva_bragg(): #REVER
    #Derivada da funcao de 3 grau obtida na funcao anterior
    a, b, c, in_a, in_b, in_c, x, energias = calibracao_distancias() #x sao as novas distancias
    
    #Plotando:
    dev_exp =  np.gradient(energias, x, edge_order = 1)
    fig5, ax5 = plt.subplots()
    ax5.scatter(x, -dev_exp)
    ax5.set_title("Perda de energia")
    ax5.set_xlabel("Distância (mm)")
    ax5.set_ylabel("dE/dx (MeV/mm")
    return -dev_exp, energias



def alcance_empirico():
    E_Am241 = ((5443*13.0 + 5486*84.5)/(13 + 84.5))/1000
    return 10*(0.005*E_Am241 + 0.285)*E_Am241**(3/2)


def programa():
    file = np.loadtxt("programa.txt")
    energia, stopping_power = file[:,0], file[:,1]*(1.1839/1000)
    
    #Plotando:
    Y, X = curva_bragg()
    fig7, ax7 = plt.subplots()
    ax7.plot(energia,stopping_power/10, label = "Dados da simulação")
    ax7.scatter(X,Y, 4, 'g', label = "Dados experimentais")
    ax7.legend()
    ax7.set_title("Simulação com o ASTAR  e experimento")
    ax7.set_xlabel("Energia MeV")
    ax7.set_ylabel("Perda de energia por comprimento (MeV/mm)")

    

    
