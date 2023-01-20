# Algoritmos para encontrar soluções de Sistemas de Equações Lineares
#_________________________________________________________
# Universidade Federal de Santa Catarina
# Departamento de Engenharias da Mobilidade
# Curso de Cálculo Numérico
# Prof. Alexandre Zabot
# https://zabot.paginas.ufsc.br
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
import sl
import math
import numpy as np



def norma(v):
  n = 0
  for xi in v:
    n += xi*xi
  return math.sqrt(n)



def ponto_fixo( G, x0, tol, Imax, verbose=True ):
  """
  Algoritmo de Ponto Fixo para Sistemas de Equações não lineares
  
  Recebe:
    A       => Função que calcula funções de ponto fixo
    x0      => Vetor com chute inicial
    tol     => Tolerância
    Imax    => Número máximo de iterações
    verbose => Imprime passo a passo?
    
  Retorna:
    número de iterações
    vetor com solução
    erro absoluto
  """
  
  N = len(x0)     # Número de variáveis/equações
  
  Xold = np.copy(x0)     # Chute inicial
  Xnew = np.copy(x0)     # Novo valor da iteração
  


  stop=False
  iteration = 0
  
  

  if verbose:  
    print("%2d" % iteration, end="")
    print(Xnew,"----")
  
  while not stop:
    
    Xnew = G( Xold, verbose )
    
    
    e_abs = norma( Xnew - Xold )
    
    if e_abs<tol:           stop=True
    elif iteration>=Imax:   stop=True
    
    Xold = np.copy(Xnew)
    iteration += 1
    
    if verbose:
      print("%2d" % iteration, end="")
      print(Xnew, e_abs )
  
  return iteration-1,Xnew,e_abs
    
  




def printmatrix(M):
  for i in range(len(M)):
    for j in range(len(M[i])):
      print( "%13.6g" % M[i][j], end="")
    print






def newton(matriz, x0, tol, Imax, verbose=True ):
  """
  Algoritmo de Newton para Sistemas de Equações não lineares
  
  Recebe:
    matriz  => Função que calcula matriz aumentada com coef. da jacobiana e -F
    x0      => Vetor com chute inicial
    tol     => Tolerância
    Imax    => Número máximo de iterações
    verbose => Imprime passo a passo?
    
  Retorna:
    número de iterações
    vetor com solução
    erro absoluto
  """

  
  Xold = x0.copy() 
  Xnew = Xold.copy() # Novo valor da iteração


    
  stop=False
  iteration = 1
  while not stop: 
    
    if verbose:
      print("\n---------------------------------------")
      print("i = %d " % iteration                      )
    
    # ----------------------------------------------------------------------------
    # Calcula a matriz aumentada
    A = matriz(Xold)
    printmatrix( A )
    
    
    # ----------------------------------------------------------------------------
    # RESOLVE O SISTEMA LINEAR
    if verbose:
      print("\n\n------------------------------------")
      print("Etapa de Resolução do Sistema linear para obter y a partir de Jy = -F")
      print("\n\n-----")
    y = sl.eliminacao_gaussiana(A,False)
    if verbose: print( "y_%02d" % (iteration-1), y )
    
    
    # ----------------------------------------------------------------------------
    # Encontra a nova iteração para X
    Xnew = Xold + y
    if verbose:
      print("X_%02d" % (iteration-1), Xold )
      print("X_%02d" % iteration, Xnew )
    
    
    # ----------------------------------------------------------------------------
    # Calcula o erro e atualiza os vetores
    eA = norma( Xnew - Xold )
    if verbose: print("Erro Absoluto = ",eA   )
  
  
    if   eA<tol      : stop=True
    elif iteration>=Imax: stop=True
    
    Xold = Xnew
    iteration += 1
  
  return iteration-1,Xnew,eA


