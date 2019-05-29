#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 15:50:13 2019

@author: Alan Henschel Costa
"""
import numpy as np
import copy

def gausss(A):
    n = len(A)

    for i in range(0, n):
        maxNum = abs(A[i][i])
        maxRow = i
        for k in range(i+1, n):
            if abs(A[k][i]) > maxNum:
                maxNum = abs(A[k][i])
                maxRow = k

        for k in range(i, n+1):
            tmp = A[maxRow][k]
            A[maxRow][k] = A[i][k]
            A[i][k] = tmp

        for k in range(i+1, n):
            if(A[i][i] == 0):
                continue
            c = -A[k][i]/A[i][i]
            for j in range(i, n+1):
                if i == j:
                    A[k][j] = 0
                else:
                    A[k][j] += c * A[i][j]

    x = [0 for i in range(n)]
    for i in range(n-1, -1, -1):
        if(A[i][i] == 0):
            continue
        x[i] = A[i][n]/A[i][i]
        for k in range(i-1, -1, -1):
            A[k][n] -= A[k][i] * x[i]
    return x


def gauss(A,b):
    
    gauss = copy.deepcopy(A)
    for i in range(len(gauss)):
        gauss[i].append(b[i])
    gauss = gausss(gauss)  
    return gauss

def gauss2(A):
    gauss = copy.deepcopy(A)
    gauss = gausss(gauss)
    return gauss
   
def formaPadrao(A,sinais,b,m,funcao):
    N_variaveis_adicionais = len(sinais)
    flag = 0
    for j in range(m):
        for i in range(N_variaveis_adicionais):
            if (sinais[j]=="<=" and j==i) or (sinais[j]=="<" and j==i):
                A[j].append(1.0)
            elif sinais[j]==">=" and j==i or sinais[j]==">" and j==i:
                A[j].append(-1.0)
                flag = 1
            else:
                A[j].append(0.0)
                
    for i in range(N_variaveis_adicionais):
        funcao.append(0.0)
    
    return A,N_variaveis_adicionais,flag


def transposta(matriz):
    transpose = np.zeros((len(matriz[0]),len(matriz)))
    for c in range(len(matriz)):
      for d in range(len(matriz[0])):
         transpose[d][c] = matriz[c][d]
         
    return transpose

def transposta2(matriz):
    transpose = []
    for i in range(len(matriz)-1,-1,-1):
        transpose.append(matriz[i])
         
    return transpose
         
def column(matrix, i):
    
    return [row[i] for row in matrix]

def solucaoOtima(Xb,indexacaoB,funcao,tipo):
    
    resultado = 0
    func = []
    
    for i in range(len(Xb)):
        print('\nX{} = {}'.format(i, Xb[i]))
        print('\n index={}=colunas={}'.format(i,indexacaoB[i]))
    for i in indexacaoB:
        func.append(funcao[i])
    for i in range(len(Xb)):
        resultado += Xb[i]*func[i]
    
    if tipo=="max":
        resultado = resultado * (-1)
    print("\n A  SOLUCAO EHHHH:",resultado)
    exit(0)

def passo1(B,b):
    Xb = gauss(B,b)
    print("\nXB:",Xb)
    return Xb
def passo2(B,funcao,indexacaoN,N,indexacaoB):
    transpose = transposta(B)
    multiplicador = []
    for i in range(len(transpose)):
        multiplicador.append(list(np.append(transpose[i],funcao[indexacaoB[i]])))
    multiplicador2 = gauss2(multiplicador)
    print("\nvetor multiplicador simplex:",multiplicador2)
    
    #2.2
    CN = []
    for j in range(len(N[0])):
        c = 0
        for i in range(len(multiplicador2)):
            c +=multiplicador2[i] * N[i][j]
        CN.append(funcao[indexacaoN[j]] - c)
        
    #2.3
    CNK = min(CN)
    indexCNK = CN.index(min(CN))
    print("CUSTOS",CN)
    print("VARIAVEL A ENTRAR NA BASE:",indexCNK)

    return CNK,indexCNK
def passo4(N,indexCNK,B):
    colunaN = column(N,indexCNK)
    y = gauss(B,colunaN)
    print("\nYYYYYYYYYYY",y)
    return y
def passo5(y,Xb):
    e = float("inf")
    indexY = 0         
    for i in range(len(Xb)):
        if y[i]> 0 and Xb[i]/y[i] < e:
            e = Xb[i]/y[i]
            indexY = i
    print("\nEEEEEEEEE:",e)
    print("\nEEEEEEEEEYYYYYY:",indexY)
    return e,indexY
    
def passo6(indexY,indexacaoB,indexacaoN,indexCNK,B,N,A): 
    print("INDEX SAIII",indexY)
    print("\nAQUIIIIIIIIIIIIIIIIIII\n")
    auxx = indexacaoB[indexY]
    auxx2 = indexacaoN[indexCNK]
    print("\nAUXXXXXXXXXXX",auxx)
    indexacaoB[indexY] = auxx2        
    print("\nCNKKKKKKK",indexCNK)
    indexacaoN[indexCNK] = auxx
    print("\nindexacaoB:",indexacaoB)
    print("\nindexacaoN:",indexacaoN)    
    auxx2 = []
    for i in indexacaoB:
        auxx2.append(column(A,i))        
    for i in range(len(auxx2)):
        for j in range(len(auxx2[0])):
            B[j][i] = auxx2[i][j]
            
    auxx2 = []
    for i in indexacaoN:
        auxx2.append(column(A,i))
                
    for i in range(len(auxx2)):
        for j in range(len(auxx2[0])):
            N[j][i] = auxx2[i][j]        
    print("\nMATRIZ AAAAAA",A)
    print("\nBASE:",B)
    print("\nNAO-BASE:",N)
    print("\nindexacaoB:",indexacaoB)
    print("\nindexacaoN:",indexacaoN)
    
            

def segunda_fase(A,sinais,b,funcao,m,n,N_variaveis_adicionais,tipo,indexacaoB = [],indexacaoN = [],B = [],N = []):
    print("\n\nENTRANDO NA SEGUNDA FASE IIIIIIIIII\n\n")
    if(len(indexacaoB) == 0):
        for j in range(m):
            t = []
            for i in range(n,n + N_variaveis_adicionais):
                if j ==0:
                    indexacaoB.append(i)
                t.append(A[j][i])
            B.append(t)
            l = []
            for i in range(n):
                if j == 0:
                    indexacaoN.append(i)
                l.append(A[j][i])
            N.append(l)
         
    
    print("\nBASE:",B)
    
    print("\nNAO-BASE:",N)

    print("\nindexacaoB:",indexacaoB)

    print("\nindexacaoN:",indexacaoN)
    
    condicao = True
    iteracao = 1
    while(condicao):
        #passo1
        Xb = passo1(B,b)
        #passo 2
        CNK,indexCNK = passo2(B,funcao,indexacaoN,N,indexacaoB)
        # passo 3
        if(CNK >= 0.000):
            print("\nSolucao Otima encontrada!!!!")
            condicao = False
            solucaoOtima(Xb,indexacaoB,funcao,tipo)
            break
        #passo 4 
        y = passo4(N,indexCNK,B)
        #passo5
        ver = sum(1 for number in y if number < 0)
        if(ver == len(y)):
            print("problema nao tem solucao otima finita")
            condicao = False
            solucaoOtima(Xb,indexacaoB,funcao,tipo)
            break
        else:
            e,indexY = passo5(y,Xb)
        
        #passo6
        passo6(indexY,indexacaoB,indexacaoN,indexCNK,B,N,A)
        print("\niteracao:::",iteracao)
        iteracao+=1
    
def primeira_fase(A,sinais,b,funcao,m,n,N_variaveis_adicionais,tipo):
    print("\n\nENTRANDO NA PRIMEIRA FASE IIIIIIIIIIIIIIIIII\n\n")
    A1 = copy.deepcopy(A)
    funcao2 = funcao.copy()
    n1 = len(A1[0])
    for j in range(m):
        for i in range(m):
            if i==j:
                A1[j].append(1.0)
            else:
                A1[j].append(0.0)
     
    for i in range(len(funcao)):
        funcao2[i] = 0
    for i in range(m):
        funcao2.append(1)
        
    print("MATRIZ A:",A)
    print("\nfuncao",funcao2)
    
    B = []
    N = []
    indexacaoB = []
    indexacaoN = []
    artificial = []
    for j in range(len(A1)):
        t = []
        for i in range(n1,n1 + m):
            if j ==0:
                indexacaoB.append(i)
                artificial.append(i)
            t.append(A1[j][i])
        B.append(t)
        l = []
        for i in range(n1):
            if j == 0:
                indexacaoN.append(i)
            l.append(A1[j][i])
        N.append(l)
         
    
    print("\nBASE:",B)
    
    print("\nNAO-BASE:",N)

    print("\nindexacaoB:",indexacaoB)

    print("\nindexacaoN:",indexacaoN)
    
    condicao = True
    iteracao = 1
    while(condicao):
        #passo1
        Xb = passo1(B,b)
        #passo 2
        CNK,indexCNK = passo2(B,funcao2,indexacaoN,N,indexacaoB)
        # passo 3
        if(CNK >= 0):
            condicao = False
            for i in artificial:
                print("ARITFICIAL DO IF::",artificial)
                print("B.count(i)",indexacaoB.count(i))
                if(indexacaoB.count(i) == 0):
                    print("\nBASEEE encontrada!!!!")
                    continue
                else:
                    print("\nVariavel artifical na base encontrada PROBLEMA INFACTIVEL!!!!")
                    indexacaoB = []
                    indexacaoN = []
                    B = []
                    N = []
                    return indexacaoB,indexacaoN,B,N
            break
        #passo 4
        y = passo4(N,indexCNK,B) 
        #passo5
        ver = sum(1 for number in y if number < 0)
        if(ver == len(y)):
            print("\nProblema Original Infactıvel")
            condicao = False
            break
        else:
            e,indexY = passo5(y,Xb)
        #passo6
        passo6(indexY,indexacaoB,indexacaoN,indexCNK,B,N,A1)
        print("\niteracao:::",iteracao)
        iteracao+=1
    print("ARIFICIAL",artificial)
    print("INDEXACAONN",indexacaoN )

    for j in range(len(N)):
        indexauxx = indexacaoN.copy()
        for i in artificial:
            print("COLUNA::",N[j])
            print("INDEXACAON::",indexauxx.index(i))
            print("INDEXACAON:::",indexauxx)
            N[j].remove(N[j][indexauxx.index(i)])
            indexauxx.remove(i)
        
    for i in artificial:
        indexacaoN.remove(i)
    print("INDEXACAOBB",indexacaoB)
    print("INDEXACAONN",indexacaoN )
    print("BBBBB",B)
    print("NNNN",N)
    return indexacaoB,indexacaoN,B,N
    
def main():
    while(1):
        tipo = input("O problema e de max ou min ?:")
        if(tipo=='max' or tipo=='min'):
            break

    n = int(input("quantidade de variaveis?:"))
    print()
    funcao = []
    print("entre com a função")
    for i in range(n):
        funcao.append(float(input()))
 
    if tipo == "min":
        teste = sum(1 for number in funcao if number > 0)
        if(teste == len(funcao)):
            print("\n PROBLEMA DE MIN COM TODAS AS VARIAVEIS POSTIVAS")
            exit(0)
    if tipo == "max":
        funcao = [x * (-1) for x in funcao]
    
    print("valores da funcao",funcao)
    
    m = int(input("quantidade de restrições?:"))

    A = []
    sinais = []
    print("sinais das restrições \n")
    for i in range(m):
        sinais.append(input())
    print("sinais da restrição",sinais)
    
    print("entre com os valores de A:\n")
    for i in range(m):
        a= []
        for j in range(n):
            a.append(float(input()))
        A.append(a)
        
    print("matriz A:\n",A)

    b = []
    print("entre com os valores de b:\n")
    for j in range(m):
        b.append(int(input()))
    
    print("vetor b:\n",b)

    for j in range(m):
        if b[j] < 0:
            A[j] = [x * (-1) for x in A[j]]
            b[j] = b[j] * (-1)
            if sinais[j] =="<=":
                sinais[j] =">="
            elif sinais[j] ==">=":
                sinais[j] ="<="
            elif sinais[j] ==">":
                sinais[j] ="<"
            elif sinais[j] =="<":
                sinais[j] =">"
            
    print("matriz expandida:\n",A,sinais,b)
    flag = 0
    A,N_variaveis_adicionais,flag = formaPadrao(A,sinais,b,m,funcao)
    
    if(flag):
        indexacaoB,indexacaoN,B,N = primeira_fase(A,sinais,b,funcao,m,n,N_variaveis_adicionais,tipo)
        print("INDEXACAOOOB",indexacaoB)
        print("indeaxaoN",indexacaoN)
        print("BBBBBBBBBBBB",B)
        print("NNNNNNNN",N)
        if(len(B)>0):
         segunda_fase(A,sinais,b,funcao,m,n,N_variaveis_adicionais,tipo,indexacaoB,indexacaoN,B,N)
    else:
       segunda_fase(A,sinais,b,funcao,m,n,N_variaveis_adicionais,tipo)    
    
main()
