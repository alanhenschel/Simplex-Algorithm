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


def pivotamento(A):
    for i in range(len(A) - 1):
        for j in range(i + 1,len(A)):
            n1 = A[j][i]
            if n1 == 0:
                for k in range(i,len(A)):
                    if A[k][i] != 0:
                        swap(A, i, k)
                        n1 = A[j][i]
            n2 = A[i][i]
            if n2 == 0:
                continue
            div = n1 / n2
            A[j] = sumVector(multVector(A[i], -div), A[j])
    return A


def multVector(v, x):
    newV = []
    for i in range(len(v)):
        newV += [v[i] * x]
    return newV

def sumVector(v1, v2):
    newV = []
    for i in range(len(v1)):
        newV += [v1[i] + v2[i]]
    return newV

def swap(A, i, j):
    newI = []
    newJ = []
    for index in range(len(A) + 1):
        newI += [A[j][index]]
        newJ += [A[i][index]]
    A[i] = newI
    A[j] = newJ

def printMatrix(A):
    for i in A:
        for j in i:
            print ('{}'.format(j), end=" ")
        print("\n")

def triangular(A):
    x = [0 for i in range(len(A))]
    for i in range(len(A) - 1, -1, -1):
        print(A[i][i])
        x[i] = A[i][len(A)]/A[i][i]
        for j in range(i-1, -1, -1):
            A[j][len(A)] -= A[j][i] * x[i]
    return x
    

def gauss(A,b):
    
    gauss = copy.deepcopy(A)
    for i in range(len(gauss)):
        gauss[i].append(b[i])
    gauss = gausss(gauss)  
    #gauss = pivotamento(gauss)
    #gauss = triangular(gauss)
    return gauss

def gauss2(A):
    
    gauss = copy.deepcopy(A)     
    #gauss = pivotamento(gauss)
    #gauss = triangular(gauss)
    gauss = gausss(gauss)
    return gauss
   
def formaPadrao(A,sinais,b,m,funcao):
    N_variaveis_adicionais = len(sinais)
    flag = 0
    for j in range(m):
        for i in range(N_variaveis_adicionais):
            if sinais[j]=="<=" and j==i:
                A[j].append(1.0)
            elif sinais[j]==">=" and j==i:
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
    for i in indexacaoB:
        func.append(funcao[i])
    for i in range(len(Xb)):
        resultado += Xb[i]*func[i]
    
    if tipo=="max":
        resultado = resultado * (-1)
    print("\na FUDENDO SOLUCAO EHHHH:",resultado)

def segunda_fase(A,sinais,b,funcao,m,n,N_variaveis_adicionais,tipo,indexacaoB = [],indexacaoN = [],B = [],N = []):
    print("\n\nENTRANDO NA SEGUDNA FASE IIIIIIIIII\n\n")
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
        Xb = gauss(B,b)
        print("\nXB:",Xb)
        #passo 2
        transpose = transposta(B)
    
        print("tamanho transposta",len(transpose))
        print("tamanho funcao",len(funcao))
        multiplicador = []
        for i in range(len(transpose)):
            multiplicador.append(list(np.append(transpose[i],funcao[indexacaoB[i]])))
    
        print("\nvetor multiplicador simplex:",multiplicador)
        multiplicador2 = gauss2(multiplicador)
        print("\nvetor multiplicador222GAUSS simplex:",multiplicador2)
    
        #2.2
        CN = []
        print("INDEXACAOOONN:",indexacaoN)
        print("FUNCAOOOO",funcao)
        print("NNNNNN",N)
        for j in range(len(N[0])):
            c = 0
            for i in range(len(multiplicador2)):
                c +=multiplicador2[i] * N[i][j]
    
            CN.append(funcao[indexacaoN[j]] - c)
        
        #2.3
    
        CNK = min(CN)
        indexCNK = CN.index(min(CN))
        print("CUSTOS",CN)
        print("VARIAVEL A ENTRAR NA BASE",indexCNK)
    
        # passo 3
    
        if(CNK >= 0):
            print("\nSolucao Otima encontrada!!!!")
            condicao = False
            solucaoOtima(Xb,indexacaoB,funcao,tipo)
            break
    
        #passo 4 
        colunaN = column(N,indexCNK)
        
        y = gauss(B,colunaN)
    
        print("\nYYYYYYYYYYY",y)
    
        #passo5
        ver = sum(1 for number in y if number < 0)
        if(ver == len(y)):
            print("problema nao tem solucao otima finita")
            condicao = False
            break
        else:
            e = float("inf")
            indexY = 0         
            for i in range(len(Xb)):
                if y[i]> 0 and Xb[i]/y[i] < e:
                    e = Xb[i]/y[i]
                    indexY = i
    
    
            print("\nEEEEEEEEE:",e)
            print("\nEEEEEEEEEYYYYYY:",indexY)
    
        
            #passo6
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
        Xb = gauss(B,b)
        print("\nXB:",Xb)
        #passo 2
        transpose = transposta(B)
    
        print("tamanho transposta",len(transpose))
        print("tamanho funcao",len(funcao))
        multiplicador = []
        for i in range(len(transpose)):
            multiplicador.append(list(np.append(transpose[i],funcao2[indexacaoB[i]])))
    
        print("\nvetor multiplicador simplex:",multiplicador)
        multiplicador2 = gauss2(multiplicador)
        print("\nvetor multiplicador222GAUSS simplex:",multiplicador2)
    
        #2.2
        CN = []
        for j in range(len(N[0])):
            c = 0
            for i in range(len(multiplicador2)):
                c +=multiplicador2[i] * N[i][j]
    
            CN.append(funcao2[indexacaoN[j]] - c)
        
        #2.3
    
        CNK = min(CN)
        indexCNK = CN.index(min(CN))
        print("CUSTOS",CN)
        print("VARIAVEL A ENTRAR NA BASE",indexCNK)
    
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
            
    
        #passo 4 
        colunaN = column(N,indexCNK)
        
        y = gauss(B,colunaN)
    
        print("\nYYYYYYYYYYY",y)
    
        #passo5
        ver = sum(1 for number in y if number < 0)
        if(ver == len(y)):
            print("\nProblema Original Infactıvel")
            condicao = False
            break
        else:
            e = float("inf")
            indexY = 0         
            for i in range(len(Xb)):
                if y[i]> 0 and Xb[i]/y[i] < e:
                    e = Xb[i]/y[i]
                    indexY = i
    
    
            print("\nEEEEEEEEE:",e)
            print("\nEEEEEEEEEYYYYYY:",indexY)
    
        
            #passo6
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
                auxx2.append(column(A1,i))
                
            for i in range(len(auxx2)):
                for j in range(len(auxx2[0])):
                    B[j][i] = auxx2[i][j]
            
            auxx2 = []
            for i in indexacaoN:
                auxx2.append(column(A1,i))
                
            for i in range(len(auxx2)):
                for j in range(len(auxx2[0])):
                    N[j][i] = auxx2[i][j]
                    
            print("\nMATRIZ AAAAAA",A1)
            print("\nBASE:",B)
            print("\nNAO-BASE:",N)
            print("\nindexacaoB:",indexacaoB)
            print("\nindexacaoN:",indexacaoN)
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
    tipo = input("O problema e de max ou min ?:")

    n = int(input("quantidade de variaveis?:"))
    print()
    funcao = []
    print("entre com a função")
    for i in range(n):
        funcao.append(float(input()))
 
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
