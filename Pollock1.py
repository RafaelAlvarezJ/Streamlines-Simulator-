# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 21:50:32 2021

@author: rafa_
"""
#Programa para el ejercicio 3.30 del libro Datta-Gupta

import numpy as np
import matplotlib.pyplot as plt
import AuxiliaresPollock as A
import DatosPollock as d
import pandas as pd
import Units as u

"Aquí se crean los archivos de texto"
f = open("datos.txt","w")
f.write("Streamline,CeldaX,CeldaY,Tau,Dtau,TauCentr\n")
f.close()

g = open("datos2.txt","w")
g.write(' 0.8 \n ')
g.close()

"""
CÁLCULO DE PUNTOS INICIALES DE LAS STREAMLINES
"""
Xg,Yg = A.InitSL(d.DX, d.DY, d.NSL)

"""
INICIO DEL ALGORITMO DE POLLOCK
"""
Tiempo=np.arange(0,d.TiempoTotal,d.dt)
fw_avg=[] # Lista para la gráfica de fw
TAU=[] # Lista donde se van guardando los diversos dtau de la streamline en cuestión
ejeXgraficoH=[] # Lista creada para numerar las streamlines en el gráfico de tiempo de vuelo


for k in range (d.NSL): # Para todas las  SL
    if (k != 0):
        del csp, cea
        # fw_avg.append(fw[-1])
        # if len(fw_avg)==NSL:
        #     fw_avg2.append(st.mean(fw_avg))
        #     fw_avg=[]
        
    s=0 #Contador para los loops del ciclo para cada SL
    i=d.xcell_inj; j=d.ycell_inj #es la celda de partida de cada SL (0,0) 
    
    X=[]; Y=[]; DTAU=[]  # En cada loop se reinicia el arreglo (o lista)
    aux1=Xg[k]; aux2=Yg[k] #Punto inicial de la SL en cuestión
    
    X.append(aux1) #Se inicia la trayectoria de la streamline
    Y.append(aux2) #Se inicia la trayectoria de la streamline

    # Salen muchas SL pero mal trazadas    
    #while (i != d.xcell_prod) & (j != d.ycell_prod): # Mientras que no se alcance la celda de producción
    
    # Solo aparecen las primeras SL convergiendo bien 
    #while (i*j != (d.nx_cell-1)*(d.ny_cell-2) ): 
        
    # Hola crayola    
    while (i*j != 72):    

        # Condicional para asignar cara de entrada, al inicio del algoritmo para una SL
        if s==0:
            if   X[-1]==d.DX:
                csp=2
            elif Y[-1]==d.DY:
                csp=4
                
        cea,i,j =A.cara(csp,i,j) #aquí se obtiene la cara de entrada actual (cea) de la partícula
        
        
        aux1,aux2,aux3,csp = A.algoritmoPollock(X[-1], Y[-1], 
                                                d.U_xl[i][j], d.U_xr[i][j], 
                                                d.U_yb[i][j], d.U_yt[i][j], 
                                                cea, i,j,d.DX, d.DY, d.phi)
        
        
        # "Intento (1) para corregir "
        # if (X[-1] == (d.LX-d.DX)) & (Y[-1] >= (d.LY-d.DY)): # Llega por la cara izquierda
        #     break
        
        # if (X[-1] >= (d.LX-d.DX)) & (Y[-1] == (d.LY-d.DY)): # Llega por la cara inferior
        #     break
        # "Intento (1) para corregir"
        
        X.append(aux1)
        Y.append(aux2)
        DTAU.append(aux3)
        
        if (i == 9)&(j == 10):
            break
        if (i == 10)&(j == 9):
            break
        
        s=s+1
        if s>20: #Máximo de loops por SL 
            break
        
        AuxiliarTao=sum(DTAU) # Tiempo de vuelo hasta la cara
        #AuxiliarTao=AuxiliarTao-aux3/2 # Tiempo de vuelo "en el centro de la celda"
        
        f = open("datos.txt","a")
        #f.write(f' \n Cumulative Time : {"{:.3g}".format(time.Cum+time.dt)}  Seg  \n  =============================   \n ') 
        f.write(f' {k} , {i} , {j} , {AuxiliarTao} , {aux3},{AuxiliarTao-aux3/2}\n')
        f.close()
    
    TiempoVuelo=sum(DTAU)

    
    TAU.append(TiempoVuelo)
    ejeXgraficoH.append(k+1)
    X=np.array(X)
    Y=np.array(Y)
    #plt.plot(X,Y,color="black",marker="o",linestyle="dotted")
    plt.plot(X,Y,marker="o",linestyle="dotted")
    
    "Aquí se hacen los cálculos de transporte por streamline"
    DTAU=np.array(DTAU) # La lista se hace arreglo para los cálculos de transporte
    fw=np.zeros(len(DTAU))
    Sw= np.ones(len(DTAU))*d.Sw_ini
    A.transporte(Sw, fw, d.Srw, d.Sor, d.Muo, d.Muw, DTAU, d.dt, d.TiempoTotal, k, Tiempo, fw_avg)
    
    # Para el flujo fraccional
    

"Post proceso"

# Se crea el DataFrame con los datos de las celdas: i, j, tau, dtau, Sw
DATOS=pd.read_csv("datos.txt")
DATOS['Sw']=pd.read_csv("datos2.txt")

DATOS["Sw*Dtau"]=DATOS["Sw"]*DATOS["Dtau"]

# Se crean los arrays para los gráficos de contorno
graftau = np.zeros([d.nx_cell, d.ny_cell])
grafSw  = np.zeros([d.nx_cell, d.ny_cell])

"INICIO MODIFICACIÓN"
for i in range (d.nx_cell):
    for j in range (d.ny_cell):
        
        if (i==0) and (j==0): # en el pozo inyector
            graftau[i][j] = 0
            grafSw[i][j]  = (1-d.Sor) 
            
        elif (i == d.nx_cell-1) and (j == d.ny_cell-1): # en el pozo productor
            
            # Para el tiempo de vuelo
            "Opción 1"
            # array_TAU=np.array(TAU)
            # graftau[i][j] = array_TAU.mean()
            
            "Opción 2"
            graftau[i][j] = DATOS["Tau"].max()
            
            # Para la saturación
            grafSw[i][j] = DATOS["Sw"].min()
        
        else:
            
            df_aux = DATOS[ ( DATOS["CeldaX"].isin([i]) ) & ( DATOS["CeldaY"].isin([j]) ) ]
            
            # Para el tiempo de vuelo
            graftau[i][j] = df_aux["TauCentr"].mean()
            
            # Para la saturación
            numerador = df_aux["Sw*Dtau"].sum()
            denominador = df_aux["Dtau"].sum()
            grafSw[i][j] = numerador/denominador

"FIN MODIFICACIÓN"

fw_avg2=A.promedio_fw(fw_avg,Tiempo,d.NSL)

"Gráfica de Streamlines"
plt.figure(1)    
    
    #"Se escriben los datos de saturación"
    

#plt.title("$k_x=k_y$")
plt.xlabel("$x, [m]$")
plt.ylabel("$y, [m]$")
plt.grid()

"Gráfica de Tiempo de Vuelo"
# plt.figure(2)
# plt.barh(ejeXgraficoH,TAU,color="red",edgecolor="black")
# plt.ylabel("$Streamline$")
# plt.xlabel("Tiempo de Vuelo, $[s]$")
# plt.grid()

"Gráfica de contorno 1"
holax=np.linspace(0, d.LX, d.nx_cell)
holay=np.linspace(0, d.LY, d.ny_cell)
xg,yg = np.meshgrid(holax,holay) # Generación de malla para graficar

plt.figure(3)
plt.contourf(xg,yg,graftau/u.day)
plt.colorbar()
plt.title("Tiempo de Vuelo, $[d]$")
plt.xlabel("$x, [m]$")
plt.ylabel("$y, [m]$")
#plt.grid()

"Gráfica de contorno 2"
plt.figure(4)
plt.contourf(xg,yg,grafSw)
plt.colorbar()
plt.title(f'$S_w @ {d.TiempoTotal/u.day} [d]$')
plt.xlabel("$x, [m]$")
plt.ylabel("$y, [m]$")
#plt.grid()

"Gráfica de fw"
plt.figure(5)
Tiempo=np.array(Tiempo)
plt.plot(Tiempo/u.day, fw_avg2)
plt.title("Flujo fraccional, $[1/1]$")
plt.xlabel("$Tiempo, [d]$")
plt.ylabel("$f_w, [1]$")
plt.grid()