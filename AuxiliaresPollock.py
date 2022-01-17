# -*- coding: utf-8 -*-
"""
Created on Sun Oct 17 15:29:28 2021

@author: rafa_

Archivo en el cual se encuentran todas las funciones auxiliares para la 
simulación numérica de Streamlines. De esta manera, el código principal no es 
tan extenso."""

import numpy as np
import statistics as st
from numba import jit

@jit(nopython=True)
def InitSL(dx,dy,n):
    # Función que nos da los puntos iniciales de las SL, dando como input 
    # el número de SL deseado (debe ser par) y las dimensiones de la celda
    
    # pensada para un modelo five-spot, donde de la celda del pozo inyector, 
    # solo dos de sus caras llevarán streamlines
    
    # - Se asume una celda cuadrada
    
    # n : número de streamlines (tiene que ser un número par)
    
    Xg=np.zeros(n) 
    Yg=np.zeros(n)
    
    aux1 = (n/2) # para los puntos de partida de las SL
    aux2 = (dx/aux1) # para dividir una cara en (aux1) partes iguales
    aux3 = aux2/2
    
    for i in range (n):
        if i < aux1: # para la primera mitad
            Xg[i] = dx
            Yg[i] = aux2*(i)+aux3
        else:        # para la segunda mitad
            Xg[i] = aux2*(i)+aux3-dx
            Yg[i] = dy
    
    return (Xg,Yg)

@jit(nopython=True)
def Press(Po_old, i, j, nodosx, nodosy):
    # Función que tiene como objetivo acortar el código principal. Se entregan 
    # las Sw E,W,N,S considerando las condiciones de frontera (cerradas).
    
    # Función análoga a "Sw"
    
    # ESQUINA SUROESTE
    if (i==0) and (j==0):
        P_E = Po_old[i+1][j]
        P_W = Po_old[i][j]
        P_N = Po_old[i][j+1]
        P_S = Po_old[i][j]
        
    # ESQUINA NOROESTE
    elif (i==0) and (j==nodosy-1):
        P_E = Po_old[i+1][j]
        P_W = Po_old[i][j]
        P_N = Po_old[i][j]
        P_S = Po_old[i][j-1]
        
    # ESQUINA SURESTE
    elif (i==nodosx-1) and (j==0):
        P_E = Po_old[i][j]
        P_W = Po_old[i-1][j]
        P_N = Po_old[i][j+1]
        P_S = Po_old[i][j]
 
    # ESQUINA NORESTE
    elif (i==nodosx-1) and (j==nodosy-1):
        P_E = Po_old[i][j]
        P_W = Po_old[i-1][j]
        P_N = Po_old[i][j]
        P_S = Po_old[i][j-1]
        
    # FRONTERA OESTE
    elif (i==0):
        P_E = Po_old[i+1][j]
        P_W = Po_old[i][j]
        P_N = Po_old[i][j+1]
        P_S = Po_old[i][j-1]
        
    # FRONTERA ESTE
    elif (i==nodosx-1):
        P_E = Po_old[i][j]
        P_W = Po_old[i-1][j]
        P_N = Po_old[i][j+1]
        P_S = Po_old[i][j-1]
        
    # FRONTERA SUR
    elif (j==0):
        P_E = Po_old[i+1][j]
        P_W = Po_old[i-1][j]
        P_N = Po_old[i][j+1]
        P_S = Po_old[i][j]
        
    # FRONTERA NORTE
    elif (j==nodosy-1):
        P_E = Po_old[i+1][j]
        P_W = Po_old[i-1][j]
        P_N = Po_old[i][j]
        P_S = Po_old[i][j-1]
        
    else:
        P_E = Po_old[i+1][j]
        P_W = Po_old[i-1][j]
        P_N = Po_old[i][j+1]
        P_S = Po_old[i][j-1]

    return (P_E, P_W, P_N, P_S)


def pseudotiempo(c,Dtau):
    "  Función que se ocupa para el calculo de la posición de la partícula "
    if c==0:
        y=Dtau
    else:
        y=(np.exp(c*Dtau)-1)/c
    return (y)

@jit(nopython=True)
def velocidad(vp1,cp,p1,p):
    # vp1=velocidad a la iquierda o al fondo de la celda
    # cp=coeficiente
    # p1=coordenada a la izquierda o al fondo de la celda
    # p=punto de partida, (x/y)
    v=vp1+cp*(p-p1)
    return (v)

@jit(nopython=True)
def cara(csp,i,j):
    """ Función que devuelve la cara entrada actual (cea) de la partícula, también se
    encarga de modificar los índices de celda (i,j) para el siguiente cálculo """ 
    
    # 1 = cara oeste
    # 2 = cara este
    # 3 = cara sur
    # 4 = cara norte
    
    #la cara de salida pasada (csp) es la cara de entrada actual (cea)
    if   csp==1: #Sale al oeste
        cea=2
        i=i-1
        
    elif csp==2: #Sale al este
        cea=1
        i=i+1
        
    elif csp==3: #Sale al sur
        cea=4
        j=j-1
        
    elif csp==4: #Sale al norte
        cea=3
        j=j+1
        
    return (cea,i,j)

#@jit(nopython=True)
def algoritmoPollock(x0,y0,ux1,ux2,uy1,uy2,cea,i,j,DX,DY,phi):
    #Cálculo de los coeficientes
    cx=(ux2-ux1)/DX
    cy=(uy2-uy1)/DY
    
    #print(f"i={i}, j={j}")
    #print(f"cx={cx}, cy={cy}, suma={cx+cy}")
    
    
    # cea (cara de entrada actual)
    # 1 = cara oeste
    # 2 = cara este
    # 3 = cara sur
    # 4 = cara norte
    """
    INICIO MODIFICACIÓN
    """
    "La velocidad en x,y es constante."
    if cx==0 and cy==0: 
        if cea==1:   # Cara oeste
            ux0=ux1
            uy0=uy1#velocidad(uy1, cy, j*DY, y0)
            
        elif cea==2: # Cara este
            ux0=ux2
            uy0=uy1#velocidad(uy1, cy, j*DY, y0)
        
        elif cea==3: # Cara sur
            uy0=uy1
            ux0=ux1#velocidad(ux1, cx, i*DX, x0)
            
        elif cea==4: # Cara norte
            uy0=uy2
            ux0=ux1#velocidad(ux1, cx, i*DX, x0)
            
            
        Dtau_x1=(DX*i    -x0)/ux0    
        Dtau_x2=(DX*(i+1)-x0)/ux0
        Dtau_y1=(DY*j    -y0)/uy0
        Dtau_y2=(DY*(j+1)-y0)/uy0
        
        
    
    else: #Si cx,cy son diferentes de 0
            # Cara de entrada actual
        if cea==1:   # Cara oeste
            ux0=ux1
            uy0=velocidad(uy1, cy, j*DY, y0)
            
        elif cea==2: # Cara este
            ux0=ux2
            uy0=velocidad(uy1, cy, j*DY, y0)
        
        elif cea==3: # Cara sur
            uy0=uy1
            ux0=velocidad(ux1, cx, i*DX, x0)
            
        elif cea==4: # Cara norte
            uy0=uy2
            ux0=velocidad(ux1, cx, i*DX, x0)
        
        "Cálculo de tiempo de vuelo a cada cara"
        # Para velocidades en x
        if ux0 == 0: # Para evitar la división entre cero
            Dtau_x1 = np.NaN
            Dtau_x2 = np.NaN
        else:
            if ux1 == 0: # Para evitar el log(cero)
                Dtau_x1 = np.NaN
            else:
                Dtau_x1=(1/cx)*np.log(ux1/ux0)
            
            if ux2 == 0: # Para evitar el log(cero)
                Dtau_x2 = np.NaN
            else:
                Dtau_x2=(1/cx)*np.log(ux2/ux0)        
        
        # Para velocidades en y
        if uy0 == 0: # Para evitar la división entre cero
            Dtau_y1 = np.NaN
            Dtau_y2 = np.NaN
        else:
            if uy1 == 0: # Para evitar el log(cero)
                Dtau_y1 = np.NaN
            else:
                Dtau_y1=(1/cy)*np.log(uy1/uy0)
            
            if uy2 == 0: # Para evitar el log(cero)
                Dtau_y2 = np.NaN 
            else:
                Dtau_y2=(1/cy)*np.log(uy2/uy0)
        
    "Aquí se sale de los condicionales de cx,cy"    
    #print(f"ux1={ux1}, ux2={ux2}, uy1={uy1}, uy2={uy2} ")

    Dtau_array=np.array([Dtau_x1,Dtau_x2,Dtau_y1,Dtau_y2])
    #print(f"{Dtau_array}")
    

    """
    FIN MODIFICACIÓN
    """

    for i in range (len(Dtau_array)):
        if Dtau_array[i]<=1:
            Dtau_array[i]=np.NaN
            
    
    Dtau=np.nanmin(Dtau_array)#*phi
    
    #Para devolver la cara de salida 
    if Dtau==Dtau_array[0]:
        csp=1
    elif Dtau==Dtau_array[1]:
        csp=2
    elif Dtau==Dtau_array[2]:
        csp=3
    elif Dtau==Dtau_array[3]:
        csp=4
        
    x_new=x0+ux0*pseudotiempo(cx,Dtau)
    y_new=y0+uy0*pseudotiempo(cy,Dtau)
    #print(f"x_new={x_new}, y_new={y_new}")
    #print(f"cea={cea}, csp={csp} \n")
    
    return (x_new,y_new,Dtau*phi, csp)
#    return(Dtau)
#    return(Dtau_array)

def kr(Sw, Srw, Sro):
    "Función que calcula las permeabilidades relativas cuadráticas"
    Scf=(Sw-Srw)/(1-Srw-Sro)
    
    krw=Scf**2
    kro=(1-Scf)**2
    
    return(krw, kro)

def fractionalflow(kro,krw,Muo,Muw):
    "Función que se encarga de calcular el fluj fraccional de agua"
    fw=1/(1+(kro*Muw)/(krw*Muo))
    
    return(fw)
    
#@jit(nopython=True)
def transporte(Sw, fw, Srw, Sro, Muo, Muw, DTAU, dt, TiempoTotal, k, Tiempo, fw_avg):
    "Función que calcula la saturación explícitamente"
    # Sw, DTAU, fw son arreglos. Las demás variables son escalares.
    
    # tiempodevuelo=sum(DTAU)

    # DTAU2=np.linspace(0,tiempodevuelo,100) # Segunda discretización
    
    # "Arreglos de tamaño DTAU2"
    # Sw=np.ones(len(DTAU2))*Sw_ini
    # fw=np.zeros(len(DTAU2))
    
    t=0 
    "Ciclo temporal"
    while t <= TiempoTotal:
        
        for i in range (len(DTAU)): # A lo largo de la streamline
            krw, kro = kr(Sw[i], Srw, Sro) #kro, krw se actualizan en cada ciclo
            fw[i] = fractionalflow(kro, krw, Muo, Muw)

            if i==0: # Condición de inyección
                Sw[i] = (1-Sro)
                
                
            else: # Cálculo IMPES
                aux = (dt/DTAU[i])*(fw[i]-fw[i-1])
                
                Sw[i] = Sw[i]-aux
        
            if i==(len(DTAU)-1): # Para graficar los flujos fraccionales
                fw_avg.append(fw[i])
                
            # if Sw[i] < (1-Srw): # Si es menor a la Sw mínima
            #     Sw[i] = 1-Srw
                
            # elif Sw[i] > (1-Sro):
            #     Sw[i] = 1-Sro
        
        t = t+dt
        #print(t)
    
        fw_avg.append(fw[-1])

    g = open("datos2.txt","a")
    for i in range (len(DTAU)):
        g.write(f' {Sw[i]} \n ')
    g.close()
    
#@jit(nopython=True)
def promedio_fw(fw_avg,Tiempo,NSL):
    
    fw_avg=np.array(fw_avg) # La lista se convierte en arreglo
    fw_avg2=np.zeros(len(Tiempo)) # Se crea un arreglo para fw, que coincida en tamaño con el de tiempo
    
    aux=np.zeros(NSL) # aquí se guardan los fw de cada SL en el pozo. Se va reemplazando para cada dt
    aux2=len(Tiempo) 
    
    # Para sacar el promedio de fw de las SL en el pozo
    # Primero agrupamos las N streamlines
    
    for i in range(len(Tiempo)): # Para todo el arreglo de Tiempo
        for j in range (NSL): 
            aux[j]=fw_avg[i+j*aux2]
            
        aux3=st.mean(aux)
        fw_avg2[i]=aux3 # Se guarda el fw promedio.
        
    return (fw_avg2) 