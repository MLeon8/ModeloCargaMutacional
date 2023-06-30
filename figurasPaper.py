#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 11:46:52 2023

@author: jatec
"""

import pandas
import numpy as np
import matplotlib.pyplot as plt
import colorsys
import matplotlib.colors as mc
import matplotlib as mpl

# mpl.style.use(sty)
def plotOverride(df,title):
    plt.figure()
    plt.title(title)
    plt.grid(True)
    plt.plot(df['Step'][:],df['AntiCancer'][:],label="Sistema Inmune",color = "blue" )
    plt.xlabel("Iteration(days)")
    plt.ylabel("Number of cells")
    plt.plot(df['Step'][:],df['ProCancer'][:],label="Cancer", color = "red")
    plt.legend()
    plt.savefig(f'{title}.pdf', format="pdf", bbox_inches="tight")

### Código para código que empieza inicializando todo a 0
df = pandas.read_csv('/home/jatec/Desktop/InmunoedicionDelCancer-Ising/DatosSimulaciones/Figuras Paper/Escape/model_data[Thu Jun  8 23_54_04 2023].csv')
plt.figure()
title = "Cancer: Medio, IS: Medio"
plt.title(title)
# dataM = np.matrix(df)
# promedioCC = np.zeros([100,df.shape[0]])
# promedioCC = np.array(df['AntiCancer'][len(df)-101+1:])
# promedioIS = np.array(df['ProCancer'][len(df)-101+1:])
for i in range(0,len(df)-101,101):
    # promedioCC = np.insert(promedioCC,0,np.array(df['AntiCancer'][i+1:i+101]),axis=1)
    # promedioIS = np.append(promedioIS,np.array(df['ProCancer'][i+1:i+101]),axis=1)
    plt.plot(df['Step'][i+1:i+101],df['AntiCancer'][i+1:i+101],color = "blue" )
    plt.plot(df['Step'][i+1:i+101],df['ProCancer'][i+1:i+101], color = "red")
plt.plot(df['Step'][len(df)-101 + 1:],df['AntiCancer'][len(df)-101+1:],label="Sistema Inmune",color = "blue" )
plt.plot(df['Step'][len(df)-101 + 1:],df['ProCancer'][len(df)-101+1:],label="Cancer", color = "red")
# promedioCC.append(df['AntiCancer'][len(df)-101 + 1:])
# promedioIS.append(df['ProCancer'][len(df)-101 + 1:])
df_filtered = df[df['Step'] != 0]

promedio_pc = df_filtered.groupby('Step')['ProCancer'].mean()
promedio_ac = df_filtered.groupby('Step')['AntiCancer'].mean()

maximos = df_filtered.groupby('Step')[['ProCancer', 'AntiCancer']].max()
minimos = df_filtered.groupby('Step')[['ProCancer', 'AntiCancer']].min()

plt.grid(True)
plt.xlabel("Iteration(days)")
plt.ylabel("Number of cells")
plt.legend()
plt.savefig(f'{title}.pdf', format="pdf", bbox_inches="tight")


plt.figure(figsize=(10, 5))
plt.plot(promedio_pc.index, promedio_pc, label='Promedio ProCancer', color='blue')
plt.plot(promedio_ac.index, promedio_ac, label='Promedio Sistema Inmune', color='red')
plt.fill_between(maximos.index, maximos['ProCancer'], minimos['ProCancer'], alpha=0.3, color='blue')
plt.fill_between(maximos.index, maximos['AntiCancer'], minimos['AntiCancer'], alpha=0.3, color='red')
plt.xlabel('Tiempo')
plt.ylabel('Valor')
plt.grid(True)
plt.title('Promedio de las poblaciones y su variabilidad')
plt.legend()
plt.show()