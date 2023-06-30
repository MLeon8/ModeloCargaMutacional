#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 11:16:32 2023

@author: jatec
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np

# Ruta del archivo XLS
ruta_archivo = "/home/jatec/Downloads/41586_2013_BFnature12213_MOESM256_ESM.xls"

# Pedir al usuario el valor del tumor_type
# tumor_type_usuario = input("Ingrese el valor de tumor_type: ")

# # Verificar si el tumor_type existe en el archivo
# if tumor_type_usuario not in datos['tumor_type'].values:
#     print(f"El tumor_type '{tumor_type_usuario}' no existe en el archivo.")
#     exit()

tumor_type_usuario = "Breast"
# "Glioblastoma multiforme"

# Leer la segunda página del archivo XLS
datos = pd.read_excel(ruta_archivo, sheet_name=1)

# Filtrar los datos por el valor del tumor_type dado por el usuario
datos_filtrados = datos[datos['tumor_type'] == tumor_type_usuario]

# Ordenar los datos por la columna 'n_coding_mutations' en orden ascendente
datos_ordenados = datos_filtrados.sort_values('n_coding_mutations')

# Obtener los valores de 'n_coding_mutations' y 'coding_mutation_rate'
n_coding_mutations = datos_ordenados['n_coding_mutations']
coding_mutation_rate = datos_ordenados['coding_mutation_rate']

# Configurar la escala logarítmica en el eje y
# plt.yscale('log')

# Graficar los datos
plt.plot(np.arange(len(n_coding_mutations)),n_coding_mutations, '.')

# Configurar etiquetas y título del gráfico
plt.xlabel('n')
plt.ylabel('n_coding_mutations')
plt.title(f'Datos para: {tumor_type_usuario}')

# Formatear los valores del eje y en notación estándar
ax = plt.gca()
ax.yaxis.set_major_formatter(ScalarFormatter())

# Activar el grid y los subgrid
plt.grid(True, which='both', linestyle='dotted')


# # Etiquetas para las separaciones de los subgrids
# yticks = ax.yaxis.get_major_ticks()
# for tick in yticks:
#     tick.label1.set_visible(False)  # Ocultar las etiquetas por defecto

# # Agregar etiquetas personalizadas a las separaciones de los subgrids
# for y_value in ax.yaxis.get_minorticklocs():
#     ax.text(0.5, y_value, f'{y_value:.2f}', color='gray', ha='left', va='center')


# Mostrar el gráfico
plt.show()
