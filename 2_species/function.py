# %%

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

u = np.array([0, 0])
x1 = np.array([50, 30])
x2 = np.array([200, 80])
g = np.array([0.025, 0.025, 0.025]) # growth rate
m = np.array([0.017, 0.017, 0.001]) # mortality rate
b = np.array([0.75, 0.75, 0.90]) # birth rate
cg = np.array([0.0067, 0.0067, 0.0067]) # competition growth
cb = np.array([0.0125, 0.0125, 0.0125]) # competition birth
cm = np.array([0.0008, 0.0008, 0.0008]) # competition mortality
g1 = np.array([0.013, 0.013])
g2 = np.array([0.16, 0.16])
v1 = np.array([0.066, 0.066, 0.066])
v2 = np.array([2.29, 2.29, 2.29])

def forest_new_list(x1, x2, u):
    x1_new = []
    x2_new = []
    Mod_2 = np.sum(x2 * g2)
    Mod_1 = np.sum(x1 * g1)
    for i in range(len(x1)):
        x1_new.append(max(0, x1[i] -g[i] * x1[i] * (1 - cg[i] * Mod_2) + b[i] * g2[i] * x2[i] * (1 - cb[i] * (Mod_1 + Mod_2)) - x1[i] * (cm[i] * Mod_2 + m[i])))
        x2_new.append(max(0, x2[i] + g[i] * x1[i] * (1 - cg[i] * Mod_2) - m[i] * x2[i] - u[i] * x2[i]))
    return [x1_new, x2_new]

def forest_simul_list(T, x1, x2, u):
    forest = pd.DataFrame({'x1': x1, 'x2': x2, 'Prod': 0, 'u': 0, 't': 1, 'sp' : np.arange(1, len(x1) + 1)})
    for t in range(2, T + 1):
        h = [0, 0]
        if t % 5 == 1:
            h = u[(t + 4) // 5 - 1]
        new_forest = forest_new_list(forest['x1'].tail(len(x1)).values, forest['x2'].tail(len(x2)).values, h)
        new_forest = pd.DataFrame({'x1': new_forest[0], 'x2': new_forest[1], 'Prod': 0, 'u': 0, 't': t, 'sp': np.arange(1, len(x1) + 1)})
        forest = pd.concat([forest, new_forest], ignore_index=True)
    forest['sp'] = forest['sp'].astype('category')
    return forest

def plot_population(forest) :
    fig, axs = plt.subplots(nrows=len(forest['sp'].unique()), sharex=True, figsize=(8, 6))

    for i, sp in enumerate(forest['sp'].unique()):
        df = forest[forest['sp'] == sp]
        axs[i].plot(df['t'], df['x1'], color='C0')
        axs[i].plot(df['t'], df['x2'], color='C1', linestyle='--')
        axs[i].set_ylabel(f'Species {sp}')
        axs[i].set_ylim(bottom=0)

    axs[-1].set_xlabel('Time')
    plt.show()

def forest_5(x1, x2, u):
    year = forest_new_list(x1, x2, u)
    for i in range(4):
        year = forest_new_list(year[0], year[1], u)
    return year

import itertools

f = 50
n = 2  # number of species
x1 = np.arange(0, 501, f)
x2 = np.arange(0, 501, f)
fc = 0.25
u = np.arange(0, 1+fc, fc)

# Créer le tableau de toutes les combinaisons possibles
table = list(itertools.product(x1, x1,x1,x1,u,u))
table = pd.DataFrame(table, columns=['x1_1', 'x2_1', 'x1_2', 'x2_2', 'u_1', 'u_2'])

# Faire des colonnes avec des listes genre colonne x1 avec dedans [x1_1, x1_2]
table['x1'] = table[['x1_1', 'x1_2']].values.tolist()
table['x2'] = table[['x2_1', 'x2_2']].values.tolist()
table['u'] = table[['u_1', 'u_2']].values.tolist()

# calculer les populations 5 ans après
table['x1_5'] = table.apply(lambda row: forest_5(row['x1'], row['x2'], row['u'])[0], axis=1)


import numpy as np

f = 50
x1 = np.arange(0, 501, f)
x2 = np.arange(0, 501, f)
fc = 50
u = np.arange(0, 201, fc)

M = np.zeros((len(x1), len(x1), len(x1), len(x1),len(x1), len(x1), len(x1), len(x1),len(u), len(u)))

# Now you can use M in your Python code.
