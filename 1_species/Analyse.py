import numpy as np

def x1_new(x1, x2, u):
    return x1 + (g1 - cg * x2) * x1 - m1 * x1 - u[1]

def x2_new(x1, x2, u):
    return x2 + (g2 - cb * x1 - cm * x2) * x2 - m2 * x2 + b * x2 * (1 - np.exp(-v1 * x1 - v2 * x2)) - u[2]

def forest_5(x1, x2, u1, u2):
    Year = np.array([x1_new(x1, x2, u1), x2_new(x1, x2, u2)])
    for i in range(4):
        Year = np.array([x1_new(Year[0], Year[1], 0), x2_new(Year[0], Year[1], 0)])
    return Year

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
    x1_new = np.zeros(len(x1))
    x2_new = np.zeros(len(x2))
    for i in range(len(x1)):
        x1_new[i] = x1_new(x1[i], x2[i], u)
        x2_new[i] = x2_new(x1[i], x2[i], u)
    return [x1_new, x2_new]

print("Hello world!")