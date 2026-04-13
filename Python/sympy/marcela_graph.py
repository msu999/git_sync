import csv
from matplotlib import pyplot as plt
import numpy as np
from numpy import genfromtxt

import pandas as pd

df = pd.read_csv("compression.csv")
data = df.to_numpy()

print(data)
print(data.shape)
dataT = data.T
F_1 = dataT[1][2:].astype(float)
stroke_1 = dataT[2][2:].astype(float)
A_1 = np.pi*(5.87/2)**2
stress_1 = F_1/A_1  # MPa
strain_1 = (stroke_1 - 11.91)/11.91 + 1  # Percento

F_2 = dataT[4][2:].astype(float)
stroke_2 = dataT[5][2:].astype(float)
A_2 = np.pi*(5.91/2)**2
stress_2 = F_2/A_2  # MPa
strain_2 = (stroke_2 - 12.32)/12.32 + 1  # Percento

print(F_2)
print(stress_1)

fig, ax = plt.subplots()

#ax.legend()
ax.set(xlabel="Strain [%]", ylabel="Stress [MPa]")

ax.plot(strain_1, stress_1, "b")
ax.plot(strain_2, stress_2, "y")

plt.show()
