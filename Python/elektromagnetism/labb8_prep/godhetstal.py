import numpy as np
import matplotlib.pyplot as plt

h = 0.01

x = np.arange(0.1, 3, h)

def f(Q_0):
    return 1/(np.sqrt(1 + ((x - 1/x)*Q_0)**2))

Q_0s = np.array([1, 10, 100])

fig, ax = plt.subplots()

for i in range(Q_0s.size):
    ax.plot(x, f(Q_0s[i]), label=f"Q_0 = {Q_0s[i]}")

ax.legend()

plt.show()

# Resultat: skarp resonans <=> Q_0 stor 
# <=> R liten / C liten / w_0 liten <=> R liten / C liten / L stor 
