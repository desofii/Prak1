import matplotlib
import matplotlib.pyplot as plt
import numpy as np

x, y = np.loadtxt("output.txt", unpack=True)

plt.plot(x, y, label="график функции")
plt.legend()

plt.show()
