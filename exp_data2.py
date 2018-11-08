import matplotlib.pyplot as plt
import numpy as np

table_87 = np.loadtxt('Rb-87.txt', dtype=float)
table_87 = np.transpose(table_87)

names = ['F=1 - F=0', 'F=2 - F=1', 'F=3 - F=2']
         
fig, axs = plt.subplots()

x_axis = np.arange(3)

lit_werte = table_87[1][0:3]
exp_werte1 = table_87[2][0:3]
yerr1 = table_87[3][0:3]
exp_werte2 = table_87[2][3]
yerr2 = table_87[3][3]

axs.scatter(x_axis, lit_werte, label='Literaturwerte')
axs.errorbar(x_axis, exp_werte1, yerr = yerr1, fmt = 'ro', barsabove = True,
             label='Experimentelle Werte')
axs.errorbar(x_axis[1], exp_werte2, yerr = yerr2, fmt = 'ro', barsabove = True)
axs.legend(loc=4)
axs.set_ylabel('Frequenz (in MHz)')
plt.xticks(x_axis, names)
plt.show()
