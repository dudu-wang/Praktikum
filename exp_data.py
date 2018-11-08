import matplotlib.pyplot as plt
import numpy as np

table_85 = np.loadtxt('Rb-85.txt', dtype=float)
table_85 = np.transpose(table_85)

names = ['F=2 - F=1', 'F=3 - F=2', 'F=4 - F=3']
         
fig, axs = plt.subplots()

x_axis = np.arange(3)

lit_werte = table_85[1][0:3]
exp_werte1 = table_85[2][0:3]
yerr1 = table_85[3][0:3]
exp_werte2 = table_85[2][3]
yerr2 = table_85[3][3]

axs.scatter(x_axis, lit_werte, label='Literaturwerte')
axs.errorbar(x_axis, exp_werte1, yerr = yerr1, fmt = 'ro', barsabove = True,
             label='Experimentelle Werte')
axs.errorbar(x_axis[1], exp_werte2, yerr = yerr2, fmt = 'ro', barsabove = True)
axs.legend(loc=4)
axs.set_ylabel('Frequenz (in MHz)')
plt.xticks(x_axis, names)
plt.show()
