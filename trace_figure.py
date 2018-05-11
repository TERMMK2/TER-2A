#coding : utf8

from numpy import *
from matplotlib.pyplot import *
from matplotlib import rc
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

z = loadtxt ('points.txt')
# z = loadtxt ('ro_thermo_couple_mat_variable')
#
#plot (z[:,0],z[:,2],"k", label = r' $\rho (t)$ \`a la paroi')
#plot (z[:,0],z[:,4],"r", label = r' $\rho (t)$ \`a 1 mm')
#plot (z[:,0],z[:,6],"b", label = r' $\rho (t)$ \`a 2 mm')
#plot (z[:,0],z[:,8],"g", label = r' $\rho (t)$ \`a 3 mm')
#plot (z[:,0],z[:,10],"c", label = r'$\rho (t)$ \`a 4 mm')
#plot (z[:,0],z[:,12],"y", label = r'$\rho (t)$ \`a 5 mm')


u = loadtxt ('ro_thermo_couple_mat-constant')
#v = loadtxt ('EC_PyrolyseMC_explicite.txt')

#plot (u[:,0],u[:,1],"k--", label = r'$\rho (t)$ \`a la paroi mod\`ele')
#plot (u[:,0],u[:,2],"r--", label = r'$\rho (t)$ \`a 1 mm mod\`ele')
#plot (u[:,0],u[:,3],"b--", label = r'$\rho (t)$ \`a 2 mm mod\`ele')
#plot (u[:,0],u[:,4],"g--", label = r'$\rho (t)$ \`a 3 mm mod\`ele')
#plot (u[:,0],u[:,5],"c--", label = r'$\rho (t)$ \`a 4 mm mod\`ele')
#plot (u[:,0],u[:,6],"y--", label = r'$\rho (t)$ \`a 5 mm mod\`ele')

#plot (u[:,0],u[:,1]-z[:,1],"k", label = r'\Delta T$ \`a la paroi')
#plot (u[:,0],u[:,2]-z[:,3],"r", label = r'\Delta T$ \`a 1 mm')
#plot (u[:,0],u[:,3]-z[:,5],"b", label = r'\Delta T$ \`a 2 mm')
#plot (u[:,0],u[:,4]-z[:,7],"g", label = r'\Delta T$ \`a 3 mm')
#plot (u[:,0],u[:,5]-z[:,9],"c", label = r'\Delta T$ \`a 4 mm')
#plot (u[:,0],u[:,6]-z[:,11],"y", label = r'$\Delta T$ \`a 5 mm')

plot (z[:,0],z[:,2],"k", label = r'$\rho$ \`a la paroi')
plot (z[:,0],z[:,4],"r", label = r'$\rho$ \`a 1 mm')
plot (z[:,0],z[:,6],"b", label = r'$\rho$ \`a 2 mm')
plot (z[:,0],z[:,8],"g", label = r'$\rho$ \`a 3 mm')
plot (z[:,0],z[:,10],"c", label = r'$\rho$ \`a 4 mm')
plot (z[:,0],z[:,12],"y", label = r'$\rho$ \`a 5 mm')

plot (u[:,0],u[:,1],"k--", label = r'$\rho$ \`a la paroi mod\`ele')
plot (u[:,0],u[:,2],"r--", label = r'$\rho$ \`a 1 mm mod\`ele')
plot (u[:,0],u[:,3],"b--", label = r'$\rho$ \`a 2 mm mod\`ele')
plot (u[:,0],u[:,4],"g--", label = r'$\rho$ \`a 3 mm mod\`ele')
plot (u[:,0],u[:,5],"c--", label = r'$\rho$ \`a 4 mm mod\`ele')
plot (u[:,0],u[:,6],"y--", label = r'$\rho$ \`a 5 mm mod\`ele')

#title (r"Ecart de la masse volumique (avec pyrolyse, materiau constant)")
# plot (u[:,0],u[:,2],"k", label = r'$\rho (t)$ \`a la paroi')
# plot (u[:,0],u[:,4],"r", label = r'$\rho (t)$ \`a 1 mm')
# plot (u[:,0],u[:,6],"b", label = r'$\rho (t)$ \`a 2 mm')
# plot (u[:,0],u[:,8],"g", label = r'$\rho (t)$ \`a 3 mm')
# plot (u[:,0],u[:,10],"c", label = r'$\rho (t)$ \`a 4 mm')
# plot (u[:,0],u[:,12],"y", label = r'$\rho (t)$ \`a 5 mm')
#
axes = gca()
axes.set_xlim(0, 100)
#



xlabel ('t')
ylabel (r'$\rho$')

legend (loc=0)
show()
