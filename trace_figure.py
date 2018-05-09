#coding : utf8

from numpy import *
from matplotlib.pyplot import *
from matplotlib import rc
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

z = loadtxt ('ro_thermo_couple_mat_constant')
# z = loadtxt ('ro_thermo_couple_mat_variable')
#
plot (z[:,0],z[:,1],"k--", label = r'$\rho (t)$ \`a la paroi prof')
plot (z[:,0],z[:,2],"r--", label = r'$\rho (t)$ \`a 1 mm prof')
plot (z[:,0],z[:,3],"b--", label = r'$\rho (t)$ \`a 2 mm prof')
plot (z[:,0],z[:,4],"g--", label = r'$\rho (t)$ \`a 3 mm prof')
plot (z[:,0],z[:,5],"c--", label = r'$\rho (t)$ \`a 4 mm prof')
plot (z[:,0],z[:,6],"y--", label = r'$\rho (t)$ \`a 5 mm prof')


# u = loadtxt ('EC_PyrolyseMC_implicit.txt')
v = loadtxt ('EC_PyrolyseMC_explicite.txt')

plot (v[:,0],v[:,2],"k--", label = r'$\rho (t)$ \`a la paroi')
plot (v[:,0],v[:,4],"r--", label = r'$\rho (t)$ \`a 1 mm')
plot (v[:,0],v[:,6],"b--", label = r'$\rho (t)$ \`a 2 mm')
plot (v[:,0],v[:,8],"g--", label = r'$\rho (t)$ \`a 3 mm')
plot (v[:,0],v[:,10],"c--", label = r'$\rho (t)$ \`a 4 mm')
plot (v[:,0],v[:,12],"y--", label = r'$\rho (t)$ \`a 5 mm')

# plot (u[:,0],u[:,2],"k", label = r'$\rho (t)$ \`a la paroi')
# plot (u[:,0],u[:,4],"r", label = r'$\rho (t)$ \`a 1 mm')
# plot (u[:,0],u[:,6],"b", label = r'$\rho (t)$ \`a 2 mm')
# plot (u[:,0],u[:,8],"g", label = r'$\rho (t)$ \`a 3 mm')
# plot (u[:,0],u[:,10],"c", label = r'$\rho (t)$ \`a 4 mm')
# plot (u[:,0],u[:,12],"y", label = r'$\rho (t)$ \`a 5 mm')
#
# axes = gca()
# axes.set_ylim(950, 1550)
#
# title (r"\'Evolution de Rho dans une plaque chauff\'ee en un bord \`a la paroi ainsi  qu \`a 1mm, 2mm, 3mm et 4mm")

z = loadtxt ('temp_thermo_couple_mat_constant')
# z = loadtxt ('temp_thermo_couple_mat_variable')
# z = loadtxt ('temp_thermo_couple_sans_pyrolyse')
#
# plot (z[:,0],z[:,1],"k--", label = r'$T(t)$ \`a la paroi prof')
# plot (z[:,0],z[:,2],"r--", label = r'$T(t)$ \`a 1 mm prof')
# plot (z[:,0],z[:,3],"b--", label = r'$T(t)$ \`a 2 mm prof')
# plot (z[:,0],z[:,4],"g--", label = r'$T(t)$ \`a 3 mm prof')
# plot (z[:,0],z[:,5],"c--", label = r'$T(t)$ \`a 4 mm prof')
# plot (z[:,0],z[:,6],"y--", label = r'$T(t)$ \`a 5 mm prof')
# #
# plot (v[:,0],v[:,1],"k", label = r'$T(t)$ \`a la paroi')
# plot (v[:,0],v[:,3],"r", label = r'$T(t)$ \`a 1 mm')
# plot (v[:,0],v[:,5],"b", label = r'$T(t)$ \`a 2 mm')
# plot (v[:,0],v[:,7],"g", label = r'$T(t)$ \`a 3 mm')
# plot (v[:,0],v[:,9],"c", label = r'$T(t)$ \`a 4 mm')
# plot (v[:,0],v[:,11],"y", label = r'$T(t)$ \`a 5 mm')
#
# plot (u[:,0],u[:,1],"k", label = r'$T(t)$ \`a la paroi')
# plot (u[:,0],u[:,3],"r", label = r'$T(t)$ \`a 1 mm')
# plot (u[:,0],u[:,5],"b", label = r'$T(t)$ \`a 2 mm')
# plot (u[:,0],u[:,7],"g", label = r'$T(t)$ \`a 3 mm')
# plot (u[:,0],u[:,9],"c", label = r'$T(t)$ \`a 4 mm')
# plot (u[:,0],u[:,11],"y", label = r'$T(t)$ \`a 5 mm')
#
# title (r"\'Evolution de la temp\'erature dans une plaque chauff\'ee en un bord \`a la paroi ainsi  qu \`a 1mm, 2mm, 3mm, 4mm et 5mm")

# plot (u[:,1],u[:,2],"k", label = r'$\rho (t)$ \`a la paroi')
# plot (u[:,3],u[:,4],"r", label = r'$\rho (t)$ \`a 1 mm')
# plot (u[:,5],u[:,6],"b", label = r'$\rho (t)$ \`a 2 mm')
# plot (u[:,7],u[:,8],"g", label = r'$\rho (t)$ \`a 3 mm')
# plot (u[:,9],u[:,10],"c", label = r'$\rho (t)$ \`a 4 mm')
# plot (u[:,11],u[:,12],"y", label = r'$\rho (t)$ \`a 5 mm')
#
# title (r"\'Evolution de Rho en fonction de T dans une plaque chauff\'ee en un bord \`a la paroi ainsi  qu \`a 1mm, 2mm, 3mm et 4mm")


xlabel ('T')
ylabel (r'Rho')

legend (loc=0)
show()
