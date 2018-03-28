#coding : utf8

from numpy import *
from matplotlib.pyplot import *
from matplotlib import rc
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})



u = loadtxt ('sol_pts_1/point_0')
v = loadtxt ('sol_pts_1/point_1')
w = loadtxt ('sol_pts_1/point_2')
x = loadtxt ('sol_pts_1/point_3')
y = loadtxt ('sol_pts_1/point_4')

# plot (u[:,0],u[:,1], label = r'$x(t)$ Eul\_exp')
# plot (v[:,0],v[:,1], label = r'$x(t)$ RK3')
# plot (w[:,0],w[:,1], label = r'$x(t)$ RK4')
# plot (x[:,0],x[:,1], label = r'$x(t)$ AB')

plot (u[:,0],u[:,1],"k", label = r'$T(t)$ \`a la paroi')
plot (v[:,0],v[:,1],"r", label = r'$T(t)$ \`a 1 mm')
plot (w[:,0],w[:,1],"b", label = r'$T(t)$ \`a 2 mm')
plot (x[:,0],x[:,1],"g", label = r'$T(t)$ \`a 3 mm')
plot (y[:,0],y[:,1],"c", label = r'$T(t)$ \`a 4 mm')


# plot (u[:,0],u[:,2],"r--", label = r"$\Theta'(t)$ Eul\_exp")
# plot (v[:,0],v[:,2],"b--", label = r"$\Theta'(t)$ RK3")
# plot (w[:,0],w[:,2],"g--", label = r"$\Theta'(t)$ RK4")
# plot (x[:,0],x[:,2],"c--", label = r"$\Theta'(t)$ AB")




# plot (u[:,0],u[:,2], label = r'$x(t)$ Eul\_exp')
# plot (v[:,0],v[:,2], label = r'$x(t)$ RK3')
# plot (w[:,0],w[:,2], label = r'$x(t)$ RK4')
# plot (x[:,0],x[:,2], label = r'$x(t)$ AB')

# plot (u[:,1],u[:,2],"r", label = r'$x(t)$ Eul\_exp')
# plot (v[:,1],v[:,2],"b", label = r'$x(t)$ RK3')
# plot (w[:,1],w[:,2],"g", label = r'$x(t)$ RK4')
# plot (x[:,1],x[:,2],"c", label = r'$x(t)$ AB3')


xlabel ('t')
ylabel (r'$\Theta$')


title (r"\'Evolution de la temp\'erature dans une plaque chauff\'ee en un bord \`a la paroi ainsi  qu\`a 1mm, 2mm, 3mm et 4mm") 

legend (loc=0)
show()

