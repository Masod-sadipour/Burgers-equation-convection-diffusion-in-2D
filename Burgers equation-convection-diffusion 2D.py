## Solving burgers equation in 2D
## May 2020.

# importing necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# variables and discritization parameters

nt=500
nx=41
ny=41
c=1

nu = .1
dt = .001

dx=2/(nx-1)
dy=2/(ny-1)

x = np.linspace(0, 2, nx)
y = np.linspace(0, 2, ny)

comb = np.ones((ny, nx))
u = np.ones((ny, nx))
v = np.ones((ny, nx))
un = np.ones((ny, nx))
vn = np.ones((ny, nx))
uf=np.ones((nt,nx,ny))
vf=np.ones((nt,nx,ny))
# assigning initial conditions
u[int(0.75 / dy):int(1.25 / dy + 1),int(0.75 / dx):int(1.25 / dx + 1)] = 3
v[int(0.75 / dy):int(1.25 / dy + 1),int(0.75 / dx):int(1.25 / dx + 1)] = 3

###(ploting ICs)
fig = plt.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x, y)
ax.plot_surface(X, Y, u[:], cmap=cm.jet)

plt.title('U')
plt.xlabel('X')
plt.ylabel('Y');
fig.savefig('U-I.C.png', bbox_inches='tight')

##loop across number of time steps
for n in range(nt):
    un = u.copy()
    vn = v.copy()
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            u[i,j] = (un[i, j] -(un[i, j] * dt / dx * (un[i, j] - un[i-1, j])) -vn[i, j] * dt / dy * (un[i, j] - un[i, j-1])) + (nu*dt/(dx**2))*(un[i+1,j]-2*un[i,j]+un[i-1,j])+(nu*dt/(dx**2))*(un[i,j-1]-2*un[i,j]+un[i,j+1])
            v[i,j] = (vn[i, j] -(un[i, j] * dt / dx * (vn[i, j] - vn[i-1, j]))-vn[i, j] * dt / dy * (vn[i, j] - vn[i, j-1])) + (nu*dt/(dx**2))*(vn[i+1,j]-2*vn[i,j]+vn[i-1,j])+(nu*dt/(dx**2))*(vn[i,j-1]-2*vn[i,j]+vn[i,j+1])
            uf[n,i,j]=u[i,j]  # U in every time-step
            vf[n, i, j] = v[i, j]  # V in every time-step
    # Velocity boundary conditions
    u[:,0]=1
    u[:,-1]=1
    u[0,:]=1
    u[-1,:]=1
    v[0, :] = 1
    v[-1, :] = 1
    v[:, 0] = 1
    v[:, -1] = 1

# plotting U field as a Surface
fig = plt.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x, y)
ax.plot_surface(X, Y, uf[499], cmap=cm.jet, rstride=1, cstride=1)
plt.title('U')
plt.xlabel('X')
plt.ylabel('Y');
fig.savefig('U.png', bbox_inches='tight')