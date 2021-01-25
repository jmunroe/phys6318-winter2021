import numpy as np
import xarray as xr
import netCDF4
import time
import os
import glob


## global, fixed parameters

# grid parameters
xmax = 1010
dx = 10 # grid spacing
nx = 101
nz = 10
dt = 0.1

# physical parameters
g = 9.81 # acceleration due to gravity
hmin = 0.01 # minimum layer thickness

eps = 0.05


x_eta = np.empty(nx+2) # horizontal grid for depths and displacements
x_u = np.empty(nx+2) # horizontal grid for horizontal velocities
z = np.empty(nz+2) # vertical grid

htotal = np.empty(nx+2) # initial bathymetry
hzero = np.empty((nz+2, nx+2)) # initial layer thicknesses
h = np.empty((nz+2, nx+2)) # actual layer thicnesses
dp = np.empty((nz+2, nx+2)) # dynamic pressure

eta = np.empty((nz+2, nx+2)) # actual interface displacements
etan = np.empty((nz+2, nx+2)) # interface displacements at time level n+1
eta0 = np.empty((nz+2, nx+2)) # initial interface displacements
dhdt = np.empty((nz+2, nx+2)) # actual layer-thickness change
u = np.empty((nz+2, nx+2)) # actual lateral velocity
un = np.empty((nz+2, nx+2)) # lateral velocity at time level n+1
rho = np.empty(nz+2) # layer densities

wet = np.empty((nz+2, nx+2), dtype=int)


# In[383]:


def init_output(basename='', **attrs):
      
    # ensure we are writing to a new file each time
    previous_files = sorted(glob.glob(f'OUTPUT/{basename}_*.nc'))
    if len(previous_files) == 0:
        counter = 1    
    else:
        counter = int(previous_files[-1][-7:-3]) + 1
    filename = f'OUTPUT/{basename}_{counter:04d}.nc'
        
    nc = netCDF4.Dataset(filename, "w")

    # NetCDF files have 'dimensions'
    nc.createDimension("time")
    nc.createDimension("x", nx+2)
    nc.createDimension("z", nz+2)

    # NetCDF files have 'variables'
    nc.createVariable("time", "f8", ("time",))
    # NetCDF variables can have 'attributes'
    nc.variables["time"].units = "seconds since 2000-01-01"
    nc.variables["time"].calendar = "gregorian"
    
    nc.createVariable("x_eta", "f8", ("x"))
    nc.variables["x_eta"].units = "m"
    nc.variables["x_eta"][:] = x_eta
    
    nc.createVariable("x_u", "f8", ("x"))
    nc.variables["x_u"].units = "m"
    nc.variables["x_u"][:] = x_u
      
    nc.createVariable("htotal", "f8", ("x"))
    nc.variables["htotal"].units = "m"
    nc.variables["htotal"][:] = htotal
    
    nc.createVariable("eta", "f8", ("time", "z", "x"))
    nc.variables["eta"].units = "m"
    nc.createVariable("u", "f8", ("time", "z", "x"))
    nc.variables["u"].units = "m s-1"
    nc.createVariable("h", "f8", ("time", "z", "x"))
    nc.variables["h"].units = "m"

    # NetCDF files have also have global attributes
    nc.setncatts(attrs)
    nc.history = "Created " + time.ctime(time.time())
    nc.source = "OMB Exercise 7"

    # It is important to close a NetCDF file
    nc.close()
    
    return filename

def write_output(filename, t, n=None):
    # open a NetCDF file in 'append' mode
    nc = netCDF4.Dataset(filename, mode='a')
    
    # if n is not provided, place values in the last position
    if n is None:
        n = len(nc.variables['time'])
        
    nc.variables["time"][n] = t
    nc.variables["eta"][n, :, :] = eta
    nc.variables["h"][n, :, :] = h
    nc.variables["u"][n, :, :] = u
    
    # It is important to close a NetCDF file
    nc.close()


# ### Subroutines

# In[384]:


def init():

    # calculate horizontal grids
    x_eta[:] = np.arange(-0.75 * dx, xmax+dx, dx)
    x_u[:] = x_eta + 0.5 * dx
    
    # bathymetry
    for k in range(1, nx+1):
        htotal[k] = 100
        
    # triangle-shaped island
    for k in range(31, 52):
        htotal[k] = 100 - 95*(k-30)/21
    for k in range(52, 72):
        htotal[k] = 100 - 95*(71 - k + 1)/20

    htotal[0] = -10
    htotal[nx+1] = -10
    
    # undisturbed layer thicknesses & interface displacements
    hini = np.ones(nz+2)*10
        
    for k in range(0, nx+2):
        htot = htotal[k]
        for i in range(1, nz+1):
            hzero[i, k] = max( min( hini[i], htot), 0)
            eta[i, k] = max(0, -htot)
            htot = htot - hini[1]
            
    # layer densities
    rho[0] = 0 # air density ignored
    rho[1] = 1025
    for i in range(2, nz+1):
        rho[i] = 1026 + (i-2)/(nz-2)*0.5
        
    # boundary values for dp and eta
    for k in range(0, nx+2):
        dp[0, k] = 0 # air pressure ignored
        eta[nz+1, k] = 0 # sea floor is rigid
            
    # store initial interface displacements
    for k in range(0, nx+2):
        for i in range(1, nz+2):
            eta0[i, k] = eta[i, k]
            
    # layer thicknesses, wet\dry pointers and velocities
    for i in range(1, nz+1):
        for k in range(0, nx+2):
            h[i, k] = hzero[i, k]
            wet[i, k] = 1
            if h[i, k] < hmin:
                wet[i, k] = 0
            u[i, k] = 0
            un[i, k] = 0

        
def dyn():
    
    # calculate dynamic pressure
    for k in range(0, nx+2):
        for i in range(1, nz+1):
            dp[i, k] = dp[i-1, k] + (rho[i] - rho[i-1])*g*eta[i,k]
            
    for k in range(1, nx+1):
        for i in range(1, nz+1):
            
            # velocity predictor for wet grid cells
            pgradx = -(dp[i, k+1]- dp[i, k])/rho[i]/dx
            un[i, k] = 0
            
            if wet[i, k]:
                if wet[i, k+1] or (pgradx>0):
                    un[i, k] = u[i, k] + dt*pgradx
            else:
                if wet[i, k+1] and (pgradx<0):
                    un[i, k] = u[i, k] + dt*pgradx

    # layer-thickness change predictor
    for k in range(1, nx+1):
        for i in range(1, nz+1):
            hep = 0.5*(un[i,k]+abs(un[i,k]))*h[i,k]
            hen = 0.5*(un[i,k]-abs(un[i,k]))*h[i,k+1]
            hue = hep+hen
            hwp = 0.5*(un[i,k-1]+abs(un[i, k-1]))*h[i,k-1]
            hwn = 0.5*(un[i,k-1]-abs(un[i, k-1]))*h[i,k]
            huw = hwp+hwn

            dhdt[i,k] = -(hue-huw)/dx

            
    # update interface displacements
    for k in range(1, nx+1):
        deta = 0
        for i in range(nz, 0, -1):
            deta = deta + dhdt[i, k]
            etan[i, k] = eta[i, k] + dt*deta
            
    # apply Shapiro filter
    shapiro()
    
    # update layer thicknesses, lateral velocities and wet/dry pointers
    for k in range(1, nx+1):
        for i in range(1, nz+1):
            h[i, k] = hzero[i, k] + eta[i, k] - eta[i+1, k] - eta0[i, k] + eta0[i+1, k]
            u[i, k] = un[i, k]
            wet[i, k] = 1
            if h[i, k] < hmin:
                wet[i, k] = 0
    
def shapiro():
    for i in range(1, nz+1):
        for k in range(1, nx+1):
            if wet[i, k]:
                term1 = (1.0-0.5*eps*(wet[i, k+1]+wet[i, k-1]))*etan[i, k]
                term2 = 0.5*eps*(wet[i, k+1]*etan[i, k+1]+wet[i, k-1]*etan[i, k-1])
                eta[i, k] = term1 + term2
            else:
                eta[i, k] = etan[i, k]


# In[385]:


def multi():
    
    # initialize arrays to the initial values
    init()
    
    # determine maximum water depth
    # automatic setting of time step (10% below CFL threshold)
    hmax = 100 # total water depth
    dt = 0.9 * dx / np.sqrt(g*hmax)
    print(f"Time step = {dt:.1f} seconds")
    
    # set epsilon for Shapiro filter
    eps = 0.05
        
    # parameters for wave paddle
    Apaddle = 1 # amplitude in metres
    Tpaddle = 10 # period in seconds CASE 1
    #Tpaddle = 2*3600 # period in seconds CASE 2
    
    
    # runtime parameters
    tmax = 10 * Tpaddle
    ntot = int(tmax/dt)

    # output parameter
    tout = Tpaddle/10
    nout = int(tout / dt)
    
    filename = init_output('multi')
     
    for n in range(ntot):
        t = n*dt
        
        for i in range(1, nz+1):
            eta[i, 1] = Apaddle*np.sin(2*np.pi*t/Tpaddle)

        #---- prognostic equations ----
        dyn()
        #------------------------------

        if n % nout == 0:
            write_output(filename, t)
            
    return filename


filename = multi()
print(filename)


