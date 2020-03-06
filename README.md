## Miniconda
First a python distribution must be installed. I found it best to use Miniconda, a lightweight version of Anaconda's python distribution with just conda and its dependncies, which was straightforward to set up following the excellent [guide](https://medium.com/@rabernat/custom-conda-environments-for-data-science-on-hpc-clusters-32d58c63aa95) by @rabernat. Instead of using his example environment.yml, I've provided a parcels.yml file in this repository that will create an environment for Parcels along with some useful dependencies.

## The Parcels python script
This script must convert output files from ROMS into a Parcels Fieldset, initialize the location and timing of particles to be released, and of course execute the integration itself. I use the packages [xarray](http://xarray.pydata.org/en/stable/) and [xgcm](https://xgcm.readthedocs.io/en/latest/) to facilitate reading netCDF data and interpolation on a structured grid.

1) Import packages

```
import sys 
import xarray as xr
import xgcm
import numpy as np
from datetime import timedelta as delta
from parcels import AdvectionRK4, ErrorCode, FieldSet, JITParticle, ParticleFile, ParticleSet
```

2) Load in multiple netcdf files as an xarray dataset

```
filenames =  "/work/wtorres/particles/*nc"
ds = xr.open_mfdataset(filenames, chunks={'ocean_time': 12}, combine="by_coords", parallel=True, decode_times = False)

x = ds['x_psi'].values
y = ds['y_psi'].values
t = ds['ocean_time'].values
```

where chunk size can be adjusted to improve performance. From the xarray [docs](http://xarray.pydata.org/en/stable/dask.html#chunking-and-performance)...
> A good rule of thumb is to create arrays with a minimum chunksize of at least one million elements (e.g., a 1000x1000 matrix).

3) Create an xgcm grid object to facilitate interpolation

```
#non-redundant dims: xi_rho (outer), eta_rho (outer), xi_psi (inner), eta_psi (inner)
ds = ds.rename({'eta_u': 'eta_rho', 'xi_v': 'xi_rho', 'xi_psi': 'xi_u', 'eta_psi': 'eta_v'})

coords={'xi':{'center':'xi_rho', 'inner':'xi_u'}, 
        'eta':{'center':'eta_rho', 'inner':'eta_v'}, 
        's':{'center':'s_rho', 'outer':'s_w'}}

grid = xgcm.Grid(ds, coords=coords)
```

4) Interpolate the horizontal velocity components to ψ points (check out https://www.myroms.org/wiki/Numerical_Solution_Technique for more details on the discretization), and compute the Lagrangian velocity as a sum of the Eulerian and Stokes velocities. The depth-averaged velocity is used here to keep the example in 2D.

```
#velocity to psi points
ds['ubar_lagrangian'] = ds['ubar'] + ds['ubar_stokes']
ds['vbar_lagrangian'] = ds['vbar'] + ds['vbar_stokes']

ds['ubar_lagrangian_psi'] = grid.interp(ds.u_lagrangian, 'eta')
ds['vbar_lagrangian_psi'] = grid.interp(ds.v_lagrangian, 'xi')
```

5) Define a Parcels [Fieldset](http://oceanparcels.org/gh-pages/html/#module-parcels.fieldset) from the lagrangian velocity field and horizontal coordinates of ψ points. Even if the grid is in meters, the dimensions must be named 'lon' and 'lat' for consistency.

```
data = {'U': da.ubar_lagrangian_psi.values, 'V': ds.vbar_lagrangian_psi.values}
dimensions = {'lon': x, 'lat': y, 'time': t }
fieldset = FieldSet.from_data(data, dimensions, transpose = False, mesh = 'flat') #mesh = 'flat' for cartesian, 'spherical' for curvilinear 
```

6) Define a recovery kernel to handle particles going out of bounds

```
def DeleteParticle(particle, fieldset, time):
    print("Deleting particle")
    particle.delete()

recovery = {ErrorCode.ErrorOutOfBounds: DeleteParticle,
            ErrorCode.ErrorThroughSurface: DeleteParticle}
```

7) Initialize a Parcels [Particleset](http://oceanparcels.org/gh-pages/html/#module-parcels.particleset). Here I seed 1000 particles along the line y = 1000m, although there are several other ways to do this.
```
pset = ParticleSet.from_line(fieldset = fieldset, size = 1000, start = (1000, 1000), finish = (8000,1000), pclass = JITParticle )
```

8) Specify advection kernel (4th order Runge-Kutta is selected here), choose name and output time step for particle locations, and execute Parcels

```
kernels = AdvectionRK4
output_file = pset.ParticleFile(name= "/work/wtorres/temp/reefParticles", outputdt = delta(seconds = 60) )
pset.execute( kernels, runtime = delta(hours = 12.0), dt = delta(seconds = 60), output_file = output_file, recovery = recovery )
output_file.export()
output_file.close()
```

## Submitting a Parcels job
I have provided a sample .pbs script in this repository that runs the python script described above on the HPC resource with the command
```
qsub parcels.pbs
``` 
Make sure the allocation name, paths to the error and output files, email, and script location match your configuration.
