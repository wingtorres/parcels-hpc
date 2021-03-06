## Miniconda
First a python distribution must be installed. I found it best to use Miniconda, a lightweight version of Anaconda's python distribution with just conda and its dependencies, which was straightforward to set up following the excellent [guide](https://medium.com/@rabernat/custom-conda-environments-for-data-science-on-hpc-clusters-32d58c63aa95) by @rabernat. Instead of using his example environment.yml, I've provided a parcels.yml file in this repository that will create an environment for Parcels along with some useful dependencies.

## The Parcels python script
This script must convert output files from ROMS into a Parcels [Fieldset](http://oceanparcels.org/gh-pages/html/#module-parcels.fieldset), initialize the location and timing of particles to be released through Parcels' [Particleset](http://oceanparcels.org/gh-pages/html/#module-parcels.particleset) module, and of course execute the integration itself. I use the packages [xarray](http://xarray.pydata.org/en/stable/) and [xgcm](https://xgcm.readthedocs.io/en/latest/) to facilitate reading netCDF data and interpolation on a structured grid.

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
```

where chunk size can be adjusted to improve performance. From the xarray [docs](http://xarray.pydata.org/en/stable/dask.html#chunking-and-performance)...
> A good rule of thumb is to create arrays with a minimum chunksize of at least one million elements (e.g., a 1000x1000 matrix).

3) Create an xgcm grid object to facilitate interpolation

```
#non-redundant dims: xi_rho (center), eta_rho (center), xi_u (inner), eta_v (inner)
ds = ds.rename({'eta_u': 'eta_rho', 'xi_v': 'xi_rho', 'xi_psi': 'xi_u', 'eta_psi': 'eta_v'})

coords={'xi':{'center':'xi_rho', 'inner':'xi_u'}, 
        'eta':{'center':'eta_rho', 'inner':'eta_v'}, 
        's':{'center':'s_rho', 'outer':'s_w'}}

grid = xgcm.Grid(ds, coords=coords)
```

4) Calculate Lagrangian velocities. The depth-averaged velocity is used here to keep the example in 2D.

```
#Calculate Lagrangian velocity
ds['ubar_lagrangian'] = ds['ubar'] + ds['ubar_stokes']
ds['vbar_lagrangian'] = ds['vbar'] + ds['vbar_stokes']
```

4.5) If applicable, orient Lagrangian velocities eastward/northward, requiring an interpolation to a common coordinate if necessary. On the staggered Arakawa-C grid used in ROMS, it is best to interpolate the velocities to ψ-points - see https://www.myroms.org/wiki/Numerical_Solution_Technique for more details on the discretization. Here, ds.angle is the grid angle.

```
ds['ubar_lagrangian_psi'] = grid.interp(ds.ubar_lagrangian, 'eta')
ds['vbar_lagrangian_psi'] = grid.interp(ds.vbar_lagrangian, 'xi')
ds['angle_psi'] = grid.interp(grid.interp(ds.angle,'eta'), 'psi')
ds['uveitheta'] = (ds.ubar_lagrangian_psi + 1j*ds.vbar_lagrangian_psi)*np.exp(1j*ds.angle_psi) 
ds['ubar_lagrangian_psi'] = np.real(ds.uveitheta)
ds['vbar_lagrangian_psi'] = np.imag(ds.uveitheta)
```


5) Define a Parcels [Fieldset](http://oceanparcels.org/gh-pages/html/#module-parcels.fieldset) from the xarray dataset object. The velocities can be on different grid like for an Arakawa-C grid, just make sure the dimensions dictionary is consistent with how the variables are referenced.

```
variables = {'U': 'ubar_lagrangian',
             'V': 'vbar_lagrangian'}

dimensions = {'U': {'lon': 'lon_u', 'lat': 'lat_u', 'depth': 's_rho', 'time': 'ocean_time'},
              'V': {'lon': 'lon_v', 'lat': 'lat_v', 'depth': 's_rho', 'time': 'ocean_time'}}
fieldset = FieldSet.from_xarray_dataset(ds, variables = variables, dimensions = dimensions, mesh = 'spherical')
```

6) Define a recovery kernel to handle particles going out of bounds

```
def DeleteParticle(particle, fieldset, time):
    print("Deleting particle")
    particle.delete()

recovery = {ErrorCode.ErrorOutOfBounds: DeleteParticle,
            ErrorCode.ErrorThroughSurface: DeleteParticle}
```

7) Initialize a Parcels [Particleset](http://oceanparcels.org/gh-pages/html/#module-parcels.particleset). Here I seed the whole domain with equally spaced particles

```
lon = np.linspace(ds.lon_rho.min(), ds.lon_rho.max(), num=2**5)
lat = np.linspace(ds.lat_rho.min(), ds.lat_rho.max(), num=2**5) 
lons, lats = np.meshgrid(lon,lat)

pset = ParticleSet.from_list(fieldset = fieldset, pclass = JITParticle, time = ds.ocean_time.values[0], lon = lons, lat = lats, depth = depth )
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
I have provided a sample .pbs script in this repository that activates the Parcels python environment and runs the python script described above on the HPC resource with the command
```
qsub parcels.pbs
``` 
Make sure the environment name, allocation, paths to the error and output files, email, and script location match your configuration.
