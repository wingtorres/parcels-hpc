import sys
import xarray as xr
import xgcm
import numpy as np
from datetime import timedelta as delta
from parcels import AdvectionRK4, ErrorCode, FieldSet, JITParticle, ParticleFile, ParticleSet

filenames =  "/work/wtorres/particles/*nc"
ds = xr.open_mfdataset(filenames, chunks={'ocean_time': 12}, combine="by_coords", parallel=True, decode_times = False)

#non-redundant dims: xi_rho (center), eta_rho (center), xi_u (inner), eta_v (inner)
ds = ds.rename({'eta_u': 'eta_rho', 'xi_v': 'xi_rho', 'xi_psi': 'xi_u', 'eta_psi': 'eta_v'})

coords={'xi':{'center':'xi_rho', 'inner':'xi_u'}, 
        'eta':{'center':'eta_rho', 'inner':'eta_v'}, 
        's':{'center':'s_rho', 'outer':'s_w'}}

grid = xgcm.Grid(ds, coords=coords)

#Calculate Lagrangian velocity
ds['u_lagrangian'] = ds['u'] + ds['u_stokes']
ds['v_lagrangian'] = ds['v'] + ds['v_stokes']

#If applicable, orient Lagrangian velocity eastward/northward. ds.angle is the grid angle.
# ds['angle_psi'] = grid.interp(grid.interp(ds.angle,'eta'), 'psi')
# ds['uveitheta'] = (ds.u_lagrangian_psi + 1j*ds.v_lagrangian_psi)*np.exp(1j*ds.angle_psi) 
# ds['u_lagrangian '] = np.real(ds.uveitheta)
# ds['v_lagrangian'] = np.imag(ds.uveitheta)

variables = {'U': 'u_lagrangian',
             'V': 'v_lagrangian'}

dimensions = {'U': {'lon': 'lon_u', 'lat': 'lat_u', 'depth': 's_rho', 'time': 'ocean_time'},
              'V': {'lon': 'lon_v', 'lat': 'lat_v', 'depth': 's_rho', 'time': 'ocean_time'}}

fieldset = FieldSet.from_xarray_dataset(ds, variables = variables, dimensions = dimensions, mesh = 'spherical', allow_time_extrapolation = True)

def DeleteParticle(particle, fieldset, time):
    print("Deleting particle")
    particle.delete()

recovery = {ErrorCode.ErrorOutOfBounds: DeleteParticle,
            ErrorCode.ErrorThroughSurface: DeleteParticle}

#cover domain in particles
lon = np.linspace(ds.lon_rho.min(), ds.lon_rho.max(), num=2**6)
lat = np.linspace(ds.lat_rho.min(), ds.lat_rho.max(), num=2**6) 
lons, lats = np.meshgrid(lon,lat)

pset = ParticleSet.from_list(fieldset = fieldset, pclass = JITParticle, time = ds.ocean_time.values[0], lon = lons, lat = lats, depth = depth )

kernels = AdvectionRK4 
output_file = pset.ParticleFile(name= "/work/wtorres/temp/reefParticles", outputdt = delta(seconds = 60) )
pset.execute( kernels, runtime = delta(hours = 12.0), dt = delta(seconds = 60), output_file = output_file, recovery = recovery )
output_file.export()
output_file.close()

