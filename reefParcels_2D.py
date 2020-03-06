import sys
import xarray as xr
import xgcm
import numpy as np
from datetime import timedelta as delta
from parcels import AdvectionRK4, ErrorCode, FieldSet, JITParticle, ParticleFile, ParticleSet

filenames =  "/work/wtorres/particles/*nc"
ds = xr.open_mfdataset(filenames, chunks={'ocean_time': 12}, combine="by_coords", parallel=True, decode_times = False)

x = Ds['x_psi'].values
y = Ds['y_psi'].values
t = Ds['ocean_time'].values

#non-redundant dims: xi_rho (center), eta_rho (center), xi_u (inner), eta_v (inner)
ds = ds.rename({'eta_u': 'eta_rho', 'xi_v': 'xi_rho', 'xi_psi': 'xi_u', 'eta_psi': 'eta_v'})

coords={'xi':{'center':'xi_rho', 'inner':'xi_u'}, 
        'eta':{'center':'eta_rho', 'inner':'eta_v'}, 
        's':{'center':'s_rho', 'outer':'s_w'}}

grid = xgcm.Grid(ds, coords=coords)

#velocity to psi points
ds['ubar_lagrangian'] = ds['ubar'] + ds['ubar_stokes']
ds['vbar_lagrangian'] = ds['vbar'] + ds['vbar_stokes']

ds['ubar_lagrangian_psi'] = grid.interp(ds.u_lagrangian, 'eta')
ds['vbar_lagrangian_psi'] = grid.interp(ds.v_lagrangian, 'xi')

dimensions = {'lon': x, 'lat': y, 'time': t }
fieldset = FieldSet.from_data(data, dimensions, transpose = False, mesh = 'flat') 

def DeleteParticle(particle, fieldset, time):
    print("Deleting particle")
    particle.delete()

recovery = {ErrorCode.ErrorOutOfBounds: DeleteParticle,
            ErrorCode.ErrorThroughSurface: DeleteParticle}

pset = ParticleSet.from_line(fieldset = fieldset, size = 1000, start = (1000, 1000), finish = (8000,1000), pclass = JITParticle )

kernels = AdvectionRK4 
output_file = pset.ParticleFile(name= "/work/wtorres/temp/reefParticles", outputdt = delta(seconds = 60) )
pset.execute( kernels, runtime = delta(hours = 12.0), dt = delta(seconds = 60), output_file = output_file, recovery = recovery )
output_file.export()
output_file.close()

