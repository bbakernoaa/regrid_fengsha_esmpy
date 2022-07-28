#!/usr/bin/env python 

from netCDF4 import Dataset
import netCDF4 as nc
import ESMF
import os
import numpy as np
from netCDF4 import date2num,num2date
import datetime as dt
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def create_ds_dict(dset):
    d = {}
    for var in dset.variables.values():
        d[var.name] = dset[var.name][:].copy()
    return d
        
        
def create_esmf_grid_from_gridspec(ds,coords=None):
    if coords != None:
        grid = ESMF.Grid(filename=ds, filetype=ESMF.FileFormat.GRIDSPEC,coord_names=['grid_latt','grid_lont'])
    else:
        grid = ESMF.Grid(filename=ds, filetype=ESMF.FileFormat.GRIDSPEC)
    return grid
    
if __name__ == '__main__':
    
    parser = ArgumentParser(description='Regrid FENGSHA dust inputs to UFS-Regional Domains', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--src', help='input FENGSHA file', type=str, required=True)
    parser.add_argument('-g', '--grid', help='input grid file', required=False)
    parser.add_argument('-o', '--output', help='output file name', default=None, required=False)
    parser.add_argument('-r', '--regrid_method', help='output file name', default='nearest', required=False)
    parser.add_argument('-w', '--weight_file', help='ESMF Weight File; if not generated will output to this', default='ESMF_Weight_file.nc', required=False)
    args = parser.parse_args()
    
    # regrid method 
    method = args.regrid_method
    if method.lower() == 'nearest':
        RegridMethod = ESMF.RegridMethod.NEAREST_STOD
    elif method.lower() == 'bilinear':
        RegridMethod = ESMF.RegridMethod.BILINEAR
    elif method.lower() == 'patch':
        RegridMethod = ESMF.RegridMethod.PATCH
        
    regrid_file_found=False
    if os.path.isfile(args.weight_file):
        regrid_file_found=True
        print('Weight File Found.... Will Re-use.')
        print('If you want to regenerate remove or choose a different filename')
    else:
        print('Weight File NOT Found.... Generating')
    # open the nexus output
    ds = nc.Dataset(args.src)
    ds_dict = create_ds_dict(ds)
    
    
    # open the grid file
    grid_ds = nc.Dataset(args.grid)
    grid_dict = create_ds_dict(grid_ds)

    # create ESMF grids
    tgt_shape = grid_dict['grid_latt'].shape
    tgt_latt = grid_dict['grid_latt']
    tgt_lont = grid_dict['grid_lont']
    tgt_lat  = grid_dict['grid_lat']
    tgt_lon  = grid_dict['grid_lon']

    target = ESMF.Grid(np.array(tgt_shape), staggerloc=[ESMF.StaggerLoc.CENTER, ESMF.StaggerLoc.CORNER],coord_sys=ESMF.CoordSys.SPH_DEG)
    tgt_cen_lon = target.get_coords(0, staggerloc=ESMF.StaggerLoc.CENTER)
    tgt_cen_lat = target.get_coords(1, staggerloc=ESMF.StaggerLoc.CENTER)
    tgt_con_lon = target.get_coords(0, staggerloc=ESMF.StaggerLoc.CORNER)
    tgt_con_lat = target.get_coords(1, staggerloc=ESMF.StaggerLoc.CORNER)

    tgt_cen_lon[...] = tgt_lont
    tgt_cen_lat[...] = tgt_latt
    
    tgt_con_lon[...] = tgt_lon
    tgt_con_lat[...] = tgt_lat

    
    # the source is a COARDS format --- can just read 
    source = create_esmf_grid_from_gridspec(args.src)
    
    # create regridding obj 
    x,y = source.upper_bounds[0]
    #srcfield.data[...] = np.zeros((x,y))
    
    source_field = ESMF.Field(source, staggerloc=ESMF.StaggerLoc.CENTER)
    target_field = ESMF.Field(target, staggerloc=ESMF.StaggerLoc.CENTER)
    
    if regrid_file_found == False:
        regrid = ESMF.Regrid(source_field,target_field, regrid_method=RegridMethod, filename=args.weight_file) 
    
    regrid = ESMF.RegridFromFile(source_field,target_field, args.weight_file)

    # regrid fields
    regrid_field=False
    regridded_dict = {} 
    for v in ['uthres','clayfrac','ssm','sandfrac','albedo_drag']:
        if (v is 'ssm') | (v =='albedo_drag'):
            # this has multiple time slices fields
            if v in ds_dict.keys():
                timeslices = 12
                regrid_field=True
            else:
                regrid_field = False
        else:
            if v in ds_dict.keys():
                timeslices = 1
                regrid_field=True
            else:
                regrid_field = False
        if regrid_field:
            if timeslices == 1:
                srcField = ESMF.Field(source, staggerloc=ESMF.StaggerLoc.CENTER)
                srcField.data[...] = ds_dict[v][:,:].T
                regridded_dict[v] = regrid(srcField,target_field).data.copy()
            else:
                x,y = target_field.data.shape
                tmp_array = np.zeros((12,x,y))
                srcField = ESMF.Field(source, staggerloc=ESMF.StaggerLoc.CENTER)
                for i in range(12):
                    srcField.data[...] = ds_dict[v][i,:,:].T
                    tmp_array[i,:,:] = regrid(srcField,target_field).data.copy()
                    
                regridded_dict[v] = tmp_array

    # open file
    ncfile = nc.Dataset(args.output,mode='w',format='NETCDF4')
    
    #create dims
    x_dim = ncfile.createDimension('x',grid_ds['grid_xt'].shape[0])
    y_dim = ncfile.createDimension('y',grid_ds['grid_yt'].shape[0])
    time_dim = ncfile.createDimension('time',None)
    
    # add attributes
    ncfile.title='FENGSHA Dust Emission Input Data'
    ncfile.contact='barry.baker@noaa.gov'
    ncfile.created=dt.date.today().strftime('Created %Y-%m-%d')
    
    # create variables - lat and lon and time
    lat = ncfile.createVariable('latitude',np.float32,('y','x',))
    lat.long_name = 'latitude'
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('longitude',np.float32,('y','x',))
    lon.long_name='longitude'
    lon.units = 'degrees_east'
    
    time = ncfile.createVariable('time',np.float64,('time',))
    time.long_name = 'time'
    time.units = ds['time'].units
    time[:] = ds['time'][:]

    # create all other variables
    nc_var_dict = {}
    for var in regridded_dict.keys():
        if (var is 'ssm') | (var =='albedo_drag'):
            nc_var_dict[var] = ncfile.createVariable(var, np.float32, ('time','y','x',), zlib=True)
        else:
            nc_var_dict[var] = ncfile.createVariable(var, np.float32, ('y','x',), zlib=True)
        if var == 'uthres':
            nc_var_dict[var].units = 'm s-1'
        else:
            nc_var_dict[var].units='-'
            
    #populate variables

    # # first with lat and lon and time
    lat[:] = grid_dict['grid_latt'][:]
    lon[:] = grid_dict['grid_lont'][:]

    # # now the other variables
    for var in nc_var_dict.keys():
        if var == 'uthres':
            nc_var_dict[var][:] = np.nan_to_num(regridded_dict[var][:], posinf=1.0e8, neginf=1.0e8, nan=1.0e8 )
        else:
            nc_var_dict[var][:] = np.nan_to_num(regridded_dict[var][:], posinf=-1.0, neginf=-1.0, nan=-1.0 )    
    
    ncfile.close()
