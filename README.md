# regrid_fengsha_esmpy

This is a simple utility to use ESMPy to regrid the FENGSHA inputs for the UFS Weather Model Regional Grids

## Example Use


```bash
./esmpy_interp_fengsha.py -s /scratch1/RDARCH/rda-arl-gpu/Barry.Baker/emissions/NASA/ExtData/Dust/FENGSHA_p81_10km_inputs.nc -g grid_spec_RRFS_CONUS_13km.nc -o testc.nc
```
## Requirements 

- netCDF4 
- ESMPy
- os
- numpy 