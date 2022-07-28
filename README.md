# regrid_fengsha_esmpy

This is a simple utility to use ESMPy to regrid the FENGSHA inputs for the UFS Weather Model Regional Grids.  
The purpose is to be able to support regridding on NOAA's production machines with available python environments.

```bash
./esmpy_interp_fengsha.py -h
usage: esmpy_interp_fengsha.py [-h] -s SRC [-g GRID] [-o OUTPUT]
                               [-r REGRID_METHOD] [-w WEIGHT_FILE]

Regrid FENGSHA dust inputs to UFS-Regional Domains

optional arguments:
  -h, --help            show this help message and exit
  -s SRC, --src SRC     input FENGSHA file (default: None)
  -g GRID, --grid GRID  input grid file (default: None)
  -o OUTPUT, --output OUTPUT
                        output file name (default: None)
  -r REGRID_METHOD, --regrid_method REGRID_METHOD
                        output file name (default: nearest)
  -w WEIGHT_FILE, --weight_file WEIGHT_FILE
                        ESMF Weight File; if not generated will output to this
                        (default: ESMF_Weight_file.nc)
```


## Example Use


```bash
./esmpy_interp_fengsha.py -s /scratch1/RDARCH/rda-arl-gpu/Barry.Baker/emissions/NASA/ExtData/Dust/FENGSHA_p81_10km_inputs.nc -g grid_spec_RRFS_CONUS_13km.nc -o testc.nc
```
## Requirements 

- netCDF4 
- ESMPy
- os
- numpy 
- datetime
