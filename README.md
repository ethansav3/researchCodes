# Simba lightcone generation

## Requirements
- numpy
- astropy
- caesar

## Running

`python generate_lightcone.py area zmin zmax`
```
$ python generate_lightcone.py -h
usage: generate_lightcone.py [-h] [--N N] [-f FILENAME] [--overwrite] [--halos] [-v] area z_min z_max

Generate a Simba lightcone.

positional arguments:
  area                  Survey area (deg^2)
  z_min                 Minimum redshift
  z_max                 Maximum redshift

optional arguments:
  -h, --help            show this help message and exit
  --N N                 Number of lightcones. (default: 1)
  -f FILENAME, --filename FILENAME
                        Output filename (Full path from current directory. Must include string
                        formatting to iterate over lightcone number) (default:
                        output/halo_lightcone_%03d.h5)
  --overwrite           Overwrite existing datasets (default: False)
  --halos               Output halos, rather than galaxies. (default: False)
  -v, --verbose         Increase output verbosity (default: False)

```

`N` is the number of lightcones you wish to generate, `A` is the sky area in deg<sup>2</sup> (must be small enough that it covers the whole box, as lateral tiling is not yet implemented; 0.5 deg<sup>2</sup> works), and `zmin` and `zmax` are the lower and upper redshift limits. You can generate lightcones for both galaxy and halo catalogues (warning: a halo lightcone takes longer).

## Output

Each lightcone is saved in its own HDF5 file in the `output` folder, specified by `filename`. The HDF5 is structured so each top group is a snapshot key, within which are datasets for RA, DEC, index and redshift ('z'). The index refers back to the caesar galaxy (halo) arrays for that snapshot.
