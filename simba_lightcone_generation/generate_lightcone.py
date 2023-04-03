import os
import sys
import argparse

import h5py
import numpy as np
from astropy.cosmology import z_at_value
import astropy.units as u

from simba import simba


sb = simba()
verbose = True

parser = argparse.ArgumentParser(description='Generate a Simba lightcone.', add_help=True,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('area', type=float, help="Survey area (deg^2)")
parser.add_argument('z_min', type=float, help="Minimum redshift")
parser.add_argument('z_max', type=float, help="Maximum redshift")

parser.add_argument('--N', type=int, default=1, action='store', 
                    help="Number of lightcones.")

parser.add_argument('-f', '--filename', 
                    action='store', 
                    dest='filename',
                    type=str, 
                    #required=False, 
                    default='output/halo_lightcone_%03d.h5', 
                    help="Output filename (Full path from current directory. \
                          Must include string formatting to iterate over lightcone number)")

parser.add_argument("--overwrite", help="Overwrite existing datasets", 
                    action='store_true')

parser.add_argument("--halos", help="Output halos, rather than galaxies.", 
                    action='store_true')

parser.add_argument("--mass", help="Save halo/galaxy mass information.", 
                    action='store_true')

parser.add_argument("--sfr", help="Save halo/galaxy star formation rate information.", 
                    action='store_true')

parser.add_argument("--midsnap", help="Defines how snaps are arranged in the lightcone. \
                                        By default they are arranged from the leading edge, so the z = 0 snapshot is used first (for a lightcone starting at z = 0).\
                                        If mid-snap is set each snap is represented by the mid-point of the redshift range, \
                                        so the second output snapshot is used from z = 0 (for a lightcone starting at z=0).", 
                    action='store_true')

parser.add_argument("-v", "--verbose", help="Increase output verbosity", 
                    action='store_true')

args = parser.parse_args()

verbose = args.verbose
area = args.area
N_lcs = args.N
z_min = args.z_min
z_max = args.z_max
lc_fname = args.filename
overwrite = args.overwrite
midsnap = args.midsnap
halo_flag = args.halos

print("############################")
print("Area: %.2f deg^2"%area)
print("N: %i"%N_lcs)
print("z_min: %.2f"%z_min)
print("z_max: %.2f"%z_max)
print("filename: %s"%lc_fname)
print("Halos?: %d"%halo_flag)
print("############################\n\n")


if N_lcs < 1: raise ValueError('Need at least one lightcone to generate!')

if overwrite: 
    if sb.yesno('Overwrite is set. Are you sure you wish to overwrite any existing files?\
                 (all lightcones will be overwritten)'):
        for _lc in np.arange(N_lcs):
            if os.path.exists(lc_fname%_lc): os.remove(lc_fname%_lc)
    else:
        print("Exiting...")
        sys.exit()

outs = np.loadtxt(sb.output_file)

if midsnap:
    snaps = np.array([str(s).zfill(3) for s in np.arange(0,296,2)[::-1]])
else:
    snaps = np.array([str(s).zfill(3) for s in np.arange(1,296,2)[::-1]])

zeds = np.array([1./outs[int(snap)] - 1 for snap in snaps])

## ---- check area can be covered by box
cs = sb.get_caesar(snaps[0])
L = cs.simulation.boxsize.to('kpccm').value/10**3

zeds_mask = np.where((zeds >= z_min) & (zeds <= z_max))[0]
## ---- extend zed array 
if np.min(zeds_mask) != 0:
    zeds_mask = np.append(np.min(zeds_mask) - 1, zeds_mask)
if np.max(zeds_mask) < (len(snaps)-2):
    zeds_mask = np.append(zeds_mask, np.max(zeds_mask) + 1)
if np.max(zeds_mask) > (len(snaps)-2):
    zeds_mask = zeds_mask[zeds_mask < (len(snaps) - 1)] 


L_unit = sb.cosmo.kpc_comoving_per_arcmin(zeds[zeds_mask]).to('Mpc / degree')
A = (L_unit * area**0.5).value
print(np.max(np.asarray(L_unit)), np.max(np.asarray(A)), L, area)

# if np.any(A > L):
#     raise ValueError('Specified area too large for simulation box (lateral tiling not yet implemented)')

## ---- start lightcone creation
#lc_out = {str(_lc): {} for _lc in np.arange(N_lcs)}
lc_out = {str(_lc): {} for _lc in np.linspace(N_lcs, N_lcs, 1)}

_N_all = 0
A_A = 0.
for i,snapA in enumerate(snaps[zeds_mask]):
    
    z   = 1./outs[int(snapA)] - 1
    z_B = 1./outs[int(snapA) - 2] - 1
    print("\nz:",z,snapA)

    cs = sb.get_caesar(snapA)
    a = cs.simulation.scale_factor    

    L_unit = sb.cosmo.kpc_comoving_per_arcmin(z_B).to('Mpc / degree')
    A = (L_unit * area**0.5).value

    z_offset = sb.cosmo.comoving_distance(z).value
    if midsnap:
        z_offset = sb.cosmo.comoving_distance(1./outs[int(snapA)+1]-1).value
        # z_offset -= (L/2)

    if verbose: print("z_offset: %.2f"%z_offset)
    
    if halo_flag:
        coods_cMpc = np.array([h.pos.to('kpccm').value/10**3 for h in cs.halos])
    else:
        coods_cMpc = np.array([g.pos.to('kpccm').value/10**3 for g in cs.galaxies])
   
    if args.mass:
        if halo_flag:
            mass = np.array([h.mass.value for h in cs.halos])
        else:
            mass = np.array([g.mass.value for g in cs.galaxies])
    
    if args.sfr:
        if halo_flag:
            sfr = np.array([h.sfr.value for h in cs.halos])
        else:
            sfr = np.array([g.sfr.value for g in cs.galaxies])

    if len(coods_cMpc) == 0:
        print("No galaxies left! Exiting...")
        break

    _N_all += len(coods_cMpc)
    
    if verbose: print("Generating lightcone selection...")

  #  for _lc in np.arange(N_lcs):
    for _lc in np.linspace(N_lcs, N_lcs, 1):
        lc_out[str(_lc)][snapA] = {}
        if verbose: print("N lightcone:", _lc)

        i = np.random.randint(0,3)  # randomly choose axes
        j = i
        while j == i: j = np.random.randint(0,3)
        k = np.where((np.arange(0,3) != i) & (np.arange(0,3) != j))[0][0]
    
        xmin,ymin = np.random.rand(2) * (L-A)   # get minimum box mask coordinate
        if verbose: print("xmin:",xmin, "\nymin:",ymin, "\nA:",A, "\nL:",L)
  
        # lightcone 'frustum' angle
        theta = np.arctan((A - A_A) / (2*L))
        # dx == distance along cood x between top and bottom of frustum edge
        dx = np.abs(L - coods_cMpc[:,k]) * np.tan(theta)

        lc_idx_arr = ((coods_cMpc[:,i] > (xmin + dx)) &\
                      (coods_cMpc[:,i] < ((xmin+A) - dx)) &\
                      (coods_cMpc[:,j] > (ymin + dx)) &\
                      (coods_cMpc[:,j] < ((ymin+A) - dx)))
        
        if verbose: print("N(lightcone cut):",np.sum(lc_idx_arr))
    
        lc_idx_arr = np.where(lc_idx_arr)[0]        
        _coods = coods_cMpc[lc_idx_arr]

        _frac = np.abs(_coods[:,i] - xmin - (A/2)) / ((A/2) - dx[lc_idx_arr])
        _coods[:,i] = _frac * ((A * u.Mpc) / L_unit)
        
        _frac = np.abs(_coods[:,j] - ymin - (A/2)) / ((A/2) - dx[lc_idx_arr])
        _coods[:,j] = _frac * ((A * u.Mpc) / L_unit)

        _coods = _coods[:,[i,j,k]]

        _ra  = _coods[:,0]
        _dec = _coods[:,1]
        _redshift = np.array([z_at_value(sb.cosmo.comoving_distance, _c * u.Mpc) \
                              for _c in (_coods[:,2] + z_offset)])

        ## re-filter by _redshift
        redshift_mask = (_redshift > z_min) & (_redshift < z_max)
        _ra = _ra[redshift_mask]
        _dec = _dec[redshift_mask]
        _redshift = _redshift[redshift_mask]
        lc_idx_arr = lc_idx_arr[redshift_mask]

        with h5py.File(lc_fname%_lc,'a') as h5file:
            h5file.require_group(snapA)

        sb.create_dataset(lc_fname%_lc, lc_idx_arr, 'index', group=snapA, overwrite=True)
        sb.create_dataset(lc_fname%_lc, _ra, 'RA', group=snapA, overwrite=True)
        sb.create_dataset(lc_fname%_lc, _dec, 'DEC', group=snapA, overwrite=True)
        sb.create_dataset(lc_fname%_lc, _redshift, 'z', group=snapA, overwrite=True)

        if args.mass: 
            sb.create_dataset(lc_fname%_lc, mass[lc_idx_arr], 'Stellar Mass', group=snapA, overwrite=True)
        
        if args.sfr: 
            sb.create_dataset(lc_fname%_lc, sfr[lc_idx_arr], 'SFR', group=snapA, overwrite=True)
        
        with h5py.File(lc_fname%_lc,'a') as h5file:
            h5file.require_group(snapA)


    A_A = A


print("_N_all:",_N_all)

