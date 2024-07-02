#!/usr/bin/env python

import argparse
from cryo2grid import Cryo2Grid

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('-C', '--convert2mrc', type=str, help='Convert density map to mrc file name')
args = parser.parse_known_args()
if args[0].convert2mrc is not None:
    c2g = Cryo2Grid()
    c2g.Convert2MRC(args[0].convert2mrc)
    exit(0)

parser = argparse.ArgumentParser(description='Convert CryoEM/X-ray maps to AutoDock grid maps.', add_help=True)
parser.add_argument('-C', '--convert2mrc', type=str, help='Convert density map to mrc file name')
parser.add_argument('-m', '--map-file', nargs='+', type=str, help='Density map file name', required=True)
parser.add_argument('-r', '--map-receptor', nargs='+', type=str, help='Receptor associated with density map')
parser.add_argument('-g', '--grid-maps', nargs='+', type=str, help='Grid map file name(s)', required=True)
parser.add_argument('-b', '--gaussian_filter_sigma', type=float, help='Width of Gaussian filter in Angstrom (0 means no filter)', default=0)
parser.add_argument('-n', '--noise_std_range', type=float, help='Fraction of std.dev. of density map to add as normal-distributed noise (0 means no filter)', default=0)
parser.add_argument('-l', '--log_min', type=float, help='Min for logistic curve', default=-3.0)
parser.add_argument('-w', '--log_width', type=float, help='width for logistic curve (<0 multiplier of data std.dev.)', default=-4)
parser.add_argument('-x0', '--log_x0', type=float, help='relative inflection point for logistic curve (<0 set automatically based on data median)', default=-1)
parser.add_argument('-c', '--rmsd_cutoff', type=float, help='RMSD cutoff for automatic alignment', default=2.0)
parser.add_argument('-s', '--mask-file', type=str, help='PDB file to use to mask density')
parser.add_argument('-p', '--repeat_unit_cell', type=bool, help='Allow density map unit cell repetition', default=True)
parser.add_argument('-a', '--output_aligned_receptors', help='Output aligned receptors as pdb for inspection', action='store_true')
parser.add_argument('-i', '--output_intermediate_densities', help='Output intermediate densities for inspection', action='store_true')
args = parser.parse_args()

c2g = Cryo2Grid(args.map_file, args.map_receptor, args.rmsd_cutoff, args.gaussian_filter_sigma, args.noise_std_range, args.repeat_unit_cell, args.output_aligned_receptors)
grid_maps      = c2g.ReadGridMaps(args.grid_maps)
density        = c2g.ReadMapFiles()
if args.output_intermediate_densities is True:
    c2g.WriteDensity(density, "read_density")
mask           = c2g.CreateMask(density, args.mask_file)
masked_density = c2g.ApplyMask(density, mask)
if args.output_intermediate_densities is True and args.mask_file is not None:
    c2g.WriteDensity(masked_density, "masked_density")
modified       = c2g.ModifyDensity(masked_density, args.log_min, args.log_width, args.log_x0)
if args.output_intermediate_densities is True:
    c2g.WriteDensity(modified, "modified_density")
c2g.WriteGridMaps(modified, 1)
