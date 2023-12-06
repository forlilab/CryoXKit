#!/usr/bin/env python

import argparse
from cryo2grid import Cryo2Grid

parser = argparse.ArgumentParser(description='Convert CryoEM/X-ray maps to AutoDock grid maps.')
parser.add_argument('-m', '--map-file', nargs='+', type=str, help='Density map file name', required=True)
parser.add_argument('-r', '--map-receptor', nargs='+', type=str, help='Receptor associated with density map')
parser.add_argument('-g', '--grid-maps', nargs='+', type=str, help='Grid map file name(s)', required=True)
parser.add_argument('-l', '--log_min', type=float, help='Min for logistic curve', default=-3.0)
parser.add_argument('-w', '--log_width', type=float, help='width for logistic curve', default=4.0)
parser.add_argument('-x0', '--log_x0', type=float, help='relative inflection point for logistic curve (<0 set automatically based on data median)', default=-1)
parser.add_argument('-p', '--repeat_unit_cell', type=bool, help='Allow density map unit cell repetition', default=True)
parser.add_argument('-a', '--output_aligned_receptors', type=bool, help='Output aligned receptors as pdb for inspection', default=False)
parser.add_argument('-s', '--mask-file', type=str, help='PDB file to use to mask density')
args = parser.parse_args()

c2g = Cryo2Grid(args.map_file, args.map_receptor, args.repeat_unit_cell, args.output_aligned_receptors)
grid_maps = c2g.ReadGridMaps(args.grid_maps)
density   = c2g.ReadMapFiles()
mask      = c2g.CreateMask(grid_maps, args.mask_file)
c2g.ApplyMask(density, mask)
modified  = c2g.ModifyDensity(density, args.log_min, args.log_width, args.log_x0)
c2g.WriteGridMaps(modified, 1)
