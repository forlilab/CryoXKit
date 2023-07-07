#!/usr/bin/env python

import argparse
from cryo2grid import Cryo2Grid

parser = argparse.ArgumentParser(description='Convert CryoEM/X-ray maps to AutoDock grid maps.')
parser.add_argument('-m', '--map-file', nargs='+', type=str, help='Density map file name', required=True)
parser.add_argument('-r', '--map-receptor', nargs='+', type=str, help='Receptor associated with density map')
parser.add_argument('-g', '--grid-maps', nargs='+', type=str, help='Grid map file name(s)', required=True)
parser.add_argument('-l', '--log_min', type=float, help='Min for logistic curve', default=-3.0)
parser.add_argument('-k', '--log_slope', type=float, help='slope for logistic curve', default=2.0)
parser.add_argument('-x0', '--log_x0', type=float, help='Inflection point for logistic curve', default=0.5)
args = parser.parse_args()

c2g = Cryo2Grid(args.map_file, args.map_receptor)
grid_maps = c2g.ReadGridMaps(args.grid_maps)
density   = c2g.ReadMapFiles()
modified  = c2g.ModifyDensity(density, args.log_min, args.log_slope, args.log_x0)
c2g.WriteGridMaps(modified, 1)

