#!/usr/bin/env python

import argparse
from cryo2grid import Cryo2Grid

parser = argparse.ArgumentParser(description='Convert CryoEM/X-ray maps to AutoDock grid maps.')
parser.add_argument('-m', '--map-file', type=str, help='Map file name', required=True)
parser.add_argument('-g', '--grid-maps', nargs='+', type=str, help='Grid map file name(s)', required=True)
args = parser.parse_args()

c2g = Cryo2Grid(args.map_file)
grid_maps = c2g.ReadGridMaps(args.grid_maps)
density   = c2g.ReadMapFile()
modified  = c2g.ModifyDensity(density)
c2g.WriteGridMaps(modified, grid_maps, args.grid_maps, 1)

