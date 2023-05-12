#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Cryo2Grid
#

import os
import glob
import stat
import sys

from .cryo2grid_wrapper import *

class Cryo2Grid:
    def __init__(self, map_file=None):
        """Initialize the Cryo2Grid object.

        Arguments:
            map_file (str): Map (X-ray or CryoEM) file name

        """
        self._map_file     = map_file
        self._map_x_dim    = None
        self._map_y_dim    = None
        self._map_z_dim    = None
        self._map_x_center = None
        self._map_y_center = None
        self._map_z_center = None
        self._grid_spacing = None
    
    def ReadGridMaps(self, gridmap_files):
        """Read AD4 grid maps.

        Args:
            gridmap_files (list of strings): Grid map file names

        """
        gridmaps = read_grid_maps(gridmap_files);
        self._map_x_dim    = (gridmaps[0])[0]
        self._map_y_dim    = (gridmaps[0])[1]
        self._map_z_dim    = (gridmaps[0])[2]
        self._map_x_center = (gridmaps[0])[3]
        self._map_y_center = (gridmaps[0])[4]
        self._map_z_center = (gridmaps[0])[5]
        self._grid_spacing = (gridmaps[0])[6]
        return gridmaps;
    
    def WriteGridMaps(self, density, gridmaps, gridmap_files, write_type=1):
        """Read AD4 grid maps.

        Args:
            density: density to add to grid maps (usually the modified density)
            gridmaps: the grid map data as read by the function above
            gridmap_files (list of strings): Grid map file names
            write_type: type of grid map file to write (0 .. nothing, 1 .. AD4, 2 .. MRC)

        """
        write_grid_maps(density, gridmaps, gridmap_files, write_type);
    
    def ReadMapFile(self, map_file=None, map_x_dim=None, map_y_dim=None, map_z_dim=None, map_x_center=None, map_y_center=None, map_z_center=None, grid_spacing=None):
        """Read Map (X-ray or CryoEM) file

        Arguments:
            map_file (str): Map file name (code autodetects CCP4/MRC and DSN6/Brix formats)
            map_{x,y,z}_dim (int): Number of grid points in x,y,z-direction
            map_{x,y,z}_center (float): center vector of grid map
            grid_spacing(float): grid map spacing
        Notes:
            - the map file name can be set at class initialization already
            - the grid map parameters will automatically be set when reading grid maps (see above) but will have to be specified otherwise

        """
        if map_file is None:
            map_file=self._map_file
        if map_file is None:
            raise RuntimeError('Error: No Map file specified.')
        if map_x_dim is None:
            map_x_dim=int(self._map_x_dim)
        if map_x_dim is None:
            raise RuntimeError('Error: Number of grid points in x-direction not set.')
        if map_y_dim is None:
            map_y_dim=int(self._map_y_dim)
        if map_y_dim is None:
            raise RuntimeError('Error: Number of grid points in y-direction not set.')
        if map_z_dim is None:
            map_z_dim=int(self._map_z_dim)
        if map_z_dim is None:
            raise RuntimeError('Error: Number of grid points in z-direction not set.')
        if map_x_center is None:
            map_x_center=self._map_x_center
        if map_x_center is None:
            raise RuntimeError('Error: Grid map center_x not set.')
        if map_y_center is None:
            map_y_center=self._map_y_center
        if map_y_center is None:
            raise RuntimeError('Error: Grid map center_y not set.')
        if map_z_center is None:
            map_z_center=self._map_z_center
        if map_z_center is None:
            raise RuntimeError('Error: Grid map center_z not set.')
        if grid_spacing is None:
            grid_spacing=self._grid_spacing
        if grid_spacing is None:
            raise RuntimeError('Error: Grid map spacing not set.')
        
        density = read_map_to_grid(map_file, 0, map_x_dim, map_y_dim, map_z_dim, map_x_center, map_y_center, map_z_center, grid_spacing)
        return density;

    def ModifyDensity(self, density, log_max=-3.0, rate=2.0, x0=0.5):
        """Modify map density using logistics function

        Args:
            density: density to add to grid maps (usually the modified density)
            log_max, rate, x0: parameters for logistics function used to modify density

        """
        modified = modify_densities(density, 1, log_max, rate, x0)
        return modified;

