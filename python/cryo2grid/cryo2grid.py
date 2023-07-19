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
    def __init__(self, map_files=None, map_receptors=None, repeat_unit_cell=True):
        """Initialize the Cryo2Grid object.

        Arguments:
            map_file (str): Map (X-ray or CryoEM) file name

        """
        self._map_files     = map_files
        self._map_receptors = map_receptors
        self._repeat_unit   = repeat_unit_cell
        self._align_rec     = ''
        self._map_x_dim     = None
        self._map_y_dim     = None
        self._map_z_dim     = None
        self._map_x_center  = None
        self._map_y_center  = None
        self._map_z_center  = None
        self._grid_spacing  = None
        self._grid_files    = None
        self._gridmaps      = None
    
    def ReadGridMaps(self, gridmap_files):
        """Read AD4 grid maps.

        Args:
            gridmap_files (list of strings): Grid map file names

        """
        self._grid_files   = filter_grid_files(gridmap_files)
        self._gridmaps     = read_grid_maps(self._grid_files, self._align_rec)
        self._align_rec    = get_grid_receptor_filename(self._gridmaps, self._grid_files)
        self._map_x_dim    = (self._gridmaps[0])[1]
        self._map_y_dim    = (self._gridmaps[0])[2]
        self._map_z_dim    = (self._gridmaps[0])[3]
        self._map_x_center = (self._gridmaps[0])[4]
        self._map_y_center = (self._gridmaps[0])[5]
        self._map_z_center = (self._gridmaps[0])[6]
        self._grid_spacing = (self._gridmaps[0])[7]
        return self._gridmaps;
    
    def WriteGridMaps(self, density, write_type=1):
        """Read AD4 grid maps.

        Args:
            density: density to add to grid maps (usually the modified density)
            gridmaps: the grid map data as read by the function above
            gridmap_files (list of strings): Grid map file names
            write_type: type of grid map file to write (0 .. nothing, 1 .. AD4, 2 .. MRC)

        """
        if self._grid_files is None:
            raise RuntimeError('ERROR: No grid maps have been read yet.')
        write_grid_maps(density, self._gridmaps, self._grid_files, write_type);
    
    def ReadMapFiles(self, map_files=None, map_receptors=None, align_rec=None, map_x_dim=None, map_y_dim=None, map_z_dim=None, map_x_center=None, map_y_center=None, map_z_center=None, grid_spacing=None):
        """Read Map (X-ray or CryoEM) file(s)

        Arguments:
            map_files (str list): Density map filename(s) (code autodetects CCP4/MRC and DSN6/Brix formats)
            map_receptors (str list): Density map receptor filename(s) (if empty no alignement)
            align_rec (str): Receptor to align density map receptor(s) to
            map_{x,y,z}_dim (int): Number of grid points in x,y,z-direction
            map_{x,y,z}_center (float): center vector of grid map
            grid_spacing(float): grid map spacing
        Notes:
            - the map file name can be set at class initialization already
            - the grid map parameters will automatically be set when reading grid maps (see above) but will have to be specified otherwise

        """
        if map_files is None:
            map_files=self._map_files
        if map_files is None:
            raise RuntimeError('ERROR: No Map file specified.')
        if map_receptors is None:
            map_receptors=self._map_receptors
        if map_receptors is None: # if it's still None, we need to pass an empty string in the list
            map_receptors = ['']
        if align_rec is None:
            align_rec=self._align_rec
        if align_rec is None: # if it's still None, we need to pass an empty string in the list
            align_rec = ''
        if map_x_dim is None:
            map_x_dim=int(self._map_x_dim)
        if map_x_dim is None:
            raise RuntimeError('ERROR: Number of grid points in x-direction not set.')
        if map_y_dim is None:
            map_y_dim=int(self._map_y_dim)
        if map_y_dim is None:
            raise RuntimeError('ERROR: Number of grid points in y-direction not set.')
        if map_z_dim is None:
            map_z_dim=int(self._map_z_dim)
        if map_z_dim is None:
            raise RuntimeError('ERROR: Number of grid points in z-direction not set.')
        if map_x_center is None:
            map_x_center=self._map_x_center
        if map_x_center is None:
            raise RuntimeError('ERROR: Grid map center_x not set.')
        if map_y_center is None:
            map_y_center=self._map_y_center
        if map_y_center is None:
            raise RuntimeError('ERROR: Grid map center_y not set.')
        if map_z_center is None:
            map_z_center=self._map_z_center
        if map_z_center is None:
            raise RuntimeError('ERROR: Grid map center_z not set.')
        if grid_spacing is None:
            grid_spacing=self._grid_spacing
        if grid_spacing is None:
            raise RuntimeError('ERROR: Grid map spacing not set.')
        density = average_densities_to_grid(map_files, map_receptors, align_rec, 0, map_x_dim, map_y_dim, map_z_dim, map_x_center, map_y_center, map_z_center, grid_spacing, self._repeat_unit)
        return density;

    def ModifyDensity(self, density, log_max=-3.0, rate=2.0, x0=0.5):
        """Modify map density using logistics function

        Args:
            density: density to add to grid maps (usually the modified density)
            log_max, rate, x0: parameters for logistics function used to modify density

        """
        modified = modify_densities(density, 1, log_max, rate, x0)
        return modified;

