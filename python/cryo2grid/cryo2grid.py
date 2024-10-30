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
    def __init__(self, map_files=None, map_receptors=None, rmsd_cutoff=-1.0, gaussian_filter_sigma=0, noise_std_range=0, repeat_unit_cell=True, output_align_rec=False):
        """Initialize the Cryo2Grid object.

        Arguments:
            map_file (str): Map (X-ray or CryoEM) file name

        """
        print_version_info()
        
        self._map_files             = map_files
        self._map_receptors         = map_receptors
        self._repeat_unit           = repeat_unit_cell
        self._out_align_rec         = output_align_rec
        self._align_rec             = ''
        self._map_x_dim             = None
        self._map_y_dim             = None
        self._map_z_dim             = None
        self._map_x_center          = None
        self._map_y_center          = None
        self._map_z_center          = None
        self._grid_spacing          = None
        self._grid_files            = None
        self._gridmaps              = None
        self._rmsd_cutoff           = rmsd_cutoff
        self._gaussian_filter_sigma = gaussian_filter_sigma
        self._noise_std_range       = noise_std_range
    
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
    
    def Convert2MRC(self, mapfile):
        convert_map_to_mrc(mapfile);
    
    def WriteDensity(self, density, basename, write_type=2):
        """Write AD4 grid maps with density added

        Args:
            density: density to add to grid maps (usually the modified density)
            basename: base filename
            write_type: type of grid map file to write (0 .. nothing, 1 .. AD4, 2 .. MRC)

        """
        write_density(density, basename, write_type);
    
    def WriteGridMaps(self, density, write_type=1):
        """Write AD4 grid maps with density added

        Args:
            density: density to add to grid maps (usually the modified density)
            write_type: type of grid map file to write (0 .. nothing, 1 .. AD4, 2 .. MRC)

        """
        if self._grid_files is None:
            raise RuntimeError('ERROR: No grid maps have been read yet.')
        write_grid_maps(density, self._gridmaps, self._grid_files, write_type);
    
    def ReadMapFiles(self, map_files=None, map_receptors=None, align_rec=None, map_x_dim=None, map_y_dim=None, map_z_dim=None, map_x_center=None, map_y_center=None, map_z_center=None, grid_spacing=None, rmsd_cutoff=None, gaussian_filter_sigma=None, noise_std_range=None):
        """Read Map (X-ray or CryoEM) file(s)

        Arguments:
            map_files (str list): Density map filename(s) (code autodetects CCP4/MRC and DSN6/Brix formats)
            map_receptors (str list): Density map receptor filename(s) (if empty no alignement)
            align_rec (str): Receptor to align density map receptor(s) to
            map_{x,y,z}_dim (int): Number of grid points in x,y,z-direction
            map_{x,y,z}_center (float): center vector of grid map
            grid_spacing (float): grid map spacing
            rmsd_cutoff (float): rmsd cutoff for aligments
            gaussian_filter_sigma (float): width of gaussian filter to be applied to all density maps (0 means no filter)
            noise_std_range (float): width as a fraction of std.dev. of noise to add (negative means add before gaussian convolution)
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
            if(self._map_x_dim) is None:
                raise RuntimeError('ERROR: Number of grid points in x-direction not set.')
            map_x_dim=int(self._map_x_dim)
        if map_y_dim is None:
            if(self._map_y_dim) is None:
                raise RuntimeError('ERROR: Number of grid points in y-direction not set.')
            map_y_dim=int(self._map_y_dim)
        if map_z_dim is None:
            if(self._map_z_dim) is None:
                raise RuntimeError('ERROR: Number of grid points in z-direction not set.')
            map_z_dim=int(self._map_z_dim)
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
        if rmsd_cutoff is None:
            rmsd_cutoff=self._rmsd_cutoff
        if gaussian_filter_sigma is None:
            gaussian_filter_sigma=self._gaussian_filter_sigma
        if noise_std_range is None:
            noise_std_range=self._noise_std_range
        density = average_densities_to_grid(map_files, map_receptors, align_rec, 0, map_x_dim, map_y_dim, map_z_dim, map_x_center, map_y_center, map_z_center, grid_spacing, rmsd_cutoff, gaussian_filter_sigma, noise_std_range, self._repeat_unit, self._out_align_rec)
        return density;

    def ModifyDensity(self, density, log_max=-3.0, width=-2.0, x0=-1.0):
        """Modify map density using logistics function

        Args:
            density: density to add to grid maps (read using ReadMapFiles above)
            log_max: reward value of logistics function                (default: -3 kcal/mol)
            width:   width of logistics function                       (default: -2 == use 10%*2*std.dev. of data)
            x0:      fraction determining center of logistics function (default: -1 == use median of data)

        """
        modified = modify_densities(density, 1, log_max, width, x0)
        return modified;

    def CreateMask(self, grid_or_mask, mask_pdb, rT = 2.0, subtractive=True, create_new=True):
        """Create a mask to modify map densities

        Args:
            grid_or_mask: grid maps to define grid dimensions *or* existing mask to add to
            mask_pdb:     PDB filename to use as a mask
            rT:           radius to use for atom density mask
            subtractive:  atoms in mask_pdb are used to subtract from the map
            create_new:   start a new mask or add to an existing one

        """
        if mask_pdb is None:
            return None
        mask = create_mask(grid_or_mask, mask_pdb, rT, subtractive, create_new)
        return mask

    def ApplyMask(self, density, mask):
        """Apply a multiplicative mask to a map density

        Args:
            density: density map
            mask:    mask

        """
        if mask is None:
            return density
        masked_density = apply_mask(density, mask);
        return masked_density

