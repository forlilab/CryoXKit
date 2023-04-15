# Cryo2Grid

A tool to (currently) read DSN6 (and eventually MRC etc.) crystallography and cryo-EM files and interpolate their density values onto a cartesian grid map as used by AutoDock software.

# Build instructions

On most platforms only "make" should be needed. On macOS the compiler in the Makefile may need to be changed to `clang++`.

# Syntax

Cryo2Grid needs parameters for the density map file, grid map center x, y, z coordinates, grid map x, y, z dimensions, and optionally the grid spacing (default is 0.375 A):

`./cryo2grid <map file> <center x> <center y> <center z> <x dim> <y dim> <z dim> <optional: grid spacing (default: 0.375)>`

Once built to test everything works the following command can be used:

`./cryo2grid example/*.dsn6 42.43 41.664 10.503 40 40 40`

It should produce the grid map file `1d3g_2fofc.map` in the `./example` folder and similar output:
```
Reading DSN6 map file example/1d3g_2fofc.dsn6
    -> file size: 1728512
    -> endian swap: 1, Norm: 100
    -> x_dim = 115, y_dim = 120, z_dim = 112
    -> a_unit = 90.6, b_unit = 90.6, c_unit = 122.4
    -> alpha = 90, beta = 90, gamma = 120
    -> Fractional to cartesian conversion matrix:
	  90.6000  -45.3000   -0.0000
	        0   78.4619   -0.0000
	        0         0  122.4000
    -> Cartesian to fractional conversion matrix:
	   0.0110    0.0064    0.0000
	        0    0.0127   -0.0000
	        0         0    0.0082
    -> density coordinate range: (25.0768, 13.5440, -17.7846) A to (54.4678, 69.1212, 40.2769) A
    -> density range: -2.5810 to 12.6651
<- Finished reading densities, took 7.4070 ms.

Interpolating density data for 40x40x40 grid (spacing: 0.375 A)
    -> grid start:  (34.930, 34.164, 3.003) A
    -> grid size:   (15.000, 15.000, 15.000) A
<- Finished interpolating and writing grid map, took 9.470 ms.

Done. Overall time was 16.891 ms.
```
