# Cryo2Grid

A tool to (currently) read DSN6, BRIX, and MRC/CCP4 crystallography and cryo-EM files and interpolate their density values onto a cartesian grid map as used by AutoDock software. The tool can output AD4 grid maps (default) or grid MRC files.

# Build instructions

On most platforms only "make" should be needed. On macOS the compiler in the Makefile may need to be changed to `clang++`.

# Syntax

Cryo2Grid needs parameters for the density map file, grid map center x, y, z coordinates, grid map x, y, z dimensions, and optionally grid spacing (default is 0.375 A), type of map to write (0 .. nothing, 1 .. AD4 grid map, 2 .. MRC grid map), as well as the modifier function type (0 .. none, 1 .. logistics):

`./cryo2grid mapfile center_x center_y center_z x_dim y_dim z_dim (spacing [0.375]) (write [1 = AD4 map]) (modifier fxn [1 = logistics])`

Once built to test everything works the following command can be used:

`./cryo2grid example/*.dsn6 42.43 41.664 10.503 40 40 40`

It should produce the grid map file `1d3g_2fofc.map` in the `./example` folder and similar output:
```
Reading map file example/1d3g_2fofc.dsn6
    -> DSN6 endian-swapped, file size: 1728512
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
    -> density coordinate range: (25.077, 13.544, -17.785) A to (54.468, 69.121, 40.277) A
    -> density range: -2.581 to 12.665 (average: -0.029997 +/- 1.000289)
<- Finished reading densities, took 6.419000 ms.

Interpolating density data for 40x40x40 grid (spacing: 0.375 A)
    -> grid start:  (34.930, 34.164, 3.003) A
    -> grid size:   (15.000, 15.000, 15.000) A
<- Finished interpolating and writing AD4 grid map, took 8.768 ms.

Done. Overall time was 15.202 ms.
```
