# CryoXKit

A tool to (currently) read DSN6, BRIX, and MRC/CCP4 crystallography and cryo-EM files and interpolate their density values onto a cartesian grid map as used by AutoDock software. The tool can output AD4 grid maps (default) or grid MRC files.

# Build instructions

On most platforms only "make" should be needed. On macOS the compiler in the Makefile needs to be changed to `clang++` to use Brew's version supporting OpenMP - otherwise `-fopenmp` needs to be taken out of the `CXXFLAGS` in the Makefile.

For direct Python compilation and installation (as above, we recommend LLVM from HomeBrew for macOS), from within your favorite python environment run `make python`.

# Syntax

For the installed Python version please run `./run_cryoXkit.py`.

The standalone terminal tool CryoXKit needs parameters for the density map file and AD4 grid maps to modify *or* grid map center x, y, z coordinates, grid map x, y, z dimensions, and optionally grid spacing (default is 0.375 A), type of map to write (0 .. nothing, 1 .. AD4 grid map, 2 .. MRC grid map), as well as the modifier function type (0 .. none, 1 .. logistics) when only the modifier map is of interest:
```
./cryoXkit mapfile center_x center_y center_z x_dim y_dim z_dim (spacing [0.375]) (write [1 = AD4 map]) (modifier fxn [1 = logistics])
*or* for map modification (automatically excludes e, d, and H* maps):
./cryoXkit mapfile gridfile (spacing [0.375]) (modifier fxn [1 = logistics])
```

Once built to test everything works the following command can be used:

`./cryoXkit example/*.dsn6 42.43 41.664 10.503 40 40 40`

It should produce the grid map file `1d3g_2fofc.map` in the `./example` folder and similar output:
```
Reading map file [example/1d3g_2fofc.dsn6]
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
    -> density unit cell range: (25.0768, 13.5440, -17.7846) A to (54.4678, 69.1212, 40.2769) A
    -> density range: -2.581 to 12.665 (average: -0.029997 +/- 1.000289; median: -0.118788)
<- Finished reading densities, took 10.085 ms.

Interpolating density data for 40x40x40 grid (spacing: 0.375 A)
    -> grid start:  (34.930, 34.164, 3.003) A
    -> grid size:   (15.000, 15.000, 15.000) A
<- Finished interpolating densities to grid map points, took 1.228 ms.

Calculating median
    -> range: -2.156 to 10.769 (median: -0.578)
<- Done, took 0.143 ms.

Adjusting density values using logistics function modifier (l: -3 kcal/mol, w: 2 kcal/mol, x0: 0.122)
<- Finished adjusting, took 0.908 ms.

Writing AD4 grid map file [example/1d3g_2fofc.map]
<- Finished writing, took 11 ms.

Done. Overall runtime was 23.8 ms.
```

For map conversions, the tool only requires the map file and the AD4 grid map files, an example is (note: this will modify the AD4 map files) to run `./cryoXkit example/1ac8/*.dsn6 example/1ac8/*.map`:
```
Reading grid map files:
    -> example/1ac8/protein.A.map
    -> example/1ac8/protein.C.map
    -> example/1ac8/protein.SA.map
    -> example/1ac8/protein.N.map
<- Done, took 36.376 ms.

Reading map file [example/1ac8/1ac8_2fofc.dsn6]
    -> DSN6 endian-swapped, file size: 369152
    -> x_dim = 79, y_dim = 62, z_dim = 69
    -> a_unit = 108, b_unit = 77.3, c_unit = 51.7875
    -> alpha = 90, beta = 90, gamma = 90
    -> Fractional to cartesian conversion matrix:
	 108.0000   -0.0000   -0.0000
	        0   77.3000   -0.0000
	        0         0   51.7875
    -> Cartesian to fractional conversion matrix:
	   0.0093    0.0000    0.0000
	        0    0.0129   -0.0000
	        0         0    0.0193
    -> density unit cell range: (-9.2571, 69.1240, 28.7708) A to (50.9143, 114.4635, 77.6813) A
    -> density range: -3.681 to 11.902 (average: -0.030685 +/- 1.000417; median: -0.121865)
<- Finished reading densities, took 2.301 ms.

Interpolating density data for 60x60x60 grid (spacing: 0.375 A)
    -> grid start:  (20.674, 82.194, 36.674) A
    -> grid size:   (22.500, 22.500, 22.500) A
<- Finished interpolating densities to grid map points, took 3.682 ms.

Calculating median
    -> range: -3.515 to 11.530 (median: -0.293)
<- Done, took 0.532 ms.

Adjusting density values using logistics function modifier (l: -3 kcal/mol, w: 2 kcal/mol, x0: 0.214)
<- Finished adjusting, took 2.79 ms.

Writing AD4 grid map file [example/1ac8/protein.N.map]
Writing AD4 grid map file [example/1ac8/protein.A.map]
Writing AD4 grid map file [example/1ac8/protein.C.map]
Writing AD4 grid map file [example/1ac8/protein.SA.map]

Done. Overall runtime was 88.1 ms.
```
