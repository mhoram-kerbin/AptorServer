AptorServer - A Launch Optimizer for KSP

Copyright 2014 Mhoram Kerbin

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


# Installation on Windows

## Compile it yourself

### Prerequisites

In order to compile AptorServer, you will need a version of PSOPT.
There is a ThirdParty installer available at
https://github.com/Sauermann/psopt-installer

Follow the instructions for the windows-installation. And if you want
to generate output graphs, install GnuPlot.

Download an archive or git-repository of AptorServer and extract it in
a directory.

### Compilation

Compilation happens within MSYS.
1. Start MSYS
2. Set the environment variable PSOPT_DIR to the location of your
   PSOPT installation:
   export PSOPT_DIR="/c/psopt"
   or adjust the path in the file AptorServer Makefile
3. Change to the directory of AptorServer e.g.:
   cd /c/psopt/AptorServer
4. Compile it:
   make
5. run it with a test configuration
   ./AptorServer.exe file config.h

## Binary version

There is a precompiled version of AptorServer available at
https://www.dropbox.com/sh/e9azrbsm5hpp3nq/x3ZKIL-L4r

The binary version differs from the self-compiled variant in the
following aspect:
Instead of the software package Metis, SCOTCH is used, because Metis 4
is not open source and is not allowed to be distributed.

# Usage

AptorServer is started from the command line.
It is configured by text-commands.

Allowed command line options are

AptorServer.exe file NAME

where NAME is the filename to read.
Without any command line arguments, it loads the file config.h and
performs the actions specified within.


# Configuration

The following case sensitive commands are recognized. Empty lines are
ignored and lines starting with a \# are considered to be comments.

## Configuration of the physical surroundings

*PLANET_MASS 5.2915793E22*

Set the Mass of the planet in kg

*PLANET_RADIUS 600000*

Set the planet Radius in m

*PLANET_SCALE_HEIGHT 5000*

Set the scale height in m

*PLANET_P0 1*

Set the density at sealevel in atm

*PLANET_ROTATION_PERIOD 21600*

Set the duration of one planet rotation in s

*PLANET_SOI 84159286*

Set the radius of SOI in m

*LAUNCH_LATITUDE -0.001691999747072392*

Set the launch latitude in rad

*LAUNCH_LONGITUDE 0*

Set the launch longitude in rad (please use 0 at the moment)

*LAUNCH_ALTITUDE 77.6*

Set the launch altitude above sealevel in m

*MAX_VELOCITY 10000*

Set the maximal velocity of a ship in m/s (leave it at this value)

*NAME Tangent_2*

Set the name of the ship (may not contain whitespaces)

*TARGET_PERIAPSIS 75000*

Set the target periapsis in m

*TIME_TO_ORBIT 360*

Set a guess of how long the ascent will last in s (this command is optional)


## Configuration of the rocket

The order of the commands in this section is relevant.

Stages have to be specified by order of decreasing mass, so that you
begin with the first stage that is jettisoned.
ADD_ENGINE adds an engine to the stage associated with the previous
ADD_STAGE command.
All stages must contain at least one engine.


*ADD_STAGE 100.1 32 0.2*

Add a stage with initial mass (including fuel) of 100.1 ton, fuel mass
(including oxidizer) of 32 ton and drag coefficient of 0.2

*ADD_ENGINE 175 388 390*

Add an engine to the stage corresponding to the previous ADD_STAGE
command with a thrust of 175 kN, Isp0 of 388 and IspV of 390.

*ADD_ENGINE 175 388 390 4*

Add 4 engines to the current stage with the given stats.

## Configuration of PSOPT

These values require a bit of fiddling. It took me some time to adjust
them in a way to get reasonable results. For a more detailled
description, read the PSOPT manual.

For a 2-3 stage rocket with reasonable TWR, the following setting
should work relatively fine

-ITERATIONS 500
-SET_NODES 10 30
-MESH_REFINEMENT manual
-NLP_TOLERANCE 1.0e-5

For other configurations try and error is the way to go.

*ITERATIONS 500*

Set the number of refinement iterations to 500

*SET_NODES 10 40 ...*

Set the number of mesh nodes for each refinement. (first refinement:
10, second refinement 40, ...)

*MESH_REFINEMENT manual*

Set the mesh refinement method (manual or automatic)

*NLP_TOLERANCE 1.0e-5*

Set the error tolerance to 1.0e-5.




## Perform the Optimization

*COMPUTE*

Run Psopt with the current configuration (this will take some time)

*POSTPROCESS*
*VERBOSE_POSTPROCESS*

Create png and pdf graphics of the ascent. (requires GnuPlot)
The verbose version also displays them on screen.



# Known issues

- Not tested with other planets than Kerbin (Eve most likely will not work)
- Only Launch longitudes of value 0 are supported at the moment
