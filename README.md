# FM's FireFOAM
This is the public-facing release of the Fire Modeling team at FM's customized version of the FireFOAM solver for CFD simulations of fires. 
The key external dependencies are [ESI/OpenFOAM](https://develop.openfoam.com/Development/openfoam)-v2306 and [ESI/ThirdParty](https://sourceforge.net/projects/openfoam/files/)-v2306. See also [our company website](https://www.fm.com/about-us/our-engineering-approach/engineering-methods/open-source-fire-modeling).

# Requirements
- Recent Linux distribution (e.g., Ubuntu 20.04)
- GNU compilers (gcc, g++)
- flex
- bison
- cmake

For further details, see the [OpenFOAM documentation](https://www.openfoam.com/documentation/system-requirements)

# Installation instructions
This approach will yield an 'out of source' build in 'Opt' mode (i.e. performance-optimized, no debug flags).
Run the following commands on the command line:
```
mkdir FM-FireFOAM
cd FM-FireFOAM
git clone https://github.com/fmglobal/fireFoam.git FireFOAM-v2306
cd FireFOAM-v2306
chmod +x installation.bash
cd ..
./FireFOAM-v2306/installation.bash
```

This will yield the following structure:

```
FM-FireFOAM/
  FireFOAM-v2306/
  OpenFOAM-v2306/
  ThirdParty-v2306/
```

The `fireFoam` binary file will be located in the default location:
```
$HOME/OpenFOAM/$USER-v$VERSION/platforms/<operating system><compiler>DPInt32<compilation type>/bin/fireFoam
```

For more installation options, run
```
./FireFOAM-v2306/installation.bash --help
```

# Environment setup instructions
FireFOAM is executed by invoking the binary executable `fireFoam`. To ensure it and its dynamically-linked dependencies are available, you must set up your environment.
Doing so is achieved simply by means of sourcing the pre-defined script
```
source /path/to/FM-FireFOAM/OpenFOAM-v2306/etc/bashrc
```

# Recompiling instructions
You should only ever have to run the installation script once. If you make changes to the FireFOAM source code, you only need to recompile that part of the code, as follows:

```
source /path/to/FM-FireFOAM/OpenFOAM-v2306/etc/bashrc
cd /path/to/FM-FireFOAM/FireFOAM-v2306/build
```
If you made major changes, probably best to first run
```
make clean
```
Otherwise, to  recompile FireFOAM, run
```
cmake .. && make -j<number of cores> && make install
```

Note that if you want to have debugging symbols in FireFOAM, but don't need them in OpenFOAM, you can compile FireFOAM in DEBUG mode by instead running
```
cmake -DCMAKE_BUILD_TYPE=DEBUG .. && make -j<number of cores> && make install
```

If you wish to incorporate changes made by ESI to the OpenFOAM branch you are tracking (default: `maintenance-v2306`) that have been made after you first installed FireFOAM, you can run
```
cd /path/to/FM-FireFOAM/OpenFOAM-v2306
git pull
cd ..
./FireFOAM-v2306/installation.bash --start-step=6
```

# Tutorial cases
There are a number of tutorial cases provided with this repository under the `tutorials` directory. For details of the cases, see the `README.md` file in each tutorial's sub-directory.

Each tutorial can be run in serial computation mode by building the mesh and then invoking the solver, as follows
```
source /path/to/FM-FireFOAM/OpenFOAM-v2306/etc/bashrc
cd /path/to/FM-FireFOAM/FireFOAM-v2306/tutorials/<pick a tutorial case>
./Allrun.pre && fireFoam
```

To run a tutorial case in parallel mode, you must build the mesh, decompose the domain and then invoke the solver, as follows
```
source /path/to/FM-FireFOAM/OpenFOAM-v2306/etc/bashrc
cd /path/to/FM-FireFOAM/FireFOAM-v2306/tutorials/<pick a tutorial case>
./Allrun.pre
/path/to/parallel_config.sh --num-nodes=<number of nodes> --cores-per-node=<number of physical cores per node>
mpirun -np <number of tasks> redistributePar -decompose -constant -allRegions -overwrite -parallel >& log.decomposePar.allRegions
mpirun -np <number of tasks> fireFoam -parallel >& log.fireFoam
```
- You should choose `<number of tasks>` and `<number of nodes>` such that `<number of tasks>/<number of nodes> = <number of physical cores per node>` on your machine. 
- You may choose to invoke the above `mpirun` commands in a script to be run with a scheduler, such as SLURM, prepended with appropriate scheduler configuration commands.
In this case you should modify the `-np` option in the `mpirun` commands as `-np $SLURM_NTASKS`.

The script `parallel_config.sh` is available in the top-level tutorial directory. To use it, you must make it executable by running
```
chmod +x /path/to/parallel_config.bash
```
This script will configure the FireFOAM case files and any existing SLURM job submission files appropriately. For more information on using this script, run
```
/path/to/parallel_config.bash --help
```

# Maintainers
In last-name alphabetical order:
- Alex Krisman (alex.krisman@fmglobal.com)
- Xiaoyi Lu (xiaoyi.lu@fmglobal.com)
- Danyal Mohaddes (danyal.mohaddes@fmglobal.com)
- Ning Ren (ning.ren@fmglobal.com)
- Yi Wang (yi.wang@fmglobal.com)

# Contributors
In addition to the Maintainers, in last-name alphabetical order:
- Prateep Chatterjee
- Ankur Gupta
- Karl Meredith
- Oluwayemisi Oluwole
