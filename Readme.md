Protein Docking with Discriminative Learning v.2 (PRODDL2)
==========================================================

This is a protein-protein docking package. Its main operation is to accept 3D structures
of two proteins that are known to interact and generate several models of their
3D complex.

This package incorporates codes originally developed by the author while working 
on a docking software GRAMM-X. 

The architecture has been redesigned for better portability and parallel execution 
on a wide range of distributed backends including high-throughput batch clusters, 
multicore workstations, MPI clusters as well as collections of virtual machine 
instances in the clouds with or without shared file systems.

Author
------

Andrey Tovchigrechko <andreyto AT gmail.com>

References
----------

See original GRAMM-X publications regarding the force field and search procedure:

 - Tovchigrechko A, Vakser IA. GRAMM-X public web server for protein-protein 
   docking. Nucleic Acids Res. 2006; 34:W310-4.
 - Tovchigrechko A, Vakser IA. Development and testing of an automated 
   approach to protein docking. Proteins 2005;60(2):296-301.

Regarding portable distributed execution techniques used in this version:

 - Tovchigrechko A, Venepally P, Payne SH. PGP: Parallel Prokaryotic Proteogenomics 
   Pipeline for MPI clusters, high-throughput batch clusters and multicore workstations.
   Bioinformatics 2014.

Funding 
-------

National Science Foundation award 1048199 and allocation on Azure cloud from
Microsoft. Any opinions, findings, and conclusions or recommendations 
expressed in this material are those of the author(s) and do not necessarily 
reflect the views of the National Science Foundation or Microsoft Co.

License
-------

GPLv3. See also COPYING file that accompanies the source code.

Installation
------------

PRODDL2 uses CMake (>2.8.10) for configuration and building.

In `config/<linux|win>`, edit toolchain.cmake and shell script files to help CMake find and use
PRODDL dependencies. 

 - Build-time dependencies: 
   
   Python (>=2.7), FFTW3, Boost (>1.50), HDF5 (>1.8.3), Blitz++ (>0.10).
   PRODDL build procedure will also try to build several Python packages. One of them, Numpy (>1.8)
   might fail to build automatically depending on your Linux flavor. In that case you will have to
   install Numpy before building PRODDL.
 
 - Run-time dependencies: 
   
   [IMP](https://salilab.org/imp/) (>2.2.0). There is a sample shell script 
   `config/linux/deps_build/` for building IMP on Linux. On Windows, we recommend 
   installing a pre-compiled IMP package.

   [Makeflow](http://www3.nd.edu/~ccl/software/makeflow/) (>4.2.0). On Windows, use 
   Cygwin pre-compiled binary. We recommend placing `makeflow` executable in the 
   PATH defined in the shell scripts under `config`. Alternatively,
   you can specify the full path for the executable in `config/noarch/proddl.json.in`. 

 - If you will be using batch queuing system for distributed execution on Linux (SGE, PBS, Slurm),
   you might have to edit `config/linux/proddl-deps.rc.in` to `source` your `~/.bashrc`. This
   will make sure that any environment settings from your `~/.bashrc`, including LD_LIBRARY_PATH, 
   are propagated to the cluster jobs submitted by PRODDL workflow.

On Linux, we used GCC (>4.7.0). On Windows, we used Visual Studio 2010 in 32 bit mode.

After you have configured the dependencies, use sample commands for configuring, building and installing 
PRODDL from `config/linux/build.sh` and `config/win/build.bat`

After successfully building **and** installing PRODDL, run `ctest` in the build directory 
to execute the unit tests. 

Using
-----

Assuming that you have installed PRODDL under PRODDL_HOME, you can run:

```
PRODDL_HOME/bin/proddl-wrapper proddl-dock dock \
    PRODDL_HOME/test_data/pdb/2ptc_E.pdb \
    PRODDL_HOME/test_data/pdb/2ptc_I.pdb \    
    res.pdb
```

The above will dock proteins from two input files and generate 
an output file res.pdb with 10 top ranked models in a format that
is used for multiple NMR structures in PDB (structures are
separated by MODEL/ENDMDL records).

In that simple form, the program will run locally, using all
available CPU cores for the computationally intensive parts
of the workflow.

To run on distributed backends, add a parameter `--makeflow-args` 
to the command line. Its value should be a **quoted** string 
with all necessary Makeflow parameters.

For example, to use SGE, you can add to the command above 
something like:

```
--makeflow-args "-T sge -B '-b n -P 1111'"
```

Notice the use of single quotes inside double quotes above,
done in order to properly escape the native SGE arguments 
supplied with the `-B` argument of Makeflow.

Run `PRODDL_HOME/bin/proddl-wrapper makeflow --help` for the list 
of possible Makeflow arguments, including all flavors of the batch 
systems that it can support.

For HPC MPI clusters, you can reuse jobs submission scripts
and Makeflow parameters generated by our pipeline 
[PGP](https://bitbucket.org/andreyto/proteogenomics).

For Microsoft Azure cloud, you can run PRODDL under our
client-cloud execution framework 
[PRODDL-C](https://bitbucket.org/andreyto/proddl-c).

You can also expose PRODDL as a Galaxy tool by deploying
into Galaxy the tool files from `api/galaxy/`.

Run `PRODDL_HOME/bin/proddl-wrapper proddl-dock dock --help`
to get the list of all available options for the docking
protocol.

