# AxisSPHYNX

**AxisSPHYNX** is a 2D axial-symmetric Smoothed Particle Hydrodynamics (SPH) code based on the 3D SPH code for astrophysics [SPHYNX](https://astro.physik.unibas.ch/sphynx). The main reference for AxisSPHYNX can be found [here](https://academic.oup.com/mnras/article/392/1/346/1073465).

**AxisSPHYNXMHD** extends AxisSPHYNX to the MHD realm by adding the magnetic-stress tensor to the axisymmetric SPH equations. Furthermore, the induction and dissipative equations are consistently written in such geometry.

In this repository we currently only include AxisSPHYNXMHD which is ready to simulate a Z-Pinch test out of the box.

# Parallelization

Currently, AxisSPHYNXMHD uses OpenMP.

# Compilation

Any fortran compiler with OpenMP support should suffice. Simply compile the source file, with the corresponding openmp option, and launch the executable.

Example of compilation with GNU Fortran:
```
gfortran -fopenmp AxisSPHYNXMHD.f90 -o AxisSPHYNXMHD
./AxisSPHYNXMHD
```
If you want to limit the number of OpenMP threads you can export the OMP_NUM_THREADS environment variable before executing. For example, for using 10 threads:
```
export OMP_NUM_THREADS=10
```
