High-resolution three-dimensional (3D) seismic simulation imposes severe demands on computational memory, 
and low-storage seismic simulation becomes particularly important.
Due to high-efficiency and low-storage, half-precision floating-point 16-bit format (FP16) 
is widely used in heterogeneous computing platforms,
such as Sunway series supercomputers and GPU computing platforms, etc..
Furthermore, Low-Storage Runge-Kutta (LSRK) requires lower memory resources compared with classical Runge-Kutta.
Therefore, FP16 and LSRK provide the possibility for low-storage seismic simulation.
However, the magnitude orders of physical quantities (velocity, stress and $Lam\acute{e}$ constants) 
in elastic wave equations are influenced by P-wave and S-wave velocities and densities of the elastic media.
This will result in a huge magnitude order difference between the stored values of velocity and stress, 
which have exceeds the range of stored values with FP16.
In this paper, we introduce 3 dimensionless constants Cv, Cs and Cp into elastic wave equations and
the new elastic wave equations are generated.
The 3 constants Cv, Cs and Cp can keep the magnitude orders of velocity and stress 
at a similar magnitude level in the new elastic wave equations. 
Then, the stored values of the variables in the new equations can maintain in the range of the stored values with FP16.
In addition, we introduce LSRK due to its low-storage property.
In our paper, based on low-storage technique of FP16 and LSRK, 
we develope an improved Multi-GPU solver for seismic simulation with 
curvilinear grid finite-difference method (CGFDM).
Moreover, we perform series of seismic simulations to verify the correctness and the validity of the improved solver coupled 
with the two techniques.
The verifications indicates that 
with maintaining the calculation accuracy, the computational efficiency of the solver is significantly improved, 
and the memory usage is remarkably reduced.
Expecially in the best condition, the memory usage can be reduced to nearly 1/3 of the original CGFDM solver.
