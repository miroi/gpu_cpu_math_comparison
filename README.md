gpu_cpu_math_comparison
=======================

Comparison of few GPU and CPU mathematical routines from BLAS and LAPACK libraries.

This simple project is focused at testing BLAS and LAPACK functions in GPU (cublas,cula)
and CPU (blas, lapack) realms.

Example:
--------
Comparing DGEMM, CPU_CUDA vs CPU_BLAS, matrix size 2000x2000, Tesla K20m 
```
 gpu-cpu-comparison.x  dgemm_cuda 2000
 .
 .
 .
  cuda_dgemm GPU vs CPU test
  Routine <GPU CUDA DGEMM> spent runtime:   0.05900000 seconds.
 aver. sum of diag:   2001.0000000000000      1,1:   2.0000000000000000      1,n:   2001.0000000000000      n,1:   2001.0000000000000      n,n   4000.0000000000000 
  Routine <CPU BLAS DGEMM> spent runtime:   7.13600000 seconds.
 aver. sum of diag:   2001.0000000000000      1,1:   2.0000000000000000      1,n:   2001.0000000000000      n,1:   2001.0000000000000      n,n   4000.0000000000000 
```
