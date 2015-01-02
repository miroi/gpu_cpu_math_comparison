/*

 Space for C routines of GPGPU libraries CUDA & CULA.

 Written by Miro Ilias & Daniel Kuzma in order to extend Fortran language handling
 of GPU mathematical libraries.

*/

#include "stdio.h"
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>

// we need CUDA Toolkit routines...
#include "cublas.h"
#include "cuda.h"
#include "cuda_runtime_api.h"

/* cula */
#if defined USE_CULA
#include <cula.h>
#include <cula_blas.h>
//#include <cula_pack.h>
#endif

/* cblas */
#if defined HAVE_MKL_BLAS
/* see /opt/intel/mkl/include/mkl_cblas.h */
#include "mkl_cblas.h"
#pragma message "Using Intel MKL <mkl_cblas.h> interface"
#else
/* see /usr/include/cblas.h */
#include "cblas.h"
#pragma message "Using GNU <cblas.h> interface"
#endif

/* lapack */
#if defined HAVE_MKL_LAPACK
/* see /opt/intel/mkl/include/mkl_lapack(e).h  */
//#include <mkl_lapack.h>
#include <mkl_lapacke.h>
#pragma message "Using Intel MKL <mkl_lapacke.h> interface"
#else
#include <clapack.h>
#pragma message "Using GNU <clapack.h> interface"
#endif


int NRHS = 1; //for cula sgesv


/* C routine is necessary to obtain info about CUDA library (not possible from Fortran routine) */
void cudaabouts_c_(void)
{
 int cudaDriverVersion, cudaRuntimeVersion;
 cudaError_t stat;
 stat=cudaDriverGetVersion(&cudaDriverVersion); // get CUDA driver version
 stat=cudaRuntimeGetVersion(&cudaRuntimeVersion); // get runtimeversion
 printf("--- CUDA GPU info (in C routine) ---\n");
 printf(" CUDA driver version  :%d \n",cudaDriverVersion );
 printf(" CUDA runtime version :%d \n",cudaRuntimeVersion );
 return;
}

/* flush out some info about C compiler */
void c_compiler_info_(void)
{
// http://gcc.gnu.org/onlinedocs/cpp/Common-Predefined-Macros.html
  printf("--- C compiler info (in C routine) ---\n");
#if __GNUC__
  printf(" GNU gcc version: %s \n",__VERSION__);
// printf("file: %s\n",__BASE_FILE__);
#endif

#if __ECC
  printf(" Intel version %s",__ECC);
#endif

#if __ICC
  printf (" Intel __ICC %i \n",__ICC);
#endif

#if __INTEL_COMPILER
  printf(" Intel compiler macro %i\n",__INTEL_COMPILER);
#endif

#if __INTEL_COMPILER_BUILD_DATE
 // printf("INTEL_COMPILER_BUILD_DATE:%s \n",__INTEL_COMPILER_BUILD_DATE);
#endif

}


/* remove empty space from routine name */
char *trim_(char *buf)
{
  char *i=(char*)buf, *j=(char*)buf;  
  do {  
    if (*i != ' ')  
      *(j++) = *i;  
  } while (*(i++));  
  return buf;
}

/*   
 Space for testing GPU routines written in C
*/

double getTime()
{
 struct timeval tx;
 double t;
 gettimeofday(&tx, NULL);
 t = (tx.tv_sec) * 1000.0;
 t += (tx.tv_usec) / 1000.0;
 return t;
}


void sgemm_cuda_c_test(int *n, int *verbosity){
        
        cublasStatus status;
        float* h_A;
        float* h_B;
        float* h_C;
        float* h_C_ref;
        float* d_A = 0;
        float* d_B = 0;
        float* d_C = 0;
        float alpha = 1.0f;
        float beta = 0.0f;
        int i;
        float error_norm;
        float ref_norm;
        float diff;
        double gpu_time_1,gpu_time_2,gpu_time; // zaciatocny a konecny cas na gpu
        double cpu_time_1,cpu_time_2; // zaciatocny a konecny cas na cpu
        int n2,n3;
        
        n3 = *n;
        n2 = n3*n3;
        //printf("\n\n\n N2:%i N:%i N*:%i\n\n\n",n2,n3,*n);
        /* Initialize CUBLAS */
        //printf("\nCUBLAS SGEMM test running...\n\n");

        //printf("\n\n");
        printf("-------------------\n");
        printf("   GPU CUBLAS_SGEMM  \n");
        printf("-------------------\n");

        status = cublasInit();
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! CUBLAS initialization error\n");
                exit(EXIT_FAILURE);
        } else { fprintf(stdout,"GPU cublas initialized.\n"); }


        /* Allocate host memory for the matrices */
        h_A = (float*)malloc(n2 * sizeof(h_A[0]));
        if (h_A == 0) {
                fprintf (stderr, "!!!! host memory allocation error (A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A mtx allocated \n"); }


        h_B = (float*)malloc(n2 * sizeof(h_B[0]));
        if (h_B == 0) {
                fprintf (stderr, "!!!! host memory allocation error (B)\n");
                exit(EXIT_FAILURE);
        }  else {if(*verbosity>0)fprintf(stdout,"B mtx allocated \n"); }

        h_C = (float*)malloc(n2 * sizeof(h_C[0]));
        if (h_C == 0) {
                fprintf (stderr, "!!!! host memory allocation error (C)\n");
                exit(EXIT_FAILURE);
        }  else {if(*verbosity>0)fprintf(stdout,"C mtx allocated \n"); }

        /* Fill the matrices with test data */
        for (i = 0; i < n2; i++) {
                h_A[i] = rand() / (float)RAND_MAX;
                h_B[i] = rand() / (float)RAND_MAX;
                h_C[i] = 2;//rand() / (float)RAND_MAX;
        }
        if(*verbosity>0)fprintf(stdout,"A,B,C matrixes filled \n");

        /* Allocate device memory for the matrices A,B and C */
        status = cublasAlloc(n2, sizeof(d_A[0]), (void**)&d_A);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device memory allocation error (A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A mtx allocated on GPU \n"); }

        status = cublasAlloc(n2, sizeof(d_B[0]), (void**)&d_B);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device memory allocation error (B)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"B mtx allocated on GPU \n"); }

        status = cublasAlloc(n2, sizeof(d_C[0]), (void**)&d_C);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device memory allocation error (C)\n");
                exit(EXIT_FAILURE);
        }  else {if(*verbosity>0)fprintf(stdout,"C mtx allocated on GPU \n"); }
        /* Initialize the device matrices with the host matrices */
        gpu_time_1 = getTime();
        status = cublasSetVector(n2, sizeof(h_A[0]), h_A, 1, d_A, 1);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1;
        if(*verbosity>=2)printf("Filling time of A mtx on GPU: %f sekund\n",gpu_time/1000);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device access error (write A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A mtx filled on GPU \n"); }
        gpu_time_1 = getTime();
        status = cublasSetVector(n2, sizeof(h_B[0]), h_B, 1, d_B, 1);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1 + gpu_time;
        if(*verbosity>=2)printf("Filling time of B mtx on GPU: %f sekund\n",(gpu_time_2 - gpu_time_1)/1000);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device access error (write B)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"B mtx filled on GPU \n"); }
        gpu_time_1 = getTime();
        status = cublasSetVector(n2, sizeof(h_C[0]), h_C, 1, d_C, 1);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1 + gpu_time;
        if(*verbosity>=2)printf("Filling time of C mtx on GPU: %f sekund\n",(gpu_time_2 - gpu_time_1)/1000);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device access error (write C)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"C mtx filled on GPU \n"); }
        
        
        /* Clear last error */
        cublasGetError();

        /* Performs operation using cublas */
       // matrix-matrix multiplication  na GPU ...
        printf("GPU cublasSgemm running....");
        fflush(stdout);
        gpu_time_1 = getTime();
        cublasSgemm('n', 'n', n3, n3, n3, alpha, d_A, n3, d_B, n3, beta, d_C, n3);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1 + gpu_time;
        
        status = cublasGetError();
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! kernel execution error.\n");
                exit(EXIT_FAILURE);
        } else  {fprintf(stdout,"...done.\n"); }
        if(*verbosity>=2)printf("Time of computation on GPU: %f sekund\n",(gpu_time_2 - gpu_time_1)/1000);
        /* Allocate host memory for reading back the result from the device memory */
        h_C = (float*)malloc(n2 * sizeof(h_C[0]));
        if (h_C == 0) {
                fprintf (stderr, "!!!! host memory allocation error (C)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"New C mtx ready for GPU data. \n");}

        /* Read the result back */
        gpu_time_1 = getTime();
        status = cublasGetVector(n2, sizeof(h_C[0]), d_C, 1, h_C, 1);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1 + gpu_time;
        if(*verbosity>=2)printf("Obtaining time from GPU: %f sekund\n",(gpu_time_2 - gpu_time_1)/1000);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device access error (read C)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"C mtx obtained from GPU device. \n");}
        //printf("...done\n");
        
        h_C_ref = h_C;

#if defined CBLAS
        //printf("\n\n");
        printf("-------------------\n");
        printf("   CPU CBLAS_SGEMM  \n");
        printf("-------------------\n");
        
        printf("CPU cblas_sgemm running....");
        fflush(stdout);
    
        cpu_time_1 = getTime();

/*
 void cblas_sgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const float alpha, const float *A,
                 const int lda, const float *B, const int ldb,
                 const float beta, float *C, const int ldc);
*/
        cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n3, n3,n3, alpha, h_A, n3, h_B, n3, beta, h_C, n3);
        cpu_time_2 = getTime();
        printf("...done\n");
        //printf("\n\n");
#endif
        
           printf("-------------------\n");
           printf("   CHECK RESULT SGEMM  \n");
           printf("-------------------\n");
           float pole = 0;
           float newes;
           //printf("\n");
           for (i=0;i<n3;i++)
            {
              newes = fabs(h_C_ref[i]-h_C[i]);
              pole = pole + newes;
            if(*verbosity>2)  printf("Result check :  i=%d  GPU[%d]=%f  CPU[%d]=%f  Difference :%f \n",i, i+1,h_C[i],i+1, h_C_ref[i],newes); 
     
      
            }
       if(*verbosity==2)printf("Results : CPU - First element of matrix - %f / GPU - First element of matrix - %f\n",h_C[0],h_C_ref[0]);
        if(*verbosity<=2)printf("Results : CPU - Last element of matrix - %f / GPU - Last element of matrix - %f\n",h_C[n3-1],h_C_ref[n3-1]);
       printf("Average absolute difference :%f \n", pole/n3);
        
        printf("-------------------\n");
        printf("   TIMES SGEMM  \n");
        printf("-------------------\n");
        printf("Time of computation on CPU: %f sekund\n",(cpu_time_2-cpu_time_1)/1000);
        printf("Overall time of computation and transfer on GPU: %f sekund\n\n",gpu_time/1000);
        
        /* Memory clean up */
        free(h_A);
        free(h_B);
        free(h_C);
        //free(h_C_ref);
        status = cublasFree(d_A);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! memory free error (A)\n");
                exit(EXIT_FAILURE);
        }
        status = cublasFree(d_B);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! memory free error (B)\n");
                exit(EXIT_FAILURE);
        }
        status = cublasFree(d_C);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! memory free error (C)\n");
                exit(EXIT_FAILURE);
        }

        /* Shutdown */
        status = cublasShutdown();
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! shutdown error (A)\n");
                exit(EXIT_FAILURE);
        }

}

void sgemv_cuda_c_test(int *n, int *verbosity){

cublasStatus status;
        float* h_A;
        float* h_B;
        float* h_C;
        float* h_C_ref;
        float* d_A = 0;
        float* d_B = 0;
        float* d_C = 0;
        float alpha = 1.0f;
        float beta = 0.0f;
        int i,j;
        float error_norm;
        float ref_norm;
        float diff;
        double gpu_time_1,gpu_time_2,gpu_time; // zaciatocny a konecny cas na gpu
        double cpu_time_1,cpu_time_2; // zaciatocny a konecny cas na cpu

        int n2,n3;
        
        n3 = *n;
        n2 = n3*n3;
        //printf("\n\n\n N2:%i N:%i N*:%i\n\n\n",n2,n3,*n);
        /* Initialize CUBLAS */
        //printf("\n\n");
        printf("-------------------\n");
        printf("   GPU CUBLAS_SGEMV  \n");
        printf("-------------------\n");

        status = cublasInit();
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! CUBLAS initialization error\n");
                exit(EXIT_FAILURE);
        } else { fprintf(stdout,"GPU cublas initialized.\n"); }


        /* Allocate host memory for the matrices */
        h_A = (float*)malloc(n2 * sizeof(h_A[0]));
        if (h_A == 0) {
                fprintf (stderr, "!!!! host memory allocation error (A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A mtx allocated \n"); }


        h_B = (float*)malloc(n3 * sizeof(h_B[0]));
        if (h_B == 0) {
                fprintf (stderr, "!!!! host memory allocation error (B)\n");
                exit(EXIT_FAILURE);
        }  else {if(*verbosity>0)fprintf(stdout,"B vector allocated \n"); }

        h_C = (float*)malloc(n3 * sizeof(h_C[0]));
        if (h_C == 0) {
                fprintf (stderr, "!!!! host memory allocation error (C)\n");
                exit(EXIT_FAILURE);
        }  else {if(*verbosity>0)fprintf(stdout,"C vector allocated \n"); }

        /* Fill the matrices with test data */
        for (i = 0; i < n2; i++) {
                h_A[i] =  rand() / (float)RAND_MAX;
                
        }
        for (j = 0; j < n3; j++) {
                
                h_B[j] = rand() / (float)RAND_MAX;
                h_C[j] = 2;//rand() / (float)RAND_MAX;
        }    
        if(*verbosity>0)fprintf(stdout,"A matrixes filled and B,C vector filled\n");

        /* Allocate device memory for the matrices A and vector B and C */
        status = cublasAlloc(n2, sizeof(d_A[0]), (void**)&d_A);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device memory allocation error (A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A mtx allocated on GPU \n"); }

        status = cublasAlloc(n3, sizeof(d_B[0]), (void**)&d_B);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device memory allocation error (B)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"B vector allocated on GPU \n"); }

        status = cublasAlloc(n3, sizeof(d_C[0]), (void**)&d_C);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device memory allocation error (C)\n");
                exit(EXIT_FAILURE);
        }  else {if(*verbosity>0)fprintf(stdout,"C vector allocated on GPU \n"); }
        /* Initialize the device matrices with the host matrices */
        gpu_time_1 = getTime();
        status = cublasSetVector(n2, sizeof(h_A[0]), h_A, 1, d_A, 1);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1;
        if(*verbosity>=2)printf("Filling time of A mtx on GPU: %f sekund\n",gpu_time/1000);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device access error (write A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A mtx filled on GPU \n"); }
        gpu_time_1 = getTime();
        status = cublasSetVector(n3, sizeof(h_B[0]), h_B, 1, d_B, 1);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1 + gpu_time;
        if(*verbosity>=2)printf("Filling time of B vector on GPU: %f sekund\n",(gpu_time_2 - gpu_time_1)/1000);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device access error (write B)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"B vector filled on GPU \n"); }
        gpu_time_1 = getTime();
        status = cublasSetVector(n3, sizeof(h_C[0]), h_C, 1, d_C, 1);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1 + gpu_time;
        if(*verbosity>=2)printf("Filling time of C vector on GPU: %f sekund\n",(gpu_time_2 - gpu_time_1)/1000);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device access error (write C)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"C vector filled on GPU \n"); }
        
         /* Clear last error */
         cublasGetError();

        /* Performs operation using cublas */
       // matrix-matrix multiplication  na GPU ...
        printf("GPU cublasSgemv running....");
        fflush(stdout);
        gpu_time_1 = getTime();
        cublasSgemv('t', n3, n3, alpha, d_A, n3, d_B, 1, beta, d_C, 1);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1 + gpu_time;
        
        status = cublasGetError();
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! kernel execution error.\n");
                exit(EXIT_FAILURE);
        } else  {fprintf(stdout,"...done.\n"); }
        if(*verbosity==2)printf("Time of computation on GPU: %f sekund\n",(gpu_time_2 - gpu_time_1)/1000);
        /* Allocate host memory for reading back the result from the device memory */
        h_C = (float*)malloc(n2 * sizeof(h_C[0]));
        if (h_C == 0) {
                fprintf (stderr, "!!!! host memory allocation error (C)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"New C mtx ready for GPU data. \n");}

        /* Read the result back */
        gpu_time_1 = getTime();
        status = cublasGetVector(n3, sizeof(h_C[0]), d_C, 1, h_C, 1);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1 + gpu_time;
        if(*verbosity>=2)printf("Obtaining time from GPU: %f sekund\n",(gpu_time_2 - gpu_time_1)/1000);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device access error (read C)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"C mtx obtained from GPU device. \n");}
        
        h_C_ref = h_C;


#if defined USE_CBLAS
        //printf("\n\n");
        printf("-------------------\n");
        printf("   CPU CBLAS_SGEMV  \n");
        printf("-------------------\n");
        
        printf("CPU cblas_sgemv running....");
        fflush(stdout);
    
        cpu_time_1 = getTime();

/*
void cblas_sgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 const float *X, const int incX, const float beta,
                 float *Y, const int incY);
*/
        cblas_sgemv(CblasRowMajor, CblasNoTrans, n3,n3, alpha, h_A, n3, h_B, 1, beta, h_C, 1);
        cpu_time_2 = getTime();

        printf("...done\n");
#endif        
        //printf("\n\n");
           
           printf("-------------------\n");
           printf("   CHECK RESULT SGEMV  \n");
           printf("-------------------\n");
           float pole = 0;
           float newes;
           //printf("\n");
           for (i=0;i<n3;i++)
            {
              newes = fabs(h_C_ref[i]-h_C[i]);
              pole = pole + newes;
            if(*verbosity>2)  printf("Result check :  i=%d  GPU[%d]=%f  CPU[%d]=%f  Difference :%f \n",i, i+1,h_C[i],i+1, h_C_ref[i],newes); 
     
      
            }
           if(*verbosity==2)printf("Results : CPU - First element of matrix - %f / GPU - First element of matrix - %f\n",h_C[0],h_C_ref[0]);
          if(*verbosity<=2)printf("Results : CPU - Last element of matrix - %f / GPU - Last element of matrix - %f\n",h_C[n3-1],h_C_ref[n3-1]);
          printf("Average absolute difference :%f \n", pole/n3); 

        
       
        //printf("...done\n");
        printf("-------------------\n");
        printf("   TIMES SGEMV  \n");
        printf("-------------------\n");
                
        printf("Time of computation on CPU: %f sekund\n",(cpu_time_2-cpu_time_1)/1000);
        printf("Overall time of computation and transfer on GPU: %f sekund\n",gpu_time/1000);
        
        /* Memory clean up */
        free(h_A);
        free(h_B);
        free(h_C);
        //free(h_C_ref);
        status = cublasFree(d_A);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! memory free error (A)\n");
                exit(EXIT_FAILURE);
        }
        status = cublasFree(d_B);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! memory free error (B)\n");
                exit(EXIT_FAILURE);
        }
        status = cublasFree(d_C);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! memory free error (C)\n");
                exit(EXIT_FAILURE);
        }

        /* Shutdown */
        status = cublasShutdown();
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! shutdown error (A)\n");
                exit(EXIT_FAILURE);
        }
}

#if defined USE_CULA

void sgemm_cula_c_test(int *n, int *verbosity){
        
        culaStatus status;
        culaFloat* d_A = 0;
        culaFloat* d_B = 0;
        culaFloat* d_C = 0;
        culaInt* IPIV = NULL;
        float* h_C_ref;
        float* h_A = 0;
        float* h_B = 0;
        float* h_C = 0;
        float alpha = 1.0f;
        float beta = 0.0f;
        int i;
        float error_norm;
        float ref_norm;
        float diff;
        double gpu_time_1,gpu_time_2,gpu_time; // zaciatocny a konecny cas na gpu
        double cpu_time_1,cpu_time_2; // zaciatocny a konecny cas na cpu

        int n2,n3;
        
        n3 = *n;
        n2 = n3*n3;
        //printf("\n\n\n N2:%i N:%i N*:%i\n\n\n",n2,n3,*n);
        /* Initialize CUBLAS */
        //printf("\n\n");
        printf("-------------------\n");
        printf("   GPU CULA_SGEMM  \n");
        printf("-------------------\n");
        //printf("\nCULA SGEMM test running...\n\n");

        status = culaInitialize();
        if (status != culaNoError) {
                fprintf (stderr, "!!!! CULA initialization error\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"GPU cula initialized.\n"); }


        /* Allocate host memory for the matrices */
        h_A = (float*)malloc(n2 * sizeof(h_A[0]));
        if (h_A == 0) {
                fprintf (stderr, "!!!! host memory allocation error (A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A mtx allocated \n"); }


        h_B = (float*)malloc(n2 * sizeof(h_B[0]));
        if (h_B == 0) {
                fprintf (stderr, "!!!! host memory allocation error (B)\n");
                exit(EXIT_FAILURE);
        }  else {if(*verbosity>0)fprintf(stdout,"B mtx allocated \n"); }

        h_C = (float*)malloc(n2 * sizeof(h_C[0]));
        if (h_C == 0) {
                fprintf (stderr, "!!!! host memory allocation error (C)\n");
                exit(EXIT_FAILURE);
        }  else {if(*verbosity>0)fprintf(stdout,"C mtx allocated \n"); }

        /* Allocate device memory for the matrices A,B and C */
        d_A = (culaFloat*)malloc(n2*sizeof(culaFloat));
        if (d_A == 0) {
                fprintf (stderr, "!!!! device memory allocation error (A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A mtx allocated on graphic card\n"); }
        
        d_B = (culaFloat*)malloc(n2*sizeof(culaFloat));
        if (d_B == 0) {
                fprintf (stderr, "!!!! device memory allocation error (B)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"B mtx allocated on graphic card\n"); }
        d_C = (culaFloat*)malloc(n2*sizeof(culaFloat));
        if (d_C == 0) {
                fprintf (stderr, "!!!! device memory allocation error (C)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"C mtx allocated on graphic card\n"); }


        /* Fill the matrices with test data */
        for (i = 0; i < n2; i++) {
               d_A[i] = h_A[i] = rand() / (float)RAND_MAX;
                //d_A[i] = rand() / (float)RAND_MAX;
               d_B[i] = h_B[i] = rand() / (float)RAND_MAX;
                //d_B[i] = rand() / (float)RAND_MAX;
                h_C[i] = 2;//rand() / (float)RAND_MAX;
                d_C[i] = 21;//rand() / (float)RAND_MAX;
        }
        if(*verbosity>0)fprintf(stdout,"A,B,C matrixes filled \n");

        
        /* Clear last error */
        cublasGetError();

        /* Performs operation using cublas */
       // matrix-matrix multiplication  na GPU ...
        printf("GPU culaSgemm running....");
        fflush(stdout);
        gpu_time_1 = getTime();
        status = culaSgemm('N', 'N', n3, n3, n3, alpha, d_A, n3, d_B, n3, beta, d_C, n3);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1;
        //checkStatus(status);
        //status = cublasGetError();
        if (status != culaNoError) {
                fprintf (stderr, "!!!! kernel execution error.\n");
                exit(EXIT_FAILURE);
        } else  {fprintf(stdout,"...done.\n"); }
        
        h_C_ref = h_C;
       
      // printf("\n\n");
        printf("-------------------\n");
        printf("   CPU CBLAS_SGEMM  \n");
        printf("-------------------\n");
        printf("CPU cblas_sgemm running....");
        fflush(stdout);
    
        cpu_time_1 = getTime();
        cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n3, n3,n3, alpha, h_A, n3, h_B, n3, beta, h_C, n3);
        cpu_time_2 = getTime();

        printf("...done\n");
        
        // printf("\n\n");
        
           printf("-------------------\n");
           printf("   CHECK RESULT SGEMM  \n");
           printf("-------------------\n");
           float pole = 0;
           float newes;
           //printf("\n");
           for (i=0;i<n3;i++)
            {
              newes = fabs(h_C_ref[i]-h_C[i]);
              pole = pole + newes;
            if(*verbosity>2)  printf("Result check :  i=%d  GPU[%d]=%f  CPU[%d]=%f  Difference :%f \n",i, i+1,h_C[i],i+1, h_C_ref[i],newes); 
     
      
            }
       if(*verbosity==2)printf("Results : CPU - First element of matrix - %f / GPU - First element of matrix - %f\n",h_C[0],h_C_ref[0]);
        if(*verbosity<=2)printf("Results : CPU - Last element of matrix - %f / GPU - Last element of matrix - %f\n",h_C[n3-1],h_C_ref[n3-1]);
       printf("Average absolute difference :%f \n", pole/n3);
        
        printf("-------------------\n");
        printf("   TIMES SGEMM  \n");
        printf("-------------------\n");
        printf("Time of computation on CPU: %f sekund\n",(cpu_time_2-cpu_time_1)/1000);
        printf("Overall time of computation and transfer on GPU: %f sekund\n",gpu_time/1000);

        
        //if(*verbosity==2)printf("Time of computation on GPU: %f sekund\n",(gpu_time_2 - gpu_time_1)/1000);
        /* Allocate host memory for reading back the result from the device memory */
        
        //printf("\n\n\n");
        //if(*verbosity>0)printf("Results : CPU - First element of matrix - %f / GPU - First element of matrix - %f\n",h_C_ref[0],d_C[0]);
        //printf("Results : CPU - Last element of matrix - %f / GPU - Last element of matrix - %f\n\n\n",h_C_ref[n2-1],d_C[n2-1]);
        
        //printf("Time of computation on CPU: %f sekund\n",(cpu_time_2-cpu_time_1)/1000);
        //printf("Overall time of computation and transfer on GPU: %f sekund\n\n\n",gpu_time/1000);
        
        /* Memory clean up */
        free(d_A);
        free(d_B);
        free(d_C);
        free(h_A);
        free(h_B);
        free(h_C);
        //free(h_C_ref);

        /* Shutdown */
        culaShutdown();
}

void sgemv_cula_c_test(int *n, int *verbosity){
        
        culaStatus status;
        culaFloat* d_A = 0;
        culaFloat* d_B = 0;
        culaFloat* d_C = 0;
        culaInt* IPIV = NULL;
        float* h_C_ref;
        float* h_A = 0;
        float* h_B = 0;
        float* h_C = 0;
        float alpha = 1.0f;
        float beta = 0.0f;
        int i;
        float error_norm;
        float ref_norm;
        float diff;
        double gpu_time_1,gpu_time_2,gpu_time; // zaciatocny a konecny cas na gpu
        double cpu_time_1,cpu_time_2; // zaciatocny a konecny cas na cpu

        int n2,n3;
        
        n3 = *n;
        n2 = n3*n3;
        //printf("\n\n\n N2:%i N:%i N*:%i\n\n\n",n2,n3,*n);
        /* Initialize CUBLAS */
        //printf("\n\n");
        printf("-------------------\n");
        printf("   GPU CULA_SGEMV  \n");
        printf("-------------------\n");
        //printf("\nCULA SGEMV test running...\n\n");

        status = culaInitialize();
        if (status != culaNoError) {
                fprintf (stderr, "!!!! CULA initialization error\n");
                exit(EXIT_FAILURE);
        } else { if(*verbosity>0) fprintf(stdout,"GPU cula initialized.\n"); }


        /* Allocate host memory for the matrices */
        h_A = (float*)malloc(n2 * sizeof(h_A[0]));
        if (h_A == 0) {
                fprintf (stderr, "!!!! host memory allocation error (A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A mtx allocated \n"); }


        h_B = (float*)malloc(n3 * sizeof(h_B[0]));
        if (h_B == 0) {
                fprintf (stderr, "!!!! host memory allocation error (B)\n");
                exit(EXIT_FAILURE);
        }  else {if(*verbosity>0)fprintf(stdout,"B mtx allocated \n"); }

        h_C = (float*)malloc(n3 * sizeof(h_C[0]));
        if (h_C == 0) {
                fprintf (stderr, "!!!! host memory allocation error (C)\n");
                exit(EXIT_FAILURE);
        }  else {if(*verbosity>0)fprintf(stdout,"C mtx allocated \n"); }

        /* Allocate device memory for the matrices A,B and C */
        d_A = (culaFloat*)malloc(n2*sizeof(culaFloat));
        if (d_A == 0) {
                fprintf (stderr, "!!!! device memory allocation error (A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A mtx allocated on graphic card\n"); }
        
        d_B = (culaFloat*)malloc(n3*sizeof(culaFloat));
        if (d_B == 0) {
                fprintf (stderr, "!!!! device memory allocation error (B)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"B mtx allocated on graphic card\n"); }
        d_C = (culaFloat*)malloc(n3*sizeof(culaFloat));
        if (d_C == 0) {
                fprintf (stderr, "!!!! device memory allocation error (C)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"C mtx allocated on graphic card\n"); }


        /* Fill the matrices with test data */
        for (i = 0; i < n2; i++) {
               d_A[i] = h_A[i] = rand() / (float)RAND_MAX;
        }
        for (i = 0; i < n3; i++) {
               d_B[i] = h_B[i] = rand() / (float)RAND_MAX;
                h_C[i] = 2;//rand() / (float)RAND_MAX;
                d_C[i] = 21;//rand() / (float)RAND_MAX;
        }
        if(*verbosity>0)fprintf(stdout,"A,B,C matrixes filled \n");
       
        /* Clear last error */
        cublasGetError();

        /* Performs operation using cublas */
       // matrix-matrix multiplication  na GPU ...
        printf("GPU culaSgemv running....");
        fflush(stdout);
        gpu_time_1 = getTime();
        status = culaSgemv('t',n3,n3,alpha,d_A,n3,d_B,1,beta,d_C,1);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1;
        //checkStatus(status);
        //status = cublasGetError();
        if (status != culaNoError) {
                fprintf (stderr, "!!!! kernel execution error.\n");
                exit(EXIT_FAILURE);
        } else  {fprintf(stdout,"...done.\n"); }
       
        h_C_ref = h_C;
        //printf("\n\n");
        printf("-------------------\n");
        printf("   CPU CBLAS_SGEMV  \n");
        printf("-------------------\n");
        printf("CPU cblas_sgemv running....");
        fflush(stdout);
    
        cpu_time_1 = getTime();
        cblas_sgemv(CblasRowMajor, CblasNoTrans, n3,n3, alpha, h_A, n3, h_B, 1, beta, h_C, 1);
        cpu_time_2 = getTime();

        printf("...done\n");
        
        //printf("\n\n");
        
           printf("-------------------\n");
           printf("   CHECK RESULT SGEMV  \n");
           printf("-------------------\n");
           float pole = 0;
           float newes;
           //printf("\n");
           for (i=0;i<n3;i++)
            {
              newes = fabs(h_C_ref[i]-h_C[i]);
              pole = pole + newes;
            if(*verbosity>2)  printf("Result check :  i=%d  GPU[%d]=%f  CPU[%d]=%f  Difference :%f \n",i, i+1,h_C[i],i+1, h_C_ref[i],newes); 
     
      
            }
       if(*verbosity==2)printf("Results : CPU - First element of matrix - %f / GPU - First element of matrix - %f\n",h_C[0],h_C_ref[0]);
        if(*verbosity<=2)printf("Results : CPU - Last element of matrix - %f / GPU - Last element of matrix - %f\n",h_C[n3-1],h_C_ref[n3-1]);
       printf("Average absolute difference :%f \n", pole/n3);
        
        printf("-------------------\n");
        printf("   TIMES SGEMV  \n");
        printf("-------------------\n");
        printf("Time of computation on CPU: %f sekund\n",(cpu_time_2-cpu_time_1)/1000);
        printf("Overall time of computation and transfer on GPU: %f sekund\n",gpu_time/1000);
       
        //if(*verbosity==2)printf("Time of computation on GPU: %f sekund\n",(gpu_time_2 - gpu_time_1)/1000);
        /* Allocate host memory for reading back the result from the device memory */
        
        //printf("\n\n\n");
        //if(*verbosity>0)printf("Results : CPU - First element of matrix - %f / GPU - First element of matrix - %f\n",h_C[0],d_C[0]);
        //printf("Results : CPU - Last element of matrix - %f / GPU - Last element of matrix - %f\n\n\n",h_C[n3-1],d_C[n3-1]);
        
        //printf("Time of computation on CPU: %f sekund\n",(cpu_time_2-cpu_time_1)/1000);
        //printf("Overall time of computation and transfer on GPU: %f sekund\n\n\n",gpu_time/1000);
        
        /* Memory clean up */
        free(d_A);
        free(d_B);
        free(d_C);
        free(h_A);
        free(h_B);
        free(h_C);
        //free(h_C_ref);  

        /* Shutdown */
        culaShutdown();
}

void sgesv_cula_c_test(int *n, int *verbosity){
    int i, info,j;
    int n3 = *n;

    double gpu_time_1,gpu_time_2,gpu_time; // zaciatocny a konecny cas na gpu
    double cpu_time_1,cpu_time_2; // zaciatocny a konecny cas na cpu
    //n3=5;
    culaStatus status;
    float* A1 = NULL;
    culaFloat* A = NULL;
    culaFloat* B = NULL;
    culaFloat* X = NULL;
    culaInt* IPIV = NULL;

    float h[25] = {
            6.80f, -2.11f,  5.66f,  5.97f,  8.23f,
           -6.05f, -3.30f,  5.36f, -4.44f,  1.08f,
           -0.45f,  2.58f, -2.70f,  0.27f,  9.04f,
            8.32f,  2.71f,  4.35f, -7.17f,  2.14f,
           -9.67f, -5.14f, -7.26f,  6.08f, -6.87f
        };
    float u[15] = {
            4.02f,  6.19f, -8.22f, -7.57f, -3.03f,
           -1.56f,  4.00f, -8.67f,  1.75f,  2.86f,
            9.81f, -4.09f, -4.57f, -8.61f,  8.99f
        };

              
    culaFloat one = 1.0f;
    culaFloat thresh = 1e-5f;//0.01f;//1e-6f;
    culaFloat diff;
    srand ( time(NULL) );
    
    
    
    void checkStatus(culaStatus status)
{
    char buf[80];
    if(!status) return;
    culaGetErrorInfoString(status, culaGetErrorInfo(), buf, sizeof(buf));
    printf("%s\n", buf);
    culaShutdown();
    exit(EXIT_FAILURE);
}

    void verify_result(void)
    {
    float pole = 0;
    
    printf("Verifying Result\n\n");
    for(i = 0; i < n3; ++i)
    {
        diff = X[i] - B[i];
        if(diff < 0.0f)
            diff = -diff;
        if(diff > thresh)
            if(*verbosity>2) printf("Result check failed:  i=%d  GPU[i]=%f  CPU[i]=%f\n", i, X[i], B[i]);
    }
    float newes;
    printf("\n\n");
    for (i=0;i<n3;i++)
        {
           newes = fabs(X[i]-B[i]);
           pole = pole + newes;
           if(*verbosity>2) printf("Result check :  i=%d  GPU[%d]=%f  CPU[%d]=%f  Difference :%f \n",i, i+1,X[i],i+1, B[i],newes); 
        }
    
    
    printf("\n\n Average absolute difference :%f \n\n", pole/n3);
    } 

    printf("-------------------\n");
    printf("   GPU CULA_SGESV  \n");
    printf("-------------------\n");

    A = (culaFloat*)malloc(n3*n3*sizeof(culaFloat));
        if (A == NULL) {fprintf (stderr, "!!!! device memory allocation error (A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A mtx allocated on graphic card\n"); }

    A1 = (float*)malloc(n3*n3 * sizeof(A1[0]));
    if (A1 == 0) {
                fprintf (stderr, "!!!! host memory allocation error (A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A1 mtx allocated \n"); }
    

    if(*verbosity>0)printf("Allocating Matrix A, vectors B, X and IPIV, size N=%d \n",n3);
    A = (culaFloat*)malloc(n3*n3*sizeof(culaFloat));
    B = (culaFloat*)malloc(n3*NRHS*sizeof(culaFloat));
    X = (culaFloat*)malloc(n3*NRHS*sizeof(culaFloat));
    IPIV = (culaInt*)malloc(n3*sizeof(culaInt));
    if(!A || !B || !IPIV)
        exit(EXIT_FAILURE);

    
    for (i=0;i<n3;i++)
        {
          B[i] = rand() / (float)RAND_MAX;//rand() / (float)RAND_MAX;//i+2;//rand() / (float)RAND_MAX;
          X[i] = B[i];//rand() / (float)RAND_MAX;
        }

    for (i=0;i<n3;i++) {
    for (j=0;j<n3;j++)
        {
          A1[j*n3+i] = rand() / (float)RAND_MAX;//rand() / (float)RAND_MAX;//j+1+i+1;//rand() / (float)RAND_MAX;
          A[j*n3+i] = A1[j*n3+i];
        }
       }
    
    if(*verbosity>0)fprintf(stdout,"A,A1 matrixes and X,B vectors filled \n");
    
    if (*verbosity>10){
    int poc = 0;
    int ooo = 0;
    for (i=0;i<n3;i++) {
    for (j=0;j<n3;j++)
        {
        printf(" %f ",A1[j*n3+i]);//  A1[j*n3+i] = i+1+(j+1)*(j+1);//rand() / (float)RAND_MAX;//j+1+i+1;//rand() / (float)RAND_MAX;
         // A[j*n3+i] = A1[j*n3+i];
        if ((ooo+1)%n3 == 0) {printf("| %f\n",B[poc]);poc++;}
        ooo++;
        }
       }
    
    /*for (i=0;i<n3*n3;i++){
      printf(" %f ",A1[i]);
      if ((i+1)%n3 == 0) {printf("| %f\n",B[poc]);poc++;}
    }          */
    
    
   // printf("\n\n");
    ooo=0;
    poc=0;
    
    
    for (i=0;i<n3;i++) {
    for (j=0;j<n3;j++)
        {
        printf(" %f ",A[j*n3+i]);//  A1[j*n3+i] = i+1+(j+1)*(j+1);//rand() / (float)RAND_MAX;//j+1+i+1;//rand() / (float)RAND_MAX;
         // A[j*n3+i] = A1[j*n3+i];
        if ((ooo+1)%n3 == 0) {printf("| %f\n",X[poc]);poc++;}
        ooo++;
        }
       }
    
    /*for (i=0;i<n3*n3;i++){
      printf(" %f ",A[i]);
      if ((i+1)%n3 == 0) {printf("| %f\n",B[poc]);poc++;}
    } */
    
        }
      
    
    
        
    /* *A = {
            6.80f, -2.11f,  5.66f,  5.97f,  8.23f,
           -6.05f, -3.30f,  5.36f, -4.44f,  1.08f,
           -0.45f,  2.58f, -2.70f,  0.27f,  9.04f,
            8.32f,  2.71f,  4.35f, -7.17f,  2.14f,
           -9.67f, -5.14f, -7.26f,  6.08f, -6.87f
        }  
     
     *A1 = {
            6.80f, -2.11f,  5.66f,  5.97f,  8.23f,
           -6.05f, -3.30f,  5.36f, -4.44f,  1.08f,
           -0.45f,  2.58f, -2.70f,  0.27f,  9.04f,
            8.32f,  2.71f,  4.35f, -7.17f,  2.14f,
           -9.67f, -5.14f, -7.26f,  6.08f, -6.87f
        }
      
     *B = {
            4.02f,  6.19f, -8.22f, -7.57f, -3.03f,
           -1.56f,  4.00f, -8.67f,  1.75f,  2.86f,
            9.81f, -4.09f, -4.57f, -8.61f,  8.99f
            } 
     
     X ={
            4.02f,  6.19f, -8.22f, -7.57f, -3.03f,
           -1.56f,  4.00f, -8.67f,  1.75f,  2.86f,
            9.81f, -4.09f, -4.57f, -8.61f,  8.99f
           } 
                       */
                       
     printf("Initializing CULA\n");
    status = culaInitialize();
    checkStatus(status);                  
     //printf("Results : CPU - Last element of matrix - %f / GPU - Last element of matrix - %f\n\n\n",A[1],A1[1]);
    printf("Calling culaSgesv...");
    gpu_time_1 = getTime(); 
    status = culaSgesv(n3, NRHS, A, n3, IPIV, X, n3);
    gpu_time_2 = getTime();
    checkStatus(status);
    printf("done\n");

    //verify_result();

    printf("Shutting down CULA\n");
    culaShutdown();

    printf("----------------------\n");
    printf("   CPU CLAPACK_SGESV  \n");
    printf("----------------------\n");

    

    printf("Starting clapack_sgesv ....");

    cpu_time_1 = getTime();
#if defined HAVE_MKL_LAPACK
/*
lapack_int LAPACKE_sgesv( int matrix_order, lapack_int n, lapack_int nrhs,
                          float* a, lapack_int lda, lapack_int* ipiv, float* b,
                          lapack_int ldb );
*/

    info = LAPACKE_sgesv(CblasColMajor, n3, NRHS, A1,  n3, IPIV, B, n3);
#else // we have GNU clapack.h
/*
GNU version:
 int clapack_sgesv(const enum CBLAS_ORDER Order, const int N, const int NRHS,
                   float *A, const int lda, int *ipiv, float *B, const int ldb)
*/

    info = clapack_sgesv(CblasColMajor, n3, NRHS, A1,  n3, IPIV, B, n3);
#endif
    cpu_time_2 = getTime();
    if (info != 0) { fprintf(stderr, "clapack_sgesv failure with error %d\n", info);}
    else { fprintf(stdout,"done.\n"); }


    //printf("\n\n");
        
        
        
           printf("-------------------\n");
           printf("   CHECK RESULT SGESV  \n");
           printf("-------------------\n");
           float pole = 0;
           float newes;
           //printf("\n");
           for (i=0;i<n3;i++)
            {
              newes = fabs(X[i]-B[i]);
              pole = pole + newes;
            if(*verbosity>2)  printf("Result check :  i=%d  GPU[%d]=%f  CPU[%d]=%f  Difference :%f \n",i, i+1,X[i],i+1, B[i],newes); 
     
      
            }
       if(*verbosity==2)printf("Results : CPU - First element of matrix - %f / GPU - First element of vector - %f\n",B[0],X[0]);
        if(*verbosity<=2)printf("Results : CPU - Last element of matrix - %f / GPU - Last element of vector - %f\n",B[n3-1],X[n3-1]);
       printf("Average absolute difference :%f \n", pole/n3);
        
        printf("-------------------\n");
        printf("   TIMES SGESV  \n");
        printf("-------------------\n");
        printf("Time of computation on CPU: %f sekund\n",(cpu_time_2-cpu_time_1)/1000);
        printf("Overall time of computation and transfer on GPU: %f sekund\n\n",(gpu_time_2-gpu_time_1)/1000);


    //verify_result();
    //printf("\n\n\n");
        //if(*verbosity>0)printf("Results : CPU - First element of matrix - %f / GPU - First element of matrix - %f\n",B[0],X[0]);
        //printf("Results : CPU - Last element of matrix - %f / GPU - Last element of matrix - %f\n\n\n",B[n3-1],X[n3-1]);
        
        //printf("Time of computation on CPU: %f sekund\n",(cpu_time_2-cpu_time_1)/1000);
        //printf("Overall time of computation and transfer on GPU: %f sekund\n\n\n",(gpu_time_2-gpu_time_1)/1000);

   // free(A);
    //free(B);
   // free(X);
   // free(IPIV);
}

void sgegrf_cula_c_test(int *n, int *verbosity){
        double gpu_time_1,gpu_time_2,gpu_time; // zaciatocny a konecny cas na gpu
        double cpu_time_1,cpu_time_2; // zaciatocny a konecny cas na cpu
        
        int i;
        
         culaFloat thresh = 1e-6f;
         culaFloat diff;
        
         srand ( time(NULL) );
        
         void checkStatus(culaStatus status)
{
    char buf[80];
    if(!status) return;
    culaGetErrorInfoString(status, culaGetErrorInfo(), buf, sizeof(buf));
    printf("%s\n", buf);
    culaShutdown();
    exit(EXIT_FAILURE);
}
        
        culaStatus status;
        int n2;
        int info = 0;
        int n3 = *n;
        //n3 = 4;
        culaFloat* A = NULL;
        culaFloat* TAU_r = NULL;
        float* B = 0;
        float* TAU_r1 = 0;
        float* WORK_r = 0;
        
        float h[16] = {
            1.0f, 3.0f, 2.0f, 1.0f,
            1.0f, 1.0f, 4.0f, 1.0f,
            1.0f, 3.0f, 4.0f, 1.0f,
            1.0f, 1.0f, 2.0f, -3.0f
        };
        float h1[16] = {
            1.0f, 1.0f, 1.0f, 1.0f,
            3.0f, 1.0f, 3.0f, 1.0f,
            2.0f, 4.0f, 4.0f, 2.0f,
            1.0f, 1.0f, 1.0f, -3.0f
        };
        
        n2 = n3*n3;
        
        
        void verify_result(void)
    {
    //printf("Verifying Result\n");
    int p = 0;
    int pocit = 0;
    int j,i;
        for (i=0;i<n3;i++)
        {
        for (j=0;j<n3;j++){
          
          diff =  B[pocit] + A[j*n3+i];
          
         if(diff < 0.0f)
            diff = -diff;
         // printf("Result check failed:  i=%d  B[i]=%f  A[i]=%f\n", i, B[pocit], A[j*n3+i]);  
        if(diff > thresh){
           if(*verbosity>500) printf("Result check failed:  i=%d  B[%d,%d]=%f  A[%d,%d]=%f\n", i, i+1,j+1,B[pocit],i+1,j+1, A[j*n3+i]); 
           p++;
        
        }
        
        
        
           
          //A[j*n3+i] = rand() / (float)RAND_MAX;//h[pocit];
          //B[pocit] = A[j*n3+i];
          pocit++;
        }
        }
        float newes;
        float pole = 0;
        //pocit = ;
        //printf("\n\n\n");
       /* pocit=0;
        
        for (i=0;i<n3;i++)
        {
        for (j=0;j<n3;j++){
          
           newes = fabs(B[pocit]-A[j*n3+i]);
           pole = pole + newes;
           if(*verbosity>2) printf("Result check :  i=%d  B[%d,%d]=%f  A[%d,%d]=%f  Difference :%f \n",i, i+1,j+1,B[pocit],i+1,j+1, A[j*n3+i],newes); 
     
        
           // printf("%d",pocit);
        
        
        
           
          //A[j*n3+i] = rand() / (float)RAND_MAX;//h[pocit];
          //B[pocit] = A[j*n3+i];
          pocit++;
        }
        }  */
        /*int ooo = 0;
        printf("\nGPU\n");
        
        for (i=0;i<n3;i++) {
         for (j=0;j<n3;j++)
        {
        printf(" %f ",A[j*n3+i]);//  A1[j*n3+i] = i+1+(j+1)*(j+1);//rand() / (float)RAND_MAX;//j+1+i+1;//rand() / (float)RAND_MAX;
         // A[j*n3+i] = A1[j*n3+i];
        if ((ooo+1)%n3 == 0) {printf("\n");}
        ooo++;
        }
       }
       
       printf("\nTau_r\n");
        for (i=0;i<n3;i++) {
          printf(" %f ",TAU_r[i]);
        }
       ooo = 0;
      
       printf("\nCPU\n");
       
       for (i=0;i<n3;i++) {
         for (j=0;j<n3;j++)
        {
        printf(" %f ",B[i*n3+j]);//  A1[j*n3+i] = i+1+(j+1)*(j+1);//rand() / (float)RAND_MAX;//j+1+i+1;//rand() / (float)RAND_MAX;
         // A[j*n3+i] = A1[j*n3+i];
        if ((ooo+1)%n3 == 0) {printf("\n");}
        ooo++;
        }
       }
               
        }
        
        printf("\nTau_r\n");
        for (i=0;i<n3;i++) {
          printf(" %f ",TAU_r1[i]);
        }
       */
       //if(*verbosity>2) printf("\n\n");
       int prem = 0;
       for (i=0;i<n3;i++) {
         
         for (j=prem;j<n3;j++)
         {
           newes = fabs(A[j*n3+i]-B[i*n3+j]);
           pole = pole + newes;
         if(*verbosity>2)  printf("R[GPU][%i][%i] %f R[CPU][%i][%i] %f  Diference : %f\n",i+1,j+1,A[j*n3+i],i+1,j+1,B[i*n3+j],newes);
         }
         prem ++;
        }
       //if(*verbosity>2) printf("\n\n");
       prem = 0;
        for (i=0;i<n3;i++) {
         
         for (j=0;j<prem;j++)
         {
           newes = fabs(A[j*n3+i]-B[i*n3+j]);
           pole = pole + newes;
          if(*verbosity>2) printf("Q[GPU][%i][%i] %f Q[CPU][%i][%i] %f  Diference : %f\n",i+1,j+1,A[j*n3+i],i+1,j+1,B[i*n3+j],newes);
         }
         prem ++;
        }
        //if(*verbosity>2)printf("\n\n");
        for (i=0;i<n3;i++) {
          newes = fabs(TAU_r1[i]-TAU_r[i]);
          pole = pole + newes;
         if(*verbosity>2) printf("TAU_r(GPU)[%i] %f TAU_r(CPU)[%i] %f Diference :%f\n",i,TAU_r[i],i,TAU_r1[i],newes);
        }
        
       printf("Average absolute difference :%f \n", pole/(n3*n3));
       printf("Result check failed number: %d from %d\n",p,pocit);
 
    } 
        printf("-------------------\n");
        printf("   GPU CULA_SGEQRF  \n");
        printf("-------------------\n");
        
        B = (float*)malloc(n2 * sizeof(B[0]));
        if (B == 0) {
                fprintf (stderr, "!!!! host memory allocation error (A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"B mtx allocated \n"); }
        
        TAU_r1 = (float*)malloc(n3 * sizeof(TAU_r1[0]));
        if (TAU_r1 == 0) {
                fprintf (stderr, "!!!! host memory allocation error (TAU_r1)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"TAU_r1 allocated \n"); }
        
        WORK_r = (float*)malloc(n3 * sizeof(WORK_r[0]));
        if (WORK_r == 0) {
                fprintf (stderr, "!!!! host memory allocation error (WORK_r)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"WORK_r allocated \n"); }
        
        A = (culaFloat*)malloc(n2*sizeof(culaFloat));
        TAU_r = (culaFloat*)malloc(n3*sizeof(culaFloat));
        
        if (A == NULL) {fprintf (stderr, "!!!! device memory allocation error (A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A mtx allocated on graphic card\n"); }
        
        if (TAU_r == NULL) {fprintf (stderr, "!!!! device memory allocation error (TAU_r)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"TAU_r allocated on graphic card\n"); }
        
       /* for (i=0;i<n2;i++)
        {
         // A[i] = rand() / (float)RAND_MAX;
          B[i] = h[i];//(float) i;
          //A[i] = h1[i];
        }*/
        int pocit = 0;
        int j;
        for (i=0;i<n3;i++)
        {
        for (j=0;j<n3;j++){
          
          A[j*n3+i] = i+1 +j+1;//2;//rand() / (float)RAND_MAX;//////h[pocit];
          B[pocit] = A[j*n3+i];
          pocit++;
        }
        }
           
        //printf("Results : CPU - Last element of matrix - %f / GPU - Last element of matrix - %f\n\n\n",B[n2-1],A[n2-1]);
        
        printf("Initializing CULA\n");
        status = culaInitialize();
        checkStatus(status);
        printf("Starting cula_sgeqrf ....");
        gpu_time_1 = getTime();
        status = culaSgeqrf(n3, n3, A, n3, TAU_r);
        gpu_time_2 = getTime();
        checkStatus(status);
         printf("done.\n");
        //printf("Shutting down CULA\n\n");
        culaShutdown();
        //printf("\n\n");
        printf("----------------------\n");
        printf("   CPU LAPACKE_sgeqrf (MKL)  \n");
        printf("----------------------\n");
      
        printf("Starting clapack_sgeqrf ....");
        
#if defined HAVE_MKL_LAPACK
        cpu_time_1 = getTime();
        info=LAPACKE_sgeqrf(LAPACK_ROW_MAJOR,n3, n3, B, n3, TAU_r1);
        cpu_time_2 = getTime();
#else
        printf("clapack_sgeqrf not present in GNU clapack.h ! returning\n");
        //return;
#endif
        
        printf("done\n");
        
        //for (i=0;i<n2;i++)
       // {
         // A[i] = rand() / (float)RAND_MAX;
        //  printf("%f %f\n",A[i],B[i]);
          //A[i] = h[i]//(float) i;
          //B[i] = A[i];
        //}
       // printf("\n\n");
       // for (i=0;i<n2;i++)
       // {
         // A[i] = rand() / (float)RAND_MAX;
       //   printf("%f %f\n",TAU_r1[i],TAU_r[i]);
          //A[i] = h[i]//(float) i;
          //B[i] = A[i];
       // }
       printf("-------------------\n");
       printf("   CHECK RESULT SGEMM  \n");
       printf("-------------------\n");
        verify_result();
        //printf(" ... done\n");
        
       // printf("Results : CPU - Last element of matrix - %f / GPU - Last element of matrix - %f\n\n\n",B[n2-1],A[n2-1]);
       // if(*verbosity>0)printf("Results : CPU - First element of matrix - %f / GPU - First element of matrix - %f\n",TAU_r1[0],TAU_r[0]);
       // printf("Results : CPU - Last element of matrix - %f / GPU - Last element of matrix - %f\n\n\n",TAU_r1[n3-1],TAU_r[n3-1]);
        if(*verbosity==2)printf("Results : CPU - First element of matrix - %f / GPU - First element of matrix - %f\n",B[0],A[0]);
        if(*verbosity<=2)printf("Results : CPU - Last element of matrix - %f / GPU - Last element of matrix - %f\n",B[n3-1],A[n3-1]);
       //printf("\nAverage absolute difference :%f \n\n", pole/n3);
        //printf("\n\n");
        printf("-------------------\n");
        printf("   TIMES SGEMM  \n");
        printf("-------------------\n");
        printf("Time of computation on CPU: %f sekund\n",(cpu_time_2-cpu_time_1)/1000);
        printf("Overall time of computation and transfer on GPU: %f sekund\n\n",(gpu_time_2-gpu_time_1)/1000);
}
#endif // of USE_CULA

void c_gpu_math_tests_(int *n, char *routine, int *verbosity)
{
printf("--- Testing of C based routines --- \n");
printf("entering routine name=%s, size of the problem n=%i, verbosity level=%i \n", trim_(routine), *n, *verbosity);
//printf("trimmed entering routine name=%s\n",trim(routine));
// DOPLNTE SEM VASE C-ovske RUTINY, PAN KUZMA
//printf(" DOPLNTE SEM VASE C-ovske RUTINY, PAN KUZMA !!!!\n");
char *rutina;
rutina = trim_(routine);  
/*switch (rutina){
     case "sgemm_cuda_c": sgemm_cuda_c_test(n,verbosity);
                          break;
} */
/*switch(routine){
 case 'sgemm_cuda_c': sgemm_cuda_c_test(n,verbosity);
                    break;
}         */
if (strcmp(rutina,"sgemm_cuda_c")==0) sgemm_cuda_c_test(n, verbosity);
if (strcmp(rutina,"sgemv_cuda_c")==0) sgemv_cuda_c_test(n, verbosity);
if (strcmp(rutina,"dgemm_cuda_c")==0) sgemm_cuda_c_test(n, verbosity);
#if defined USE_CULA
if (strcmp(rutina,"sgemm_cula_c")==0) sgemm_cula_c_test(n, verbosity);
if (strcmp(rutina,"sgemv_cula_c")==0) sgemv_cula_c_test(n, verbosity);
if (strcmp(rutina,"sgesv_cula_c")==0) sgesv_cula_c_test(n, verbosity);
if (strcmp(rutina,"sgeqrf_cula_c")==0) sgegrf_cula_c_test(n, verbosity);
#else
   printf("CULA not present, exiting.\n");
#endif

//if (strcmp(rutina,"dgemm_cuda_c")==0) printf("\n\n\ndgemm_cuda_c\n\n\n");

//printf("\n\n\n%s\n\n\n",rutina);

return;

}

