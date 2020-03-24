#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/* 
g++ -Wall -I/usr/local/include -c gsl_inverse.cpp
g++ -L/usr/local/lib gsl_inverse.o -lgsl -lgslcblas -lm 
 */

void get_inverse (double *mA_inv, double *mA, int const N){

    double *mP = NULL; double *mP_inv = NULL;
    gsl_matrix_view m ;//= gsl_matrix_view_array (A, N, N);

    // initialize the inverse matrix to be solved
    double vb[N];
    int s; gsl_vector_view b; gsl_vector *x; gsl_permutation * p;

    for (int j =0; j<N; j++){

    // initialize matrix A from para /*每次要把A重新初始化成参数*mA指向后面的东西！！！！*/
    // 用A 是因为 m = matrix_view 求解逆矩阵每个列向量时会改变matrix的值！！！！！！！不能用mA！！！！！！！
        double *A = NULL;
        A = new double[N*N];
        mP = mA;
        for (int ii=0; ii<N; ii++ ){ 
            for (int jj=0; jj<N; jj++){
                A[ii*N+jj] = *(mP++);
            }
        }
        m = gsl_matrix_view_array(A,N,N);
        
        /* printf ("m = \n"); // check m == matrix A: right
        gsl_matrix_fprintf (stdout, &m.matrix, "%g");
        printf ("\n\n"); */
       
         for (int ib=0; ib<N; ib++){
            vb[ib] = (j==ib) ? 1. : 0.;
        }
        
        b = gsl_vector_view_array (vb, N); //行号: i=j
        x = gsl_vector_alloc (N);
        p = gsl_permutation_alloc (N);

        gsl_linalg_LU_decomp (&m.matrix, p, &s);

        gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);

        /*  
            printf ("x = \n");
            gsl_vector_fprintf (stdout, x, "%g");
            printf ("that is x%d\n\n",j); 
        */

        mP_inv = mA_inv + j;
        for (int i=0; i<N; i++){ 
            *(mA_inv + i*N + j) = gsl_vector_get(x,i);
            //printf("\n正在赋值mA_inv[%d][%d] = *pij = %f",i,j,*(mA_inv + i*N + j));

        }

        gsl_vector_free (x);
        gsl_permutation_free (p); 
        delete [] A;
    }
}

/* int main(){
    void get_inverse (double *mA_inv, double *mA, int const N);

    int const N = 6; int i,j;
    double *mA = NULL; double *mA_inv = NULL;
    mA = new double[N*N]; mA_inv = new double[N*N];
    // double aA[N*N]= {
    //    1,0,0,0,0,0,
    //    0,1,0,0,0,0,
    //    0,0,1,0,0,0,
    //   0,0,0,1,0,0,
    //    0,0,0,0,1,0,
    //   0,0,0,0,0,1,}; 
    double aA [N*N] = {
         1,-1,0,0,0,-1,
         0,1,-1,-1,0,0,
         0,0,0,1,-1,1,
         6,6,4,0,0,0,
         0,0,-4,4,5,0,
         0,-6,0,4,0,6};

    double *P = mA;
    for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            *(P++) = aA[i*N+j];
        }
    }
  
    get_inverse(mA_inv, mA, N);
    printf("\n*******************************************************\nC++\n");
    for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            printf ( "A_inv[%d][%d] = %f\n",i,j,mA_inv[i*N+j]);
        }
    }
  
    delete [] mA; delete [] mA_inv;
    return 0;
} */