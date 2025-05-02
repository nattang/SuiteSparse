//------------------------------------------------------------------------------
// CHOLMOD/Demo/cholmod_dl_simple: simple demo program for CHOLMOD
//------------------------------------------------------------------------------

// CHOLMOD/Demo Module.  Copyright (C) 2005-2023, Timothy A. Davis,
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0+

//------------------------------------------------------------------------------

// Read in a real symmetric or complex Hermitian matrix from stdin in
// MatrixMarket format, solve Ax=b where b=[1 1 ... 1]', and print the residual.
//
// Usage: cholmod_dl_simple < matrixfile
//
// There are four versions of this demo:
// cholmod_di_simple:   double, int32
// cholmod_dl_simple:   double, int64
// cholmod_si_simple:   float, int32
// cholmod_sl_simple:   float, int64
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cholmod.h"
// #include "cuda_runtime.h"

/* parse parameters */
// copied from AMGX
int findParamIndex(char **argv, int argc, const char *parm)
{
    int count = 0;
    int index = -1;

    for (int i = 0; i < argc; i++)
    {
        if (strncmp(argv[i], parm, 100) == 0)
        {
            index = i;
            count++;
        }
    }
    if (count == 0 || count == 1)
    {
        return index;
    }
    else
    {
        printf("ERROR: parameter %s has been specified more than once, exiting\n", parm);
        exit(1);
    }

    return -1;
}

void print_cholmod_dense_values(cholmod_dense *A)
{
    if (A == NULL || A->x == NULL)
    {
        printf("Matrix is empty or uninitialized.\n");
        return;
    }

    size_t nrow = A->nrow;
    size_t ncol = A->ncol;
    size_t ld = A->d;

    double *x = (double *)A->x;
    for (size_t i = 0; i < nrow; i++)
    {
        for (size_t j = 0; j < ncol; j++)
        {
            printf("%10.6f ", x[i + j * ld]); // column-major access
        }
        printf("\n");
    }
}

void print_cholmod_sparse_values(cholmod_sparse *A)
{
    if (A == NULL)
    {
        printf("Matrix is NULL.\n");
        return;
    }

    size_t ncol = A->ncol, nrow = A->nrow;
    int is_packed = A->packed;

    printf("cholmod_sparse: %zu x %zu\n", nrow, ncol);
    printf("xtype: %d\n", A->xtype);
    printf("dtype: %d\n", A->dtype);
    printf("stype: %d\n", A->stype);

    SuiteSparse_long *p = (SuiteSparse_long *)A->p;
    SuiteSparse_long *i = (SuiteSparse_long *)A->i;
    int stype = A->stype; // +1 => upper, –1 => lower, 0 => unsymmetric

    for (size_t j = 0; j < ncol; j++)
    {
        SuiteSparse_long col_start = p[j];
        SuiteSparse_long col_end = is_packed ? p[j + 1]
                                             : p[j] + ((int64_t *)A->nz)[j];

        for (SuiteSparse_long k = col_start; k < col_end; k++)
        {
            SuiteSparse_long row = i[k];

            // skip the opposite triangle
            if (stype > 0 && row > (SuiteSparse_long)j)
                continue; // upper-only
            if (stype < 0 && row < (SuiteSparse_long)j)
                continue; // lower-only

            double *x = (double *)A->x;
            printf("(%ld, %zu): %.6f\n",
                   (long)(row + 1), // 1-based row
                   j + 1,           // 1-based col
                   x[k]);
        }
    }
}

void compute_error(cholmod_dense *X, cholmod_dense *expected)
{
    if (X->nrow != expected->nrow || X->ncol != expected->ncol)
    {
        fprintf(stderr, "Error: Matrix dimensions must match!\n");
        exit(EXIT_FAILURE);
    }

    size_t nrow = X->nrow;
    size_t ncol = X->ncol;
    size_t ld = X->d;
    double error = 0.0;

    double *x_data = (double *)(X->x);
    double *expected_data = (double *)(expected->x);

    // printf("X:\n");
    // print_cholmod_dense_values(X);

    // printf("Expected:\n");
    // print_cholmod_dense_values(expected);

    for (size_t j = 0; j < ncol; j++)
    {
        for (size_t i = 0; i < nrow; i++)
        {
            size_t idx = i + j * ld;
            error += fabs(x_data[idx] - expected_data[idx]);
        }
    }

    printf("Total error between RXMesh and CHOLMOD: %e\n", error);
}

cholmod_dense *multiply_sparse_dense(
    cholmod_sparse *A,
    cholmod_dense *X)
{
    int32_t *Ap = (int32_t *)(A->p);
    int32_t *Ai = (int32_t *)(A->i);
    double *Ax = (double *)(A->x);
    double *Xx = (double *)(X->x);

    size_t nrow = A->nrow;
    size_t ncol = A->ncol;
    size_t xcols = X->ncol;
    size_t xld = X->d;

    double *B_data = calloc(nrow * xcols, sizeof(double));
    if (!B_data)
    {
        perror("calloc");
        return NULL;
    }

    for (size_t col = 0; col < ncol; ++col)
    {
        for (int32_t k = Ap[col]; k < Ap[col + 1]; ++k)
        {
            int32_t row = Ai[k];
            double val = Ax[k];

            for (size_t j = 0; j < xcols; ++j)
            {
                // Apply symmetric contribution
                B_data[row + nrow * j] += val * Xx[col + xld * j];
                if (row != col)
                {
                    B_data[col + nrow * j] += val * Xx[row + xld * j];
                }
            }
        }
    }

    cholmod_dense *B = malloc(sizeof(cholmod_dense));
    B->nrow = nrow;
    B->ncol = xcols;
    B->nzmax = nrow * xcols;
    B->d = nrow;
    B->x = B_data;
    B->z = NULL;
    B->xtype = CHOLMOD_REAL;
    B->dtype = CHOLMOD_DOUBLE;

    printf("A times X:\n");
    print_cholmod_dense_values(B);
    free(B->x);

    return B;
}

int main(int argc, char **argv)
{
    // parameter parsing variables
    int pidx = 0;
    int pidy = 0;
    cholmod_sparse *A;
    cholmod_dense *X, *B, *r, *expected;
    cholmod_factor *L;
    double one[2] = {1, 0}, m1[2] = {-1, 0}; // basic scalars
    cholmod_common c;
    cholmod_l_start(&c);        // start CHOLMOD
    int dtype = CHOLMOD_DOUBLE; // use double precision
    int mtype = CHOLMOD_SPARSE;

    // reading from own matrices
    if ((pidx = findParamIndex(argv, argc, "-A")) != -1)
    {
        FILE *A_file = fopen(argv[pidx + 1], "r");
        if (A_file == NULL)
        {
            perror("Error opening A file");
            exit(1);
        }

        A = cholmod_l_read_sparse2(A_file, dtype, &c);
        // write to file
        FILE *f = fopen("./A.mtx", "w");
        if (f == NULL)
        {
            perror("Error opening B file for writing\n");
            exit(1);
        }
        cholmod_l_write_sparse(f, A, NULL, NULL, &c);
        fclose(f);
        printf("wrote A to file\n");

        // cholmod_triplet *T = cholmod_l_read_triplet(A_file, &c);
        fclose(A_file);

        // T->stype = 0;
        // A = cholmod_l_triplet_to_sparse(T, T->nnz, &c);
        // cholmod_l_free_triplet(&T, &c);

        // print_cholmod_sparse_values(A);
        c.precise = true;
        c.print = (A->nrow > 5) ? 3 : 5;
        cholmod_l_print_sparse(A, "A", &c); // print the matrix
        // print_cholmod_sparse_values(A);
    }
    else
    {
        printf("No matrix passed in");
        return 0;
    }

    if (A == NULL || A->stype == 0) // A must be symmetric
    {
        printf("A is not symmetric\n");
        cholmod_l_free_sparse(&A, &c);
        cholmod_l_finish(&c);
        return (0);
    }

    // TODO: populate B with all zeros, populate B with values
    B = cholmod_l_allocate_dense(A->nrow, 3, A->nrow, A->xtype + dtype, &c);
    // cudaPointerAttributes attributes;
    // cudaPointerGetAttributes(&attributes, (void *)B->x);
    if ((pidx = findParamIndex(argv, argc, "-B")) != -1)
    {
        FILE *B_file = fopen(argv[pidx + 1], "r");
        if (B_file == NULL)
        {
            perror("Error opening B file\n");
            exit(1);
        }

        printf("populating B with values from file\n");
        B = cholmod_l_read_dense2(B_file, dtype, &c);
        fclose(B_file);

        // write to file
        // FILE *f = fopen("./B.mtx", "w");
        // if (f == NULL)
        // {
        //     perror("Error opening B file for writing\n");
        //     exit(1);
        // }
        // cholmod_l_write_dense(f, B, NULL, &c);
        // fclose(f);
        // printf("wrote B to file\n");
        // print_cholmod_dense_values(B);
        // B = cholmod_l_zeros (A->nrow, 1, A->xtype + dtype, &c);
        // printf("pointer to B: %p, %p\n", (void *)B->x, &B);

        // printf("first item of B->x: %f\n", float(B->x[0]));
    }
    else
    {
        B = cholmod_l_ones(A->nrow, 1, A->xtype + dtype, &c); // b = ones(n,1)
    }

    if ((pidx = findParamIndex(argv, argc, "-EX")) != -1)
    {
        FILE *EX_file = fopen(argv[pidx + 1], "r");
        if (EX_file == NULL)
        {
            perror("Error opening EX file\n");
            exit(1);
        }

        expected = cholmod_l_read_dense2(EX_file, dtype, &c);
        fclose(EX_file);

        // FILE *f = fopen("./Expected.mtx", "w");
        // if (f == NULL)
        // {
        //     perror("Error opening B file for writing\n");
        //     exit(1);
        // }
        // cholmod_l_write_dense(f, expected, NULL, &c);
        // fclose(f);
        // printf("wrote expected to file\n");
    }

    X = cholmod_l_allocate_dense(A->nrow, 3, A->nrow, A->xtype + dtype, &c);

    double t1 = SUITESPARSE_TIME;
    L = cholmod_l_analyze(A, &c); // analyze
    t1 = SUITESPARSE_TIME - t1;
    double t2 = SUITESPARSE_TIME;
    cholmod_l_factorize(A, L, &c); // factorize
    t2 = SUITESPARSE_TIME - t2;
    double t3 = SUITESPARSE_TIME;

    X = cholmod_l_solve(CHOLMOD_A, L, B, &c); // solve Ax=b
    t3 = SUITESPARSE_TIME - t3;

    printf("analyze   time: %10.3f sec\n", t1);
    printf("factorize time: %10.3f sec\n", t2);
    printf("solve     time: %10.3f sec\n", t3);
    printf("total     time: %10.3f sec\n", t1 + t2 + t3);

    if (expected != NULL)
    {
        compute_error(X, expected);
    }

    // cholmod_dense *B_actual = multiply_sparse_dense(A, X);
    // free(B_actual->x);
    // free(B_actual);
   
    cholmod_l_print_factor(L, "L", &c); // print the factorization
    cholmod_l_print_dense(X, "X", &c);  // print the solution
    if (expected)
    {
    }
    r = cholmod_l_copy_dense(B, &c); // r = b
#ifndef NMATRIXOPS
    cholmod_l_sdmult(A, 0, m1, one, X, r, &c);      // r = r-Ax
    double rnorm = cholmod_l_norm_dense(r, 0, &c);  // compute inf-norm of r
    double anorm = cholmod_l_norm_sparse(A, 0, &c); // compute inf-norm of A
    printf("\n%s precision results:\n", dtype ? "single" : "double");
    printf("norm(b-Ax) %8.1e\n", rnorm);
    printf("norm(A)    %8.1e\n", anorm);
    double relresid = rnorm / anorm;
    printf("resid: norm(b-Ax)/norm(A) %8.1e\n", relresid);
    fprintf(stderr, "resid: norm(b-Ax)/norm(A) %8.1e\n", relresid);
#else
    printf("residual norm not computed (requires CHOLMOD/MatrixOps)\n");
#endif
    cholmod_l_free_factor(&L, &c); // free matrices
    cholmod_l_free_sparse(&A, &c);
    cholmod_l_free_dense(&r, &c);
    cholmod_l_free_dense(&X, &c);
    cholmod_l_free_dense(&B, &c);
    cholmod_l_free_dense(&expected, &c);
    cholmod_l_print_common("common", &c);
    printf("retrieving GPU stats\n");
    cholmod_l_gpu_stats(&c);
    printf("finishing CHOLMOD \n");
    double t4 = SUITESPARSE_TIME;
    cholmod_l_finish(&c); // finish CHOLMOD
    t4 = SUITESPARSE_TIME - t4;
    printf("finished CHOLMOD in : %10.3f sec\n", t4);
    return (0);
}
