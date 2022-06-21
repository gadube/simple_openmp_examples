#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// print matrix
void printm(float *m, long N) {
    for (long i = 0; i < N; i++) {
        printf("[");
        for (long j = 0; j < N; j++) {
            printf("%7.2lf ", m[i * N + j]);
        }
        printf("]");
        printf("\n");
    }
    printf("\n\n");
}

// print vector
void printv(float *v, long N) {
    printf("[");
    for (long i = 0; i < N; i++) {
        printf("%7.2lf ", v[i]);
    }
    printf("]");
    printf("\n\n");
}

int main(int argc, char** argv) {

    if (argc != 2) {
        fprintf(stderr, "Usage: ./mxv N\n");
        exit(-1);
    }

    long N = atol(argv[1]);
    if (N <= 0) {
        fprintf(stderr, "N must be non-negative and non-zero.\n");
        exit(-1);
    }

    // set up matrix and vector
    float *m = (float *)malloc(sizeof(float) * N * N);
    float *v = (float *)malloc(sizeof(float) * N);
    float *res = (float *)calloc(N, sizeof(float));

    // populate matrix and vector
    for (long i = 0; i < N; i++) {
        for (long j = 0; j < N; j++) {
            m[i * N + j] = i * N + j;
        }
        v[i] = i;
    }

    clock_t start, end;
    double cpu_time_used;

    printf("Performing Matrix-Vector Multiplication\n");
    //actual computation
    start = clock();

    #pragma omp parallel for
    for (long i = 0; i < N; i++) {
        for (long j = 0; j < N; j++) {
            res[i] += m[i * N + j] * v[j];
        }
    }

    end = clock();
    cpu_time_used = ((double) (end - start) / CLOCKS_PER_SEC);

    /* printv(res,N); */

    printf("Done.\nTotal Time: %lf", cpu_time_used);

    // clean up
    free(res);
    free(v);
    free(m);

    return 0;
}

