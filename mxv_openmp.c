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

void verify(float *res, float *m, float *v, long N) {
    const float EPSILON = 0.000001;
    float *verify_res = (float *)calloc(N, sizeof(float));
    for (long i = 0; i < N; i++) {
        for (long j = 0; j < N; j++) {
            verify_res[i] += m[i * N + j] * v[j];
        }
    }

    for (long i = 0; i < N; i++) {
        if ((verify_res[i] - res[i]) >= EPSILON) {
            printf("Invalid at pos %ld, res = %0.10f, actual = %0.10f.\n", i, res[i], verify_res[i]);
            exit(-1);
        }
    }

    free(verify_res);
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
            m[i * N + j] = (i * N + j) * 1.234;
        }
        v[i] = i * 1.234;
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

    printf("Done.\nTotal Time: %lf\n", cpu_time_used);
    printf("Beginning Verification...\n");
    verify(res, m, v, N);
    printf("Verified.\n");

    // clean up
    free(res);
    free(v);
    free(m);

    return 0;
}

