#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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
        fprintf(stderr, "Usage: ./vpv N\n");
        exit(-1);
    }

    long N = atol(argv[1]);
    if (N <= 0) {
        fprintf(stderr, "N must be non-negative and non-zero.\n");
        exit(-1);
    }

    // set up matrix and vector
    float *v1 = (float *)malloc(sizeof(float) * N);
    float *v2 = (float *)malloc(sizeof(float) * N);
    float *res = (float *)calloc(N, sizeof(float));

    // populate matrix and vector
    #pragma omp parallel for
    for (long i = 0; i < N; i++) {
        v1[i] = i * 1.2;
        v2[i] = i * 3.4;
    }

    clock_t start, end;
    double cpu_time_used;

    printf("Performing Vector Addition\n");
    //actual computation
    start = clock();

    #pragma omp parallel for
    for (long i = 0; i < N; i++) {
        res[i] += v1[i] + v2[i];
    }

    end = clock();
    cpu_time_used = ((double) (end - start) / CLOCKS_PER_SEC);

    /* printv(res,N); */

    printf("Done.\nTotal Time: %lf", cpu_time_used);

    // clean up
    free(res);
    free(v1);
    free(v2);

    return 0;
}

