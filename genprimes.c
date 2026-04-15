#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <stdint.h>

/* 
   Key data strucutre: Bit Array so that we can store each number using 1 bit instead of 1 byte or int. This takes up much less cache space and will result in better cache performance
   My implementation of representation of primes: 
    Each bit represents a number and gives the current result of prime calculation:
        0:  "still prime candidate"
        1:  "composite (not prime)" */

/*three helper methods: bit set marks i as a composite number, bit_clear marks i as a prime candidate, and bit test checks whether the number is marked as composite already*/
#define BIT_SET(arr, i)   ((arr)[(i) >> 3] |=  (1 << ((i) & 7)))
#define BIT_CLEAR(arr, i) ((arr)[(i) >> 3] &= ~(1 << ((i) & 7)))
#define BIT_TEST(arr, i)  ((arr)[(i) >> 3] &   (1 << ((i) & 7)))

/* From online resource we know Greene all cluster nodes are equipped with 24-core Intel Cascade Lake Platinum 8268 chips, we want to make sure that our data fits entirely in the l2 cache so I choose block size to be 512kb*/

#define BLOCK_SIZE (1 << 19)

/* This function computes all primes up to sqrt(N) using Sieve of Eratosthenes
   Output: base_primes[] : list of primes used to eliminate multiples*/

static int simple_sieve(long long limit, int **base_primes) {

    // Create bit array for marking composites
    size_t bytes = (limit + 8) / 8;
    uint8_t *mark = calloc(bytes, 1);

    if (!mark) {
        perror("calloc");
        exit(1);
    }

    // Mark multiples of each number: check if multiple of each prime is already marked as a composite number
    for (long long i = 2; i * i <= limit; i++) {
        if (!BIT_TEST(mark, i)) {
            for (long long j = i * i; j <= limit; j += i) {
                BIT_SET(mark, j);
            }
        }
    }

    // Count how many primes exist up to sqrt(N)
    int count = 0;
    for (long long i = 2; i <= limit; i++) {
        if (!BIT_TEST(mark, i)) count++;
    }

    // Store primes into output array
    *base_primes = malloc(count * sizeof(int));

    int idx = 0;
    for (long long i = 2; i <= limit; i++) {
        if (!BIT_TEST(mark, i)) {
            (*base_primes)[idx++] = (int)i;
        }
    }

    free(mark);
    return count;
}

/* 
   main function
   main optimizations are: 
    1. storing only odd numbers to reduce memory storage space needed
    2. efficient memory allocation using bit arrays
    3. applying the segmented sieve solution in parallel using OpenMP. */

int main(int argc, char *argv[]) {
    /*error handling if insufficient arguments passed */
    double tstart = 0.0, ttaken;
    if (argc < 4) {
        fprintf(stderr, "Usage: %s M N threads\n", argv[0]);
        return 1;
    }

    /* handling user inputs*/
    long long M = atoll(argv[1]);
    long long N = atoll(argv[2]);
    int t = atoi(argv[3]);

    /* 0 and 1 are not prime */
    if (M < 2) M = 2;

    /* Start timer */
    tstart = omp_get_wtime();

    /* build base primes up to sqrt(N) */
    long long limit = (long long)sqrt((double)N) + 1;
    int *base_primes;
    int base_count = simple_sieve(limit, &base_primes);

    /* Storing only even numbers because no even number other than 2 can be prime */

    long long M_odd = (M % 2 == 0) ? M + 1 : M;
    long long odd_count = (N - M_odd) / 2 + 1;

    /* Allocate bit array for odd numbers only */
    size_t sieve_bytes = (odd_count + 7) / 8;
    uint8_t *sieve = calloc(sieve_bytes, 1);

    if (!sieve) {
        perror("calloc sieve");
        return 1;
    }

    /* OpenMP part
       The idea is that we divide the array into blocks so threads work independently.
       Because for each iteration (one block) we have the same number of elements (BLOCK_SIZE)
        and same operations (loop over base primes), using static scheduling here is optimal with minimal overhead*/

    long long total_blocks = (odd_count + BLOCK_SIZE - 1) / BLOCK_SIZE;

    #pragma omp parallel for num_threads(t) schedule(static)
    for (long long blk = 0; blk < total_blocks; blk++) {

        /* Compute range of indices for this block */
        long long lo = blk * BLOCK_SIZE;
        long long hi = lo + BLOCK_SIZE;
        if (hi > odd_count) hi = odd_count;

        /* Convert index range to actual numbers */
        long long num_lo = M_odd + 2 * lo;
        long long num_hi = M_odd + 2 * (hi - 1);

        /* Mark all multiples */

        for (int pi = 0; pi < base_count; pi++) {

            long long p = base_primes[pi];

            if (p == 2) continue; // even nums already removed

            /* Find first multiple of p in block */
            long long start = ((num_lo + p - 1) / p) * p;

            /* Ensure we start from p squared (nums below that already checked) */
            if (start < p * p) start = p * p;

            /* Ensure we only mark odd multiples */
            if (start % 2 == 0) start += p;

            /* Mark all multiples of p in this block */
            for (long long j = start; j <= num_hi; j += 2 * p) {

                long long idx = (j - M_odd) / 2;

                BIT_SET(sieve, idx);
            }
        }
    }

    /* Counting primes */

    long long total_count = 0;

    /* Include 2 if in range */
    if (M <= 2 && N >= 2) total_count++;

    /* Count odd primes */
    for (long long i = 0; i < odd_count; i++) {
        if (!BIT_TEST(sieve, i)) {
            total_count++;
        }
    }

    /* Allocate output array */
    long long *primes = malloc(total_count * sizeof(long long));

    long long out = 0;

    if (M <= 2 && N >= 2) {
        primes[out++] = 2;
    }

    /* Collect primes */
    for (long long i = 0; i < odd_count; i++) {
        if (!BIT_TEST(sieve, i)) {
            primes[out++] = M_odd + 2 * i;
        }
    }

    ttaken = omp_get_wtime() - tstart;

    /* File Output */

    FILE *f = fopen("output.txt", "w");

    if (!f) {
        printf("Error opening file\n");
        return 1;
    }

    for (long long i = 0; i < total_count; i++) {
        fprintf(f, "%lld ", primes[i]);
    }

    fclose(f);

    /* giving the output */

    printf("Primes found in range [%lld, %lld]: %lld\n", M, N, total_count);
    printf("Time taken: %f seconds\n", ttaken);

    free(base_primes);
    free(sieve);
    free(primes);

    return 0;
}