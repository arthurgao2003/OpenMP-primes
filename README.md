# Parallel Segmented Sieve of Eratosthenes

A high-performance C implementation for generating prime numbers within a specific range $[M, N]$ using **OpenMP**. This project is engineered for maximum cache efficiency and minimal memory footprint, specifically tuned for high-performance computing (HPC) environments.

---

## 🚀 Key Features & Optimizations

### 1. Memory Efficiency (Bit Array)
Instead of using a standard `int` (4 bytes) or `char` (1 byte) array, this implementation utilizes a **Bit Array** to track primality.
* **8x Space Reduction:** Each number is represented by a single bit ($0$ for prime, $1$ for composite).
* **Odd-Only Sieve:** Only odd numbers are stored and processed. This immediately halves the memory requirement and the iteration count, as even numbers (except 2) are disregarded.

### 2. Cache-Conscious Design
The sieve is segmented to ensure that data being processed stays within the CPU's fastest memory layers.
* **Segmented Sieve:** The range is divided into small, manageable blocks.
* **L2 Cache Tuning:** The `BLOCK_SIZE` is set to **512KB**. This is specifically optimized for the L2 cache size of **Intel Cascade Lake Platinum 8268** chips, ensuring near-instant data access during the marking phase.

### 3. Parallelization with OpenMP
* **Domain Decomposition:** The range is split into independent blocks, allowing threads to work concurrently without the need for expensive locks or synchronization primitives.
* **Static Scheduling:** Given that each block contains an identical count of odd numbers, `schedule(static)` is used to minimize the overhead of thread management and workload distribution.

---

## 🛠 Prerequisites

* **Compiler:** `gcc`
* **Library Dependencies:** * `omp.h` (OpenMP)
    * `math.h` (Math library)
    * `stdint.h` (Standard integer types)

## 🔨 Compilation

Compile with the `-O3` flag for high-level compiler optimizations and `-fopenmp` to enable parallel processing:

```bash
gcc -O3 -fopenmp genprimes.c -o genprimes -lm
```

## 💻 Usage 
Compilation
The program requires three command-line arguments:
1. M: The start of the range (inclusive).
2. N: The end of the range (inclusive).
3. threads: The number of OpenMP threads to use.

Syntax:
./genprimes <M> <N> <threads>

Example (Find primes between 2 and 1M using 8 threads):
./genprimes 2 1000000 8

## 📊 Performance & Output
Execution Flow:
1. Console Output: 
   - Displays total prime count in range [M, N].
   - Displays high-resolution execution time via omp_get_wtime().

2. File Output:
   - Generates 'output.txt' in the local directory.
   - Contains a space-separated list of all primes found.

Performance Notes:
- Utilizing bitwise macros (BIT_SET, BIT_TEST) provides O(1) marking.
- Segmented approach minimizes L2 cache misses, which is critical 
  for scaling on Intel Cascade Lake architectures.
