
# Fedora Linux 
ขั้นตอนการติดตั้ง **OpenMP + OpenBLAS + LAPACK** บน Fedora แบบ "Researcher Grade" ครับ

-----

### ขั้นตอนที่ 1: เตรียมเครื่องมือ (Build Tools)

Fedora มี Group package ที่รวมของจำเป็นไว้ให้แล้ว (GCC, Make, Git, etc.) ติดตั้งทีเดียวจบครับ

```bash
# 1. อัปเดตระบบให้ล่าสุด
sudo dnf update -y

# 2. ติดตั้งชุดเครื่องมือพัฒนา (Development Tools)
# ชุดนี้จะมี GCC (OpenMP) มาให้เรียบร้อยครับ
sudo dnf group install "Development Tools" -y
```

-----

### ขั้นตอนที่ 2: ติดตั้ง Libraries (OpenBLAS & LAPACK)

ใน Fedora เราจะลงตัว **`-devel`** เพื่อให้ได้ Header files (`.h`) สำหรับการเขียนโปรแกรมครับ

```bash
# ติดตั้ง OpenBLAS และ LAPACK (พร้อม Header files)
sudo dnf install openblas-devel lapack-devel -y
```

> **เกร็ดความรู้:**
>
>   * **`openblas-devel`**: จะให้ไฟล์ `cblas.h` และ Library `libopenblas.so`
>   * **`lapack-devel`**: จะให้ไฟล์ `lapacke.h` (C Interface) และ `liblapacke.so`

-----

### ขั้นตอนที่ 3: Code ทดสอบ (Fedora HPC Test)

สร้างไฟล์ชื่อ `fedora_test.c` ครับ Code นี้จะทดสอบการทำงานร่วมกันของทั้ง 3 ส่วน:

```c
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <cblas.h>    // Header ของ OpenBLAS (จาก package openblas-devel)
#include <lapacke.h>  // Header ของ LAPACK (จาก package lapack-devel)

int main() {
    printf("==========================================\n");
    printf("   Fedora HPC Stack Test (dnf install)    \n");
    printf("==========================================\n");

    // --- 1. OpenMP Test ---
    printf("\n[1] OpenMP Check:\n");
    printf("    Detected Processors: %d\n", omp_get_num_procs());
    printf("    Max Threads: %d\n", omp_get_max_threads());

    #pragma omp parallel
    {
        #pragma omp single
        printf("    -> Parallel Region Active! (Report form thread %d)\n", omp_get_thread_num());
    }

    // --- 2. OpenBLAS Test (Matrix Multiplication) ---
    printf("\n[2] OpenBLAS Check (DGEMM):\n");
    int N = 2000;
    // จอง Memory (ใช้ malloc ธรรมดาเพื่อความง่าย)
    double *A = malloc(N*N*sizeof(double));
    double *B = malloc(N*N*sizeof(double));
    double *C = malloc(N*N*sizeof(double));

    // Initialize (Parallel)
    #pragma omp parallel for
    for(int i=0; i<N*N; i++) { A[i] = 1.0; B[i] = 2.0; C[i] = 0.0; }

    double start = omp_get_wtime();
    // C = alpha*A*B + beta*C
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                N, N, N, 1.0, A, N, B, N, 0.0, C, N);
    double end = omp_get_wtime();
    
    double gflops = (2.0*N*N*N) / (end - start) / 1e9;
    printf("    -> Matrix Size: %d x %d\n", N, N);
    printf("    -> Time: %.4f s | Perf: %.2f GFlops\n", end-start, gflops);

    // --- 3. LAPACK Test (Linear Solver Ax=b) ---
    printf("\n[3] LAPACK Check (DGESV - Solve Ax=b):\n");
    lapack_int n = 3, nrhs = 1, lda = 3, ldb = 3, info;
    lapack_int ipiv[3];
    
    // ระบบสมการ 3x3 (Row-Major เพราะใช้ LAPACKE)
    double A_sys[9] = { 
        3.0, 1.0, 1.0, 
        1.0, 3.0, 1.0, 
        1.0, 1.0, 3.0 
    };
    double b_sys[3] = { 5.0, 5.0, 5.0 }; // คำตอบต้องได้ [1, 1, 1]

    // เรียกใช้ LAPACKE_dgesv (C Interface ที่แถมมากับ lapack-devel)
    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, A_sys, lda, ipiv, b_sys, ldb);

    if(info == 0) {
        printf("    -> Solve Success!\n");
        printf("    -> Result x: [%.1f, %.1f, %.1f]\n", b_sys[0], b_sys[1], b_sys[2]);
    } else {
        printf("    -> Solve Failed! (Info: %d)\n", info);
    }

    free(A); free(B); free(C);
    return 0;
}
```

-----

### ขั้นตอนที่ 4: การ Compile และ Run ⚙️

บน Fedora ไฟล์ Library จะถูกเก็บไว้ใน Path มาตรฐาน (`/usr/lib64`) อยู่แล้ว ทำให้คำสั่ง Compile สั้นและสะอาดมากครับ

**คำสั่ง Compile:**

```bash
gcc fedora_test.c -o fedora_test \
    -O3 -fopenmp \
    -lopenblas -llapacke -lm
```

> **อธิบาย Flags:**
>
>   * `-fopenmp`: เปิดใช้งาน OpenMP
>   * `-lopenblas`: Link กับ OpenBLAS Library
>   * `-llapacke`: Link กับ LAPACK C Interface (สำคัญ\! ถ้าไม่ใส่จะเรียก `LAPACKE_dgesv` ไม่ได้)
>   * `-lm`: Link กับ Math Library

**คำสั่ง Run:**

```bash
./fedora_test
```

-----

### ข้อสังเกตสำหรับ "Fedora vs CentOS"

1.  **Header Files:** ใน CentOS คุณอาจต้องใช้ `dgesv_` (มี underscore และต้อง pass pointer) แต่ใน Fedora คุณสามารถใช้ `LAPACKE_dgesv` (ไม่มี underscore และใช้ง่ายกว่า) ได้เลย เพราะ `lapack-devel` ของ Fedora เตรียม C Wrapper มาให้ครบครับ
2.  **Path:** ไม่ต้องกำหนด `LD_LIBRARY_PATH` ให้วุ่นวาย เพราะ Fedora เอาของลง `/usr/lib64` ซึ่ง Linux หาเจอเองโดยอัตโนมัติครับ
