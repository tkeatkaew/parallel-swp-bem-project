# CentOS 7 (Core) Server 
สำหรับงาน BEM (Boundary Element Method) ที่ต้องแก้สมการ Dense Matrix (Inverse Matrix operation) ผมขอแนะนำว่า **"อย่าใช้ BLAS/LAPACK มาตรฐานที่มากับ yum (Reference BLAS)"** ครับ เพราะมันช้ามาก (Unoptimized)

การติดตั้ง **OpenBLAS** แบบ Compile จาก Source Code ครับ เพราะมันจะ Auto-detect CPU `x86-64` ของคุณแล้ว Optimize ชุดคำสั่ง (AVX/AVX2) ให้โดยอัตโนมัติ แถม OpenBLAS ตัวนี้จะรวม LAPACK มาให้ในตัวแล้วด้วย จัดการง่ายกว่าเยอะครับ

นี่คือ Step-by-Step สำหรับเครื่อง CentOS 7:

### 1\. เตรียม Compiler และ Tools ที่จำเป็น

GCC บน CentOS 7 default คือ version 4.8.5 ซึ่งรองรับ OpenMP (v3.1) ได้ดีพอสมควรสำหรับงานทั่วไป แต่ถ้าอนาคตอยากได้ฟีเจอร์ใหม่ๆ อาจต้องลง `devtoolset` เพิ่ม (แต่เอาเบื้องต้นก่อนนะครับ)

```bash
# Update ระบบและติดตั้ง Development Tools
sudo yum update -y
sudo yum groupinstall "Development Tools" -y
sudo yum install git wget -y
```

### 2\. ติดตั้ง OpenBLAS (รวม BLAS + LAPACK Optimized)

เราจะไม่ใช้ `yum install blas-devel` แต่เราจะ Compile เองเพื่อให้มันรู้จัก Architecture ของเครื่อง server คุณจริงๆ

```bash
# เข้าไปที่โฟลเดอร์ tmp หรือที่ที่คุณใช้เก็บ source
cd /usr/local/src

# Clone OpenBLAS จาก GitHub
sudo git clone https://github.com/xianyi/OpenBLAS.git
cd OpenBLAS

# Compile (จะใช้เวลาสักพัก ขึ้นอยู่กับความแรง CPU)
# USE_OPENMP=1 เพื่อบอกให้ OpenBLAS ทำงานแบบ Parallel ได้ด้วย
sudo make USE_OPENMP=1

# Install ไปที่ /opt/OpenBLAS (แยกออกมาจะได้ไม่ปนกับ system lib)
sudo make install PREFIX=/opt/OpenBLAS
```

### 3\. Config Environment Variables

เพื่อให้ Compiler (gcc) หา Header file เจอ และให้ Runtime หา Library เจอ เราต้องบอก Path มันครับ (ขั้นตอนนี้สำคัญมาก ถ้าข้ามไปตอนรันจะเจอ error: `cannot open shared object file`)

ให้เพิ่มบรรทัดเหล่านี้ลงในไฟล์ `~/.bashrc` ของคุณ (หรือ `.bash_profile`):

```bash
# เปิดไฟล์ด้วย editor ที่ถนัด เช่น nano หรือ vi
nano ~/.bashrc
```

เพิ่มเนื้อหาด้านล่างนี้ไปท้ายไฟล์:

```bash
# OpenBLAS config for Research Work
export LD_LIBRARY_PATH=/opt/OpenBLAS/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=/opt/OpenBLAS/lib/pkgconfig:$PKG_CONFIG_PATH
export C_INCLUDE_PATH=/opt/OpenBLAS/include:$C_INCLUDE_PATH

# OpenMP config
export OMP_NUM_THREADS=4  # ตั้งค่า default thread (ปรับตามจำนวน Core จริงที่มี)
```

บันทึกไฟล์และ reload ค่าใหม่:

```bash
source ~/.bashrc
```

-----

### 4\. เขียนโค้ดทดสอบ (C + OpenMP + BLAS/LAPACK)

เรามาเขียนโค้ด C สั้นๆ เพื่อทดสอบว่า:

1.  **OpenMP** ทำงาน (เรียก Thread ได้)
2.  **BLAS** ทำงาน (เรียกฟังก์ชัน `cblas_dgemm` คูณเมทริกซ์ได้)

สร้างไฟล์ `bem_test.c`:

```c
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <cblas.h> // Header ของ OpenBLAS

int main() {
    // 1. ทดสอบ OpenMP
    int max_threads = omp_get_max_threads();
    printf("--- System Check ---\n");
    printf("Detected Max Threads (OpenMP): %d\n", max_threads);

    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        if (id == 0) {
            printf("Parallel Region working! (Reported by master thread)\n");
        }
    }

    // 2. ทดสอบ BLAS (Matrix Multiplication: C = A * B)
    // สมมติเมทริกซ์ขนาดเล็ก 2x2 เพื่อทดสอบการ Link Library
    int n = 2;
    double A[4] = {1.0, 2.0, 3.0, 4.0};
    double B[4] = {5.0, 6.0, 7.0, 8.0};
    double C[4] = {0.0, 0.0, 0.0, 0.0};

    // เรียกใช้ cblas_dgemm (Double Precision General Matrix Multiplication)
    // C = alpha*A*B + beta*C
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                n, n, n, 1.0, A, n, B, n, 0.0, C, n);

    printf("\n--- BLAS Check ---\n");
    printf("Result Matrix C:\n");
    printf("[ %f  %f ]\n", C[0], C[1]);
    printf("[ %f  %f ]\n", C[2], C[3]);
    
    // คำตอบที่ถูกควรเป็น:
    // [ 1*5+2*7  1*6+2*8 ] = [ 19  22 ]
    // [ 3*5+4*7  3*6+4*8 ] = [ 43  50 ]

    return 0;
}
```

### 5\. วิธี Compile และ Run

นี่คือจุดที่หลายคนตกม้าตายครับ การ Link library ต้องเรียงลำดับให้ถูก และต้องชี้ไปที่ OpenBLAS ที่เราลงไว้

คำสั่ง Compile:

```bash
gcc -o bem_test bem_test.c -fopenmp -L/opt/OpenBLAS/lib -lopenblas -lpthread -lm
```

  * `-fopenmp`: เปิดใช้งาน OpenMP directive
  * `-L/opt/OpenBLAS/lib`: บอกให้หา library ใน folder นี้
  * `-lopenblas`: Link กับ OpenBLAS (ซึ่งรวม LAPACK ไว้แล้ว)
  * `-lpthread`: จำเป็นสำหรับ OpenBLAS ในโหมด multi-thread
  * `-lm`: Math library มาตรฐาน

ทดลอง Run:

```bash
./bem_test
```

### ข้อแนะนำเพิ่มเติมจากประสบการณ์ส่วนตัว (Pro Tip)

ในงาน BEM นั้น Inverse Matrix เป็นงานที่กินพลังเครื่องที่สุด (`O(N^3)` complexity)

1.  **Lapack Function:** คุณจะได้ใช้ `dgesv` (แก้สมการ Ax=b) หรือ `dgetrf`/`dgetri` (หา Inverse) แน่นอน
2.  **Environment Variable:** เวลาคุณรันงานจริง อย่าลืมจูน `OMP_NUM_THREADS` ดีๆ ครับ ถ้าคุณรันโค้ดเดียว ให้ตั้งเท่ากับจำนวน Core จริง แต่ถ้าคุณรัน MPI ร่วมด้วย (Hybrid Parallel) ต้องระวังอย่าให้ Thread แย่ง CPU กันเองครับ

