// #include "x25519avxdef.h"

#ifndef __GNUC__
#pragma warning(disable:4013)
#endif

#define MASK47 ((1ULL << (255 - 52 * 4)) - 1)

__ALIGN64 static const uint64_t MASK47_[8] = {MASK47, MASK47, MASK47, MASK47,
                                     MASK47, MASK47, MASK47, MASK47};
__ALIGN64 static const uint64_t MOD_2_255_[8] = {19, 19, 19, 19, 19, 19, 19, 19};
__ALIGN64 static const uint64_t MOD_2_260_[8] = {19*32, 19*32, 19*32, 19*32,
                                        19*32, 19*32, 19*32, 19*32};
#define MASK_47 _mm512_loadu_si512(MASK47_)
#define MOD_2_255 _mm512_loadu_si512(MOD_2_255_)
#define MOD_2_260 _mm512_loadu_si512(MOD_2_260_)

#define ROUND_MUL_SRC_(I, J, S_LO, R_LO, S_HI, R_HI) \
    R_LO = fma52lo(S_LO, va1[I], vb1[J]); \
    R_HI = fma52hi(S_HI, va[I], vb[J]);

#define ROUND_MUL_(I, J, M0, M1) \
    ROUND_MUL_SRC_(I, J, M0, M0, M1, M1)

#define ROUND_MUL_SRC(I, J, S_LO, R_LO, S_HI, R_HI) \
    R_LO = fma52lo(S_LO, va[I], vb[J]); \
    R_HI = fma52hi(S_HI, va[I], vb[J]);

#define ROUND_MUL(I, J, M0, M1) \
    ROUND_MUL_SRC(I, J, M0, M0, M1, M1)

#define REDUCE_ROUND(R0, R1, R5) \
    r##R0 = fma52lo(r##R0, r##R5, MOD_2_260); \
    r##R1 = fma52lo( fma52hi(r##R1, r##R5, MOD_2_260), \
        _mm512_srli_epi64(r##R5, 52), MOD_2_260);

#define NORM(I, J) \
    r##J = add64(r##J, _mm512_srli_epi64(r##I, 52)); \
    r##I = and64_const(r##I, (1ULL << 52) - 1);
    //_mm512_and_epi64(a, _mm512_set1_epi64((long long)mask))

#define ROUND_MUL_SRC_A(I, J, S_LO, R_LO, S_HI, R_HI) \
    R_LO = fma52lo(S_LO, a##I, a##J); \
    R_HI = fma52hi(S_HI, a##I, a##J);

#define ROUND_MUL_A(I, J, M0, M1) \
    ROUND_MUL_SRC_A(I, J, M0, M0, M1, M1)

#define MASK_R4 ((1ULL << (255 - 52 * 4)) - 1)
static const uint64_t VMASK_R4[8] = {MASK_R4, MASK_R4, MASK_R4, MASK_R4,
    MASK_R4, MASK_R4, MASK_R4, MASK_R4};

#define MASK52 ((1ULL << 52) - 1)
static const uint64_t VMASK52[8] = {MASK52, MASK52, MASK52, MASK52,
    MASK52, MASK52, MASK52, MASK52};

#define REDUCE_ROUND_(R, R0, R1, R5) \
    R##R0 = fma52lo(R##R0, R##R5, MOD_2_260); \
    R##R1 = fma52lo( fma52hi(R##R1, R##R5, MOD_2_260), \
        _mm512_srli_epi64(R##R5, 52), MOD_2_260);

#define NORM_(R, I, J) \
    R##J = add64(R##J, _mm512_srli_epi64(R##I, 52)); \
    R##I = and64(R##I, _mm512_loadu_si512(VMASK52));

#define REDUCE_R4_N_R9(R)                                         \
    R##4 = fma52lo(R##4, R##9, MOD_2_260);                        \
    R##0 = fma52lo(R##0, _mm512_srli_epi64(R##4, 47), MOD_2_255);            \
    R##4 = and64(R##4, _mm512_loadu_si512(VMASK_R4));


////////////////////////////////////////////////////////////

__INLINE void ed25519_mul(U64 out[], const U64 a[], const U64 b[]) {
    U64 r0, r1, r2, r3, r4, r5, r6, r7, r8, r9;

    U64 *va = (U64*) a;
    U64 *vb = (U64*) b;
    U64 a1[5], b1[5];

    a1[0] = a[0];
    a1[1] = a[1];
    a1[2] = a[2];
    a1[3] = a[3];
    a1[4] = a[4];
    b1[0] = b[0];
    b1[1] = b[1];
    b1[2] = b[2];
    b1[3] = b[3];
    b1[4] = b[4];
    U64 *va1 = (U64*) a1;
    U64 *vb1 = (U64*) b1;
    // va1[0] = _mm512_load_epi64(&a[0]);
    // va1[1] = _mm512_load_epi64(&a[1]);
    // va1[2] = _mm512_load_epi64(&a[2]);
    // va1[3] = _mm512_load_epi64(&a[3]);
    // va1[4] = _mm512_load_epi64(&a[4]);
    // vb1[0] = _mm512_load_epi64(&b[0]);
    // vb1[1] = _mm512_load_epi64(&b[1]);
    // vb1[2] = _mm512_load_epi64(&b[2]);
    // vb1[3] = _mm512_load_epi64(&b[3]);
    // vb1[4] = _mm512_load_epi64(&b[4]);
    U64 *vr = (U64*) out;

    r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = r8 = r9 = get_zero64();

    // Full multiplication
    ROUND_MUL_(4, 4, r8, r9)
//_
    ROUND_MUL_(3, 0, r3, r4)
    ROUND_MUL_(1, 2, r3, r4)
    ROUND_MUL_(0, 3, r3, r4)
    ROUND_MUL_(2, 1, r3, r4)
    ROUND_MUL_(2, 2, r4, r5)
    ROUND_MUL_(0, 4, r4, r5)
    ROUND_MUL_(1, 3, r4, r5)
    ROUND_MUL_(3, 1, r4, r5)
    ROUND_MUL_(4, 0, r4, r5)
//_
    ROUND_MUL_(1, 4, r5, r6)
    ROUND_MUL_(2, 3, r5, r6)
    ROUND_MUL_(3, 2, r5, r6)
    ROUND_MUL_(4, 1, r5, r6)
    ROUND_MUL_(2, 4, r6, r7)
    ROUND_MUL_(3, 3, r6, r7)
    ROUND_MUL_(4, 2, r6, r7)
//_
    ROUND_MUL_(0, 0, r0, r1)
    ROUND_MUL_(0, 1, r1, r2)
    ROUND_MUL_(0, 2, r2, r3)
    ROUND_MUL_(1, 0, r1, r2)
    ROUND_MUL_(1, 1, r2, r3)
    ROUND_MUL_(2, 0, r2, r3)
    ROUND_MUL_(3, 4, r7, r8)
    ROUND_MUL_(4, 3, r7, r8)

    r4 = fma52lo(r4, r9, MOD_2_260);
    r0 = fma52lo(r0, _mm512_srli_epi64(r4, 47), MOD_2_255);
    r4 = and64(r4, MASK_47);

    REDUCE_ROUND(0, 1, 5);
    REDUCE_ROUND(1, 2, 6);
    REDUCE_ROUND(2, 3, 7);
    REDUCE_ROUND(3, 4, 8);

    // Normalize result
    NORM(0,1)
    NORM(1,2)
    NORM(2,3)
    NORM(3,4)
    // r0 = fma52lo(r0, _mm512_srli_epi64(r4, 47), MOD_2_255);
    // r4 = and64(r4, MASK_47);

    storeu64(&vr[0], r0);
    storeu64(&vr[1], r1);
    storeu64(&vr[2], r2);
    storeu64(&vr[3], r3);
    storeu64(&vr[4], r4);
}

/* SQR
c=0  (0,0)  
c=1  (0,1)  
c=2  (0,2)  (1,1)  
c=3  (0,3)  (1,2)  
c=4  (0,4)  (1,3)  (2,2)  
c=5  (1,4)  (2,3)  
c=6  (2,4)  (3,3)  
c=7  (3,4)  
c=8  (4,4)
*/

__INLINE void ed25519_sqr(U64 out[], const U64 a[]) {
    U64 r0, r1, r2, r3, r4, r5, r6, r7, r8, r9;

    U64 *va = (U64*) a;
    U64 *vb = (U64*) a;
    U64 *vr = (U64*) out;

    r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = r8 = r9 = get_zero64();

    // Square
    ROUND_MUL(0, 1, r1, r2)
    ROUND_MUL(0, 2, r2, r3)
    ROUND_MUL(0, 3, r3, r4)
    ROUND_MUL(0, 4, r4, r5)
    ROUND_MUL(1, 4, r5, r6)
    ROUND_MUL(2, 4, r6, r7)
    ROUND_MUL(3, 4, r7, r8)

    ROUND_MUL(1, 2, r3, r4)
    ROUND_MUL(1, 3, r4, r5)
    ROUND_MUL(2, 3, r5, r6)

    r1 = add64(r1, r1);
    r2 = add64(r2, r2);
    r3 = add64(r3, r3);
    r4 = add64(r4, r4);
    r5 = add64(r5, r5);
    r6 = add64(r6, r6);
    r7 = add64(r7, r7);
    r8 = add64(r8, r8);

    ROUND_MUL(0, 0, r0, r1)
    ROUND_MUL(1, 1, r2, r3)
    ROUND_MUL(2, 2, r4, r5)
    ROUND_MUL(3, 3, r6, r7)
    ROUND_MUL(4, 4, r8, r9)

    // Reduce r4 upper bits
    r4 = fma52lo(r4, r9, MOD_2_260);
    r0 = fma52lo(r0, _mm512_srli_epi64(r4, 47), MOD_2_255);
    r4 = and64(r4, MASK_47);

    REDUCE_ROUND(0, 1, 5);
    REDUCE_ROUND(1, 2, 6);
    REDUCE_ROUND(2, 3, 7);
    REDUCE_ROUND(3, 4, 8);

    // Normalize result
    NORM(0,1);
    NORM(1,2);
    NORM(2,3);
    NORM(3,4);
    // r0 = fma52lo(r0, _mm512_srli_epi64(r4, 47), MOD_2_255);
    // r4 = and64(r4, MASK_47);

    storeu64(&vr[0], r0);
    storeu64(&vr[1], r1);
    storeu64(&vr[2], r2);
    storeu64(&vr[3], r3);
    storeu64(&vr[4], r4);
}

static void MB_FUNC_NAME(ed25519_sqr_latency_)(U64 out[], const U64 a[], int count) {
    U64 r0, r1, r2, r3, r4, r5, r6, r7, r8, r9;
    U64 a0, a1, a2, a3, a4;
    U64 r4_1;
    int i;

    U64 *vr = (U64*) out;

    a0 = a[0];
    a1 = a[1];
    a2 = a[2];
    a3 = a[3];
    a4 = a[4];
    for (i = 0; i < count; ++i)
    {
        r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = r8 = r9 = get_zero64();
        r4_1 = get_zero64();

        // Square
        ROUND_MUL_A(0, 1, r1, r2)
        ROUND_MUL_A(0, 2, r2, r3)
        ROUND_MUL_A(0, 3, r3, r4_1)
        ROUND_MUL_A(0, 4, r4_1, r5)
        ROUND_MUL_A(1, 4, r5, r6)
        ROUND_MUL_A(2, 4, r6, r7)
        ROUND_MUL_A(3, 4, r7, r8)

        ROUND_MUL_A(1, 2, r3, r4)
        ROUND_MUL_A(1, 3, r4, r5)
        ROUND_MUL_A(2, 3, r5, r6)

        r1 = add64(r1, r1);
        r2 = add64(r2, r2);
        r3 = add64(r3, r3);

        r4 = add64(r4, r4_1);
        r4 = add64(r4, r4);

        r5 = add64(r5, r5);
        r6 = add64(r6, r6);
        r7 = add64(r7, r7);
        r8 = add64(r8, r8);

        ROUND_MUL_A(0, 0, r0, r1)
        ROUND_MUL_A(1, 1, r2, r3)
        ROUND_MUL_A(2, 2, r4, r5)
        ROUND_MUL_A(3, 3, r6, r7)
        ROUND_MUL_A(4, 4, r8, r9)

        r4 = fma52lo(r4, r9, MOD_2_260);
        r0 = fma52lo(r0, _mm512_srli_epi64(r4, 47), MOD_2_255);
        r4 = and64(r4, MASK_47);

        REDUCE_ROUND(0, 1, 5);
        REDUCE_ROUND(1, 2, 6);
        REDUCE_ROUND(2, 3, 7);
        REDUCE_ROUND(3, 4, 8);

        // Normalize result
        NORM(0,1)
        NORM(1,2)
        NORM(2,3)
        NORM(3,4)

        a0 = r0;
        a1 = r1;
        a2 = r2;
        a3 = r3;
        a4 = r4;
    }

    storeu64(&vr[0], r0);
    storeu64(&vr[1], r1);
    storeu64(&vr[2], r2);
    storeu64(&vr[3], r3);
    storeu64(&vr[4], r4);
}

__INLINE void ed25519_mul_dual(U64 out0[], U64 out1[],
                const U64 a0[], const U64 b0[],
                const U64 a1[], const U64 b1[]) {

    U64 r00, r01, r02, r03, r04, r05, r06, r07, r08, r09;
    U64 r10, r11, r12, r13, r14, r15, r16, r17, r18, r19;

    U64 *vr0 = (U64*) out0;
    U64 *vr1 = (U64*) out1;

    r00 = r01 = r02 = r03 = r04 = r05 = r06 = r07 = r08 = r09 = get_zero64();
    r10 = r11 = r12 = r13 = r14 = r15 = r16 = r17 = r18 = r19 = get_zero64();

    // Full multiplication
    U64 *va = (U64*) a0;
    U64 *vb = (U64*) b0;
    ROUND_MUL(4, 4, r08, r09)
    ROUND_MUL(3, 0, r03, r04)
    ROUND_MUL(1, 2, r03, r04)
    ROUND_MUL(0, 3, r03, r04)
    ROUND_MUL(2, 1, r03, r04)
    ROUND_MUL(2, 2, r04, r05)
    ROUND_MUL(0, 4, r04, r05)
    ROUND_MUL(1, 3, r04, r05)
    ROUND_MUL(3, 1, r04, r05)
    ROUND_MUL(4, 0, r04, r05)
    ROUND_MUL(1, 4, r05, r06)
    ROUND_MUL(2, 3, r05, r06)
    ROUND_MUL(3, 2, r05, r06)
    ROUND_MUL(4, 1, r05, r06)
    ROUND_MUL(2, 4, r06, r07)
    ROUND_MUL(3, 3, r06, r07)
    ROUND_MUL(4, 2, r06, r07)
    ROUND_MUL(0, 0, r00, r01)
    ROUND_MUL(0, 1, r01, r02)
    ROUND_MUL(0, 2, r02, r03)
    ROUND_MUL(1, 0, r01, r02)
    ROUND_MUL(1, 1, r02, r03)
    ROUND_MUL(2, 0, r02, r03)
    ROUND_MUL(3, 4, r07, r08)
    ROUND_MUL(4, 3, r07, r08)

    va = (U64*) a1;
    vb = (U64*) b1;
    ROUND_MUL(4, 4, r18, r19)
    ROUND_MUL(3, 0, r13, r14)
    ROUND_MUL(1, 2, r13, r14)
    ROUND_MUL(0, 3, r13, r14)
    ROUND_MUL(2, 1, r13, r14)
    ROUND_MUL(2, 2, r14, r15)
    ROUND_MUL(0, 4, r14, r15)
    ROUND_MUL(1, 3, r14, r15)
    ROUND_MUL(3, 1, r14, r15)
    ROUND_MUL(4, 0, r14, r15)
    ROUND_MUL(1, 4, r15, r16)
    ROUND_MUL(2, 3, r15, r16)
    ROUND_MUL(3, 2, r15, r16)
    ROUND_MUL(4, 1, r15, r16)
    ROUND_MUL(2, 4, r16, r17)
    ROUND_MUL(3, 3, r16, r17)
    ROUND_MUL(4, 2, r16, r17)
    ROUND_MUL(0, 0, r10, r11)
    ROUND_MUL(0, 1, r11, r12)
    ROUND_MUL(0, 2, r12, r13)
    ROUND_MUL(1, 0, r11, r12)
    ROUND_MUL(1, 1, r12, r13)
    ROUND_MUL(2, 0, r12, r13)
    ROUND_MUL(3, 4, r17, r18)
    ROUND_MUL(4, 3, r17, r18)

    REDUCE_R4_N_R9(r0)
    REDUCE_R4_N_R9(r1)

    REDUCE_ROUND_(r0, 0, 1, 5);
    REDUCE_ROUND_(r0, 1, 2, 6);
    REDUCE_ROUND_(r0, 2, 3, 7);
    REDUCE_ROUND_(r0, 3, 4, 8);

    REDUCE_ROUND_(r1, 0, 1, 5);
    REDUCE_ROUND_(r1, 1, 2, 6);
    REDUCE_ROUND_(r1, 2, 3, 7);
    REDUCE_ROUND_(r1, 3, 4, 8);

    // Normalize result
    NORM_(r0, 0,1)
    NORM_(r0, 1,2)
    NORM_(r0, 2,3)
    NORM_(r0, 3,4)

    NORM_(r1, 0,1)
    NORM_(r1, 1,2)
    NORM_(r1, 2,3)
    NORM_(r1, 3,4)

    storeu64(&vr0[0], r00);
    storeu64(&vr0[1], r01);
    storeu64(&vr0[2], r02);
    storeu64(&vr0[3], r03);
    storeu64(&vr0[4], r04);

    storeu64(&vr1[0], r10);
    storeu64(&vr1[1], r11);
    storeu64(&vr1[2], r12);
    storeu64(&vr1[3], r13);
    storeu64(&vr1[4], r14);
}

__INLINE void ed25519_sqr_dual(U64 out0[], U64 out1[],
                const U64 a0[], const U64 a1[]) {

    U64 r00, r01, r02, r03, r04, r05, r06, r07, r08, r09;
    U64 r10, r11, r12, r13, r14, r15, r16, r17, r18, r19;

    U64 *vr0 = (U64*) out0;
    U64 *vr1 = (U64*) out1;

    r00 = r01 = r02 = r03 = r04 = r05 = r06 = r07 = r08 = r09 = get_zero64();
    r10 = r11 = r12 = r13 = r14 = r15 = r16 = r17 = r18 = r19 = get_zero64();

    // Square
    U64 *va = (U64*) a0;
    U64 *vb = (U64*) a0;
    ROUND_MUL(0, 1, r01, r02)
    ROUND_MUL(0, 2, r02, r03)
    ROUND_MUL(0, 3, r03, r04)
    ROUND_MUL(0, 4, r04, r05)
    ROUND_MUL(1, 4, r05, r06)
    ROUND_MUL(2, 4, r06, r07)
    ROUND_MUL(3, 4, r07, r08)
    ROUND_MUL(1, 2, r03, r04)
    ROUND_MUL(1, 3, r04, r05)
    ROUND_MUL(2, 3, r05, r06)

    r01 = add64(r01, r01);
    r02 = add64(r02, r02);
    r03 = add64(r03, r03);
    r04 = add64(r04, r04);
    r05 = add64(r05, r05);
    r06 = add64(r06, r06);
    r07 = add64(r07, r07);
    r08 = add64(r08, r08);

    ROUND_MUL(0, 0, r00, r01)
    ROUND_MUL(1, 1, r02, r03)
    ROUND_MUL(2, 2, r04, r05)
    ROUND_MUL(3, 3, r06, r07)
    ROUND_MUL(4, 4, r08, r09)

    va = (U64*) a1;
    vb = (U64*) a1;
    ROUND_MUL(0, 1, r11, r12)
    ROUND_MUL(0, 2, r12, r13)
    ROUND_MUL(0, 3, r13, r14)
    ROUND_MUL(0, 4, r14, r15)
    ROUND_MUL(1, 4, r15, r16)
    ROUND_MUL(2, 4, r16, r17)
    ROUND_MUL(3, 4, r17, r18)
    ROUND_MUL(1, 2, r13, r14)
    ROUND_MUL(1, 3, r14, r15)
    ROUND_MUL(2, 3, r15, r16)

    r11 = add64(r11, r11);
    r12 = add64(r12, r12);
    r13 = add64(r13, r13);
    r14 = add64(r14, r14);
    r15 = add64(r15, r15);
    r16 = add64(r16, r16);
    r17 = add64(r17, r17);
    r18 = add64(r18, r18);

    ROUND_MUL(0, 0, r10, r11)
    ROUND_MUL(1, 1, r12, r13)
    ROUND_MUL(2, 2, r14, r15)
    ROUND_MUL(3, 3, r16, r17)
    ROUND_MUL(4, 4, r18, r19)

    REDUCE_R4_N_R9(r0)
    REDUCE_R4_N_R9(r1)

    REDUCE_ROUND_(r0, 0, 1, 5);
    REDUCE_ROUND_(r0, 1, 2, 6);
    REDUCE_ROUND_(r0, 2, 3, 7);
    REDUCE_ROUND_(r0, 3, 4, 8);

    REDUCE_ROUND_(r1, 0, 1, 5);
    REDUCE_ROUND_(r1, 1, 2, 6);
    REDUCE_ROUND_(r1, 2, 3, 7);
    REDUCE_ROUND_(r1, 3, 4, 8);

    // Normalize result
    NORM_(r0, 0,1)
    NORM_(r0, 1,2)
    NORM_(r0, 2,3)
    NORM_(r0, 3,4)

    NORM_(r1, 0,1)
    NORM_(r1, 1,2)
    NORM_(r1, 2,3)
    NORM_(r1, 3,4)

    storeu64(&vr0[0], r00);
    storeu64(&vr0[1], r01);
    storeu64(&vr0[2], r02);
    storeu64(&vr0[3], r03);
    storeu64(&vr0[4], r04);

    storeu64(&vr1[0], r10);
    storeu64(&vr1[1], r11);
    storeu64(&vr1[2], r12);
    storeu64(&vr1[3], r13);
    storeu64(&vr1[4], r14);
}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

__INLINE void fe52mb8_set(U64 out[], uint64_t value)
{
    storeu64(&out[0], set64((long long)value));
    storeu64(&out[1], get_zero64());
    storeu64(&out[2], get_zero64());
    storeu64(&out[3], get_zero64());
    storeu64(&out[4], get_zero64());
}

__INLINE void fe52mb8_copy(U64 out[], const U64 in[])
{
    storeu64(&out[0], _mm512_loadu_si512(&in[0]));
    storeu64(&out[1], _mm512_loadu_si512(&in[1]));
    storeu64(&out[2], _mm512_loadu_si512(&in[2]));
    storeu64(&out[3], _mm512_loadu_si512(&in[3]));
    storeu64(&out[4], _mm512_loadu_si512(&in[4]));
}

// Clang warning -Wunused-function
#if(0)
__INLINE void fe52mb8_mul_mod25519(U64 vr[], const U64 va[], const U64 vb[])
{
    U64 r0, r1, r2, r3, r4, r5, r6, r7, r8, r9;
    r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = r8 = r9 = get_zero64();

    // Full multiplication
    ROUND_MUL(4, 4, r8, r9)

    ROUND_MUL(3, 0, r3, r4)
    ROUND_MUL(1, 2, r3, r4)
    ROUND_MUL(0, 3, r3, r4)
    ROUND_MUL(2, 1, r3, r4)
    ROUND_MUL(2, 2, r4, r5)
    ROUND_MUL(0, 4, r4, r5)
    ROUND_MUL(1, 3, r4, r5)
    ROUND_MUL(3, 1, r4, r5)
    ROUND_MUL(4, 0, r4, r5)
    ROUND_MUL(1, 4, r5, r6)
    ROUND_MUL(2, 3, r5, r6)

    ROUND_MUL(3, 2, r5, r6)
    ROUND_MUL(4, 1, r5, r6)
    ROUND_MUL(2, 4, r6, r7)
    ROUND_MUL(3, 3, r6, r7)
    ROUND_MUL(4, 2, r6, r7)
    
    ROUND_MUL(0, 0, r0, r1)
    ROUND_MUL(0, 1, r1, r2)
    ROUND_MUL(0, 2, r2, r3)
    ROUND_MUL(1, 0, r1, r2)
    ROUND_MUL(1, 1, r2, r3)
    ROUND_MUL(2, 0, r2, r3)
    ROUND_MUL(3, 4, r7, r8)
    ROUND_MUL(4, 3, r7, r8)

    //REDUCE_ROUND(4, 5, 9);
    r4 = fma52lo(r4, r9, MOD_2_260); //r9 always contributes 0 to r5 (if input normalized?)
    REDUCE_ROUND(3, 4, 8);
    REDUCE_ROUND(2, 3, 7);
    REDUCE_ROUND(1, 2, 6);
    REDUCE_ROUND(0, 1, 5);

    // Reduce r4 upper bits
    r0 = fma52lo(r0, srli64(r4, 47), MOD_2_255);

    // Trim top r4 bits that were already reduced above
    r4 = and64(r4, MASK_47);

    // Normalize result
    NORM(0,1)
    NORM(1,2)
    NORM(2,3)
    NORM(3,4)

    storeu64(&vr[0], r0);
    storeu64(&vr[1], r1);
    storeu64(&vr[2], r2);
    storeu64(&vr[3], r3);
    storeu64(&vr[4], r4);
}

__INLINE void fe52mb8_sqr_mod25519(U64 out[], const U64 a[])
{
   fe52mb8_mul_mod25519(out, a, a);
}
#endif

__INLINE void fe52mb8_mul121666_mod25519(U64 vr[], const U64 va[])
{
    U64 multiplier = set64(121666);

    U64 r0, r1, r2, r3, r4, r5;
    r0 = r1 = r2 = r3 = r4 = r5 = get_zero64();

    // multiply
    r0 = fma52lo(r0, va[0], multiplier);
    r1 = fma52lo(r1, va[1], multiplier);
    r2 = fma52lo(r2, va[2], multiplier);
    r3 = fma52lo(r3, va[3], multiplier);
    r4 = fma52lo(r4, va[4], multiplier);

    r5 = fma52hi(r5, va[4], multiplier);
    r1 = fma52hi(r1, va[0], multiplier);
    r2 = fma52hi(r2, va[1], multiplier);
    r3 = fma52hi(r3, va[2], multiplier);
    r4 = fma52hi(r4, va[3], multiplier);

    // reduce
    REDUCE_ROUND(0, 1, 5);
    r0 = fma52lo(r0, _mm512_srli_epi64(r4, 47), MOD_2_255);

    // trim top r4 bits that were already reduced above
    r4 = and64(r4, MASK_47);

    // normalize
    NORM(0,1)
    NORM(1,2)
    NORM(2,3)
    NORM(3,4)

    storeu64(&vr[0], r0);
    storeu64(&vr[1], r1);
    storeu64(&vr[2], r2);
    storeu64(&vr[3], r3);
    storeu64(&vr[4], r4);
}

#define PRIME25519_LO  0x000FFFFFFFFFFFED
#define PRIME25519_MID 0x000FFFFFFFFFFFFF
#define PRIME25519_HI  0x00007FFFFFFFFFFF

// __ALIGN64 static const uint64_t prime25519[5] = {
//   PRIME25519_LO, PRIME25519_MID, PRIME25519_MID, PRIME25519_MID, PRIME25519_HI};

__ALIGN64 static const uint64_t VPRIME25519_LO[8] = 
    { PRIME25519_LO, PRIME25519_LO, PRIME25519_LO, PRIME25519_LO, 
      PRIME25519_LO, PRIME25519_LO, PRIME25519_LO, PRIME25519_LO };

__ALIGN64 static const uint64_t VPRIME25519_MID[8] = 
    { PRIME25519_MID, PRIME25519_MID, PRIME25519_MID, PRIME25519_MID, 
      PRIME25519_MID, PRIME25519_MID, PRIME25519_MID, PRIME25519_MID };

__ALIGN64 static const uint64_t VPRIME25519_HI[8] = 
    { PRIME25519_HI, PRIME25519_HI, PRIME25519_HI, PRIME25519_HI, 
      PRIME25519_HI, PRIME25519_HI, PRIME25519_HI, PRIME25519_HI };


__INLINE U64 cmov_U64(U64 a, U64 b, __mb_mask kmask)
{  return mask_mov64 (a, kmask, b); }

#define NORM_ASHIFTR(R, I, J) \
    R##J = add64(R##J, srai64(R##I, DIGIT_SIZE)); \
    R##I = and64(R##I, _mm512_loadu_si512(VMASK52));

#define NORM_LSHIFTR(R, I, J) \
    R##J = add64(R##J, _mm512_srli_epi64(R##I, DIGIT_SIZE)); \
    R##I = and64(R##I, _mm512_loadu_si512(VMASK52));

__INLINE void fe52mb8_add_mod25519(U64 vr[], const U64 va[], const U64 vb[])
{
    /* r = a+b */
    U64 r0 = add64(va[0], vb[0]);
    U64 r1 = add64(va[1], vb[1]);
    U64 r2 = add64(va[2], vb[2]);
    U64 r3 = add64(va[3], vb[3]);
    U64 r4 = add64(va[4], vb[4]);

    /* t = r-modulus (2^255-19) */
    U64 t0 = sub64(r0, _mm512_loadu_si512(VPRIME25519_LO ));
    U64 t1 = sub64(r1, _mm512_loadu_si512(VPRIME25519_MID));
    U64 t2 = sub64(r2, _mm512_loadu_si512(VPRIME25519_MID));
    U64 t3 = sub64(r3, _mm512_loadu_si512(VPRIME25519_MID));
    U64 t4 = sub64(r4, _mm512_loadu_si512(VPRIME25519_HI ));

    /* normalize r0, r1, r2, r3, r4 */
    NORM_LSHIFTR(r, 0,1)
    NORM_LSHIFTR(r, 1,2)
    NORM_LSHIFTR(r, 2,3)
    NORM_LSHIFTR(r, 3,4)

    /* normalize t0, t1, t2, t3, t4 */
    NORM_ASHIFTR(t, 0,1)
    NORM_ASHIFTR(t, 1,2)
    NORM_ASHIFTR(t, 2,3)
    NORM_ASHIFTR(t, 3,4)

    /* condition mask t4<0? (-1) : 0 */
    __mb_mask cmask = cmp64_mask(t4, get_zero64(), _MM_CMPINT_LT);

    storeu64(&vr[0], cmov_U64(t0, r0, cmask));
    storeu64(&vr[1], cmov_U64(t1, r1, cmask));
    storeu64(&vr[2], cmov_U64(t2, r2, cmask));
    storeu64(&vr[3], cmov_U64(t3, r3, cmask));
    storeu64(&vr[4], cmov_U64(t4, r4, cmask));
}

__INLINE void fe52mb8_sub_mod25519(U64 vr[], const U64 va[], const U64 vb[])
{
    /* r = a-b */
    U64 r0 = sub64(va[0], vb[0]);
    U64 r1 = sub64(va[1], vb[1]);
    U64 r2 = sub64(va[2], vb[2]);
    U64 r3 = sub64(va[3], vb[3]);
    U64 r4 = sub64(va[4], vb[4]);

    /* t = r+modulus (2^255-19) */
    U64 t0 = add64(r0, _mm512_loadu_si512(VPRIME25519_LO ));
    U64 t1 = add64(r1, _mm512_loadu_si512(VPRIME25519_MID));
    U64 t2 = add64(r2, _mm512_loadu_si512(VPRIME25519_MID));
    U64 t3 = add64(r3, _mm512_loadu_si512(VPRIME25519_MID));
    U64 t4 = add64(r4, _mm512_loadu_si512(VPRIME25519_HI ));

    /* normalize r0, r1, r2, r3, r4 */
    NORM_ASHIFTR(r, 0,1)
    NORM_ASHIFTR(r, 1,2)
    NORM_ASHIFTR(r, 2,3)
    NORM_ASHIFTR(r, 3,4)

    /* normalize t0, t1, t2, t3, t4 */
    NORM_ASHIFTR(t, 0,1)
    NORM_ASHIFTR(t, 1,2)
    NORM_ASHIFTR(t, 2,3)
    NORM_ASHIFTR(t, 3,4)

    /* condition mask r4<0? (-1) : 0 */
    __mb_mask cmask = cmp64_mask(r4, get_zero64(), _MM_CMPINT_LT);

    storeu64(&vr[0], cmov_U64(r0, t0, cmask));
    storeu64(&vr[1], cmov_U64(r1, t1, cmask));
    storeu64(&vr[2], cmov_U64(r2, t2, cmask));
    storeu64(&vr[3], cmov_U64(r3, t3, cmask));
    storeu64(&vr[4], cmov_U64(r4, t4, cmask));
}

__INLINE void fe52mb8_red_p25519(U64 vr[], const U64 va[])
{
   /* r = a-p */
   U64 r0 = sub64(va[0], _mm512_loadu_si512(VPRIME25519_LO));
   U64 r1 = sub64(va[1], _mm512_loadu_si512(VPRIME25519_MID));
   U64 r2 = sub64(va[2], _mm512_loadu_si512(VPRIME25519_MID));
   U64 r3 = sub64(va[3], _mm512_loadu_si512(VPRIME25519_MID));
   U64 r4 = sub64(va[4], _mm512_loadu_si512(VPRIME25519_HI));

   /* normalize r0, r1, r2, r3, r4 */
   NORM_ASHIFTR(r, 0, 1)
   NORM_ASHIFTR(r, 1, 2)
   NORM_ASHIFTR(r, 2, 3)
   NORM_ASHIFTR(r, 3, 4)

   /* condition mask r4<0? (-1) : 0 */
   __mb_mask cmask = cmp64_mask(r4, get_zero64(), _MM_CMPINT_LT);

   storeu64(&vr[0], cmov_U64(r0, va[0], cmask));
   storeu64(&vr[1], cmov_U64(r1, va[1], cmask));
   storeu64(&vr[2], cmov_U64(r2, va[2], cmask));
   storeu64(&vr[3], cmov_U64(r3, va[3], cmask));
   storeu64(&vr[4], cmov_U64(r4, va[4], cmask));
}

// #define USE_DUAL_MUL_SQR
// #define fe52_mul  fe52mb8_mul_mod25519
// #define fe52_sqr  fe52mb8_sqr_mod25519
#define fe52_mul  ed25519_mul
#define fe52_sqr  ed25519_sqr
#define fe52_add  fe52mb8_add_mod25519
#define fe52_sub  fe52mb8_sub_mod25519
#define fe52_mul121666  fe52mb8_mul121666_mod25519
#define fe52_sqr_power  MB_FUNC_NAME(ed25519_sqr_latency_)


/*
   Compute 1/z = z^(2^255 - 19 - 2)
   considering the exponent as
   2^255 - 21 = (2^5) * (2^250 - 1) + 11.
*/
__INLINE void fe52mb8_inv_mod25519(U64 out[], const U64 z[])
{
    __ALIGN64 U64 t0[5];
    __ALIGN64 U64 t1[5];
    __ALIGN64 U64 t2[5];
    __ALIGN64 U64 t3[5];

    /* t0 = z ** 2 */
    fe52_sqr(t0, z);

    /* t1 = t0 ** (2 ** 2) = z ** 8 */
    fe52_sqr(t1, t0);
    fe52_sqr(t1, t1);

    /* t1 = z * t1 = z ** 9 */
    fe52_mul(t1, z, t1);
    /* t0 = t0 * t1 = z ** 11 -- stash t0 away for the end. */
    fe52_mul(t0, t0, t1);

    /* t2 = t0 ** 2 = z ** 22 */
    fe52_sqr(t2, t0);

    /* t1 = t1 * t2 = z ** (2 ** 5 - 1) */
    fe52_mul(t1, t1, t2);

    /* t2 = t1 ** (2 ** 5) = z ** ((2 ** 5) * (2 ** 5 - 1)) */
    fe52_sqr_power(t2, t1, 5);

    /* t1 = t1 * t2 = z ** ((2 ** 5 + 1) * (2 ** 5 - 1)) = z ** (2 ** 10 - 1) */
    fe52_mul(t1, t2, t1);

    /* Continuing similarly... */

    /* t2 = z ** (2 ** 20 - 1) */
    fe52_sqr_power(t2, t1, 10);

    fe52_mul(t2, t2, t1);

    /* t2 = z ** (2 ** 40 - 1) */
    fe52_sqr_power(t3, t2, 20);

    fe52_mul(t2, t3, t2);

    /* t2 = z ** (2 ** 10) * (2 ** 40 - 1) */
    fe52_sqr_power(t2, t2, 10);

    /* t1 = z ** (2 ** 50 - 1) */
    fe52_mul(t1, t2, t1);

    /* t2 = z ** (2 ** 100 - 1) */
    fe52_sqr_power(t2, t1, 50);

    fe52_mul(t2, t2, t1);

    /* t2 = z ** (2 ** 200 - 1) */
    fe52_sqr_power(t3, t2, 100);

    fe52_mul(t2, t3, t2);

    /* t2 = z ** (2 ** 50) * (2 ** 200 - 1) */
    fe52_sqr_power(t2, t2, 50);

    /* t1 = z ** (2 ** 250 - 1) */
    fe52_mul(t1, t2, t1);

    /* t1 = z ** ((2 ** 5) * (2 ** 250 - 1)) */
    fe52_sqr_power(t1, t1, 5);

    /* Recall t0 = z ** 11; out = z ** (2 ** 255 - 21) */
    fe52_mul(out, t1, t0);
}

#define cswap_U64(a, b, kmask) { \
    U64 ta = mask_mov64((a), (kmask), (b)); \
    (b)    = mask_mov64((b), (kmask), (a)); \
    (a)    = ta; \
}

static void fe52mb8_cswap(U64 a[], U64 b[], __mb_mask k)
{
   cswap_U64(a[0], b[0], k)
   cswap_U64(a[1], b[1], k)
   cswap_U64(a[2], b[2], k)
   cswap_U64(a[3], b[3], k)
   cswap_U64(a[4], b[4], k)
}

#if 0
static void x25519_scalar_mul(U64 out[], U64 scalar[], U64 point[])
{
    __ALIGN64 U64 x1[5], x2[5], x3[5];
    __ALIGN64 U64        z2[5], z3[5];
    __ALIGN64 U64 tmp0[5], tmp1[5];

    fe52mb8_copy(x1, point);
    fe52mb8_set(x2, 1);
    fe52mb8_set(z2, 0);
    fe52mb8_copy(x3, x1);
    fe52mb8_set(z3, 1);

    /* read high and remove (zero) bit 63 */
    U64 e = loadu64(&scalar[3]);
    e = slli64(e, 1);

    __mb_mask swap = get_mask(0);
    int bitpos;
    for (bitpos=254; bitpos>= 0; bitpos--) {
        if(63==(bitpos%64))
            e = loadu64(&scalar[bitpos/64]);

        __mb_mask b = cmp64_mask(e, get_zero64(), _MM_CMPINT_LT);

        swap = mask_xor (swap, b);
        fe52mb8_cswap(x2, x3, swap);
        fe52mb8_cswap(z2, z3, swap);
        swap = b;
        fe52_sub(tmp0, x3, z3);
        fe52_sub(tmp1, x2, z2);
        fe52_add(x2, x2, z2); 
        fe52_add(z2, x3, z3);

        #ifdef USE_DUAL_MUL_SQR
            ed25519_mul_dual(z3, z2, x2, tmp0, z2, tmp1);
        #else
            fe52_mul(z3, x2, tmp0);
            fe52_mul(z2, z2, tmp1);
        #endif

        #ifdef USE_DUAL_MUL_SQR
            ed25519_sqr_dual(tmp0, tmp1, tmp1, x2);
        #else
            fe52_sqr(tmp0, tmp1);
            fe52_sqr(tmp1, x2);
        #endif

        fe52_add(x3, z3, z2);
        fe52_sub(z2, z3, z2);
        fe52_mul(x2, tmp1, tmp0);
        fe52_sub(tmp1, tmp1, tmp0);
        fe52_sqr(z2, z2);
        fe52_mul121666(z3, tmp1);
        fe52_sqr(x3, x3);
        fe52_add(tmp0, tmp0, z3);

        #ifdef USE_DUAL_MUL_SQR
            ed25519_mul_dual(z3, z2, x1, z2, tmp1, tmp0);
        #else
            fe52_mul(z3, x1, z2);
            fe52_mul(z2, tmp1, tmp0);
        #endif

        e = slli64(e, 1);
    }

    fe52mb8_inv_mod25519(z2, z2);
    fe52_mul(out, x2, z2);
}
#endif

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

__INLINE void ed25519_mul_dual_wonorm(U64 out0[], U64 out1[],
                const U64 a0[], const U64 b0[],
                const U64 a1[], const U64 b1[]) {

    U64 r00, r01, r02, r03, r04, r05, r06, r07, r08, r09;
    U64 r10, r11, r12, r13, r14, r15, r16, r17, r18, r19;

    U64 *vr0 = (U64*) out0;
    U64 *vr1 = (U64*) out1;

    r00 = r01 = r02 = r03 = r04 = r05 = r06 = r07 = r08 = r09 = get_zero64();
    r10 = r11 = r12 = r13 = r14 = r15 = r16 = r17 = r18 = r19 = get_zero64();

    // Full multiplication
    U64 *va = (U64*) a0;
    U64 *vb = (U64*) b0;
    ROUND_MUL(4, 4, r08, r09)
    ROUND_MUL(3, 0, r03, r04)
    ROUND_MUL(1, 2, r03, r04)
    ROUND_MUL(0, 3, r03, r04)
    ROUND_MUL(2, 1, r03, r04)
    ROUND_MUL(2, 2, r04, r05)
    ROUND_MUL(0, 4, r04, r05)
    ROUND_MUL(1, 3, r04, r05)
    ROUND_MUL(3, 1, r04, r05)
    ROUND_MUL(4, 0, r04, r05)
    ROUND_MUL(1, 4, r05, r06)
    ROUND_MUL(2, 3, r05, r06)
    ROUND_MUL(3, 2, r05, r06)
    ROUND_MUL(4, 1, r05, r06)
    ROUND_MUL(2, 4, r06, r07)
    ROUND_MUL(3, 3, r06, r07)
    ROUND_MUL(4, 2, r06, r07)
    ROUND_MUL(0, 0, r00, r01)
    ROUND_MUL(0, 1, r01, r02)
    ROUND_MUL(0, 2, r02, r03)
    ROUND_MUL(1, 0, r01, r02)
    ROUND_MUL(1, 1, r02, r03)
    ROUND_MUL(2, 0, r02, r03)
    ROUND_MUL(3, 4, r07, r08)
    ROUND_MUL(4, 3, r07, r08)

    va = (U64*) a1;
    vb = (U64*) b1;
    ROUND_MUL(4, 4, r18, r19)
    ROUND_MUL(3, 0, r13, r14)
    ROUND_MUL(1, 2, r13, r14)
    ROUND_MUL(0, 3, r13, r14)
    ROUND_MUL(2, 1, r13, r14)
    ROUND_MUL(2, 2, r14, r15)
    ROUND_MUL(0, 4, r14, r15)
    ROUND_MUL(1, 3, r14, r15)
    ROUND_MUL(3, 1, r14, r15)
    ROUND_MUL(4, 0, r14, r15)
    ROUND_MUL(1, 4, r15, r16)
    ROUND_MUL(2, 3, r15, r16)
    ROUND_MUL(3, 2, r15, r16)
    ROUND_MUL(4, 1, r15, r16)
    ROUND_MUL(2, 4, r16, r17)
    ROUND_MUL(3, 3, r16, r17)
    ROUND_MUL(4, 2, r16, r17)
    ROUND_MUL(0, 0, r10, r11)
    ROUND_MUL(0, 1, r11, r12)
    ROUND_MUL(0, 2, r12, r13)
    ROUND_MUL(1, 0, r11, r12)
    ROUND_MUL(1, 1, r12, r13)
    ROUND_MUL(2, 0, r12, r13)
    ROUND_MUL(3, 4, r17, r18)
    ROUND_MUL(4, 3, r17, r18)

    REDUCE_R4_N_R9(r0)
    REDUCE_R4_N_R9(r1)

    REDUCE_ROUND_(r0, 0, 1, 5);
    REDUCE_ROUND_(r0, 1, 2, 6);
    REDUCE_ROUND_(r0, 2, 3, 7);
    REDUCE_ROUND_(r0, 3, 4, 8);

    REDUCE_ROUND_(r1, 0, 1, 5);
    REDUCE_ROUND_(r1, 1, 2, 6);
    REDUCE_ROUND_(r1, 2, 3, 7);
    REDUCE_ROUND_(r1, 3, 4, 8);

    storeu64(&vr0[0], r00);
    storeu64(&vr0[1], r01);
    storeu64(&vr0[2], r02);
    storeu64(&vr0[3], r03);
    storeu64(&vr0[4], r04);

    storeu64(&vr1[0], r10);
    storeu64(&vr1[1], r11);
    storeu64(&vr1[2], r12);
    storeu64(&vr1[3], r13);
    storeu64(&vr1[4], r14);
}

__INLINE void fe52mb8_mul_mod25519_wonorm(U64 vr[], const U64 va[], const U64 vb[])
{
    U64 r0, r1, r2, r3, r4, r5, r6, r7, r8, r9;
    r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = r8 = r9 = get_zero64();

    // Full multiplication
    ROUND_MUL(4, 4, r8, r9)

    ROUND_MUL(3, 0, r3, r4)
    ROUND_MUL(1, 2, r3, r4)
    ROUND_MUL(0, 3, r3, r4)
    ROUND_MUL(2, 1, r3, r4)
    ROUND_MUL(2, 2, r4, r5)
    ROUND_MUL(0, 4, r4, r5)
    ROUND_MUL(1, 3, r4, r5)
    ROUND_MUL(3, 1, r4, r5)
    ROUND_MUL(4, 0, r4, r5)
    ROUND_MUL(1, 4, r5, r6)
    ROUND_MUL(2, 3, r5, r6)

    ROUND_MUL(3, 2, r5, r6)
    ROUND_MUL(4, 1, r5, r6)
    ROUND_MUL(2, 4, r6, r7)
    ROUND_MUL(3, 3, r6, r7)
    ROUND_MUL(4, 2, r6, r7)
    
    ROUND_MUL(0, 0, r0, r1)
    ROUND_MUL(0, 1, r1, r2)
    ROUND_MUL(0, 2, r2, r3)
    ROUND_MUL(1, 0, r1, r2)
    ROUND_MUL(1, 1, r2, r3)
    ROUND_MUL(2, 0, r2, r3)
    ROUND_MUL(3, 4, r7, r8)
    ROUND_MUL(4, 3, r7, r8)

    //REDUCE_ROUND(4, 5, 9);
    r4 = fma52lo(r4, r9, MOD_2_260); //r9 always contributes 0 to r5 (if input normalized?)
    REDUCE_ROUND(3, 4, 8);
    REDUCE_ROUND(2, 3, 7);
    REDUCE_ROUND(1, 2, 6);
    REDUCE_ROUND(0, 1, 5);

    // Reduce r4 upper bits
    r0 = fma52lo(r0, _mm512_srli_epi64(r4, 47), MOD_2_255);

    // Trim top r4 bits that were already reduced above
    r4 = and64(r4, MASK_47);

    storeu64(&vr[0], r0);
    storeu64(&vr[1], r1);
    storeu64(&vr[2], r2);
    storeu64(&vr[3], r3);
    storeu64(&vr[4], r4);
}

__INLINE void fe52mb8_mul121666_mod25519_wonorm(U64 vr[], const U64 va[])
{
    U64 multiplier = set64(121666);

    U64 r0, r1, r2, r3, r4, r5;
    r0 = r1 = r2 = r3 = r4 = r5 = get_zero64();

    // multiply
    r0 = fma52lo(r0, va[0], multiplier);
    r1 = fma52lo(r1, va[1], multiplier);
    r2 = fma52lo(r2, va[2], multiplier);
    r3 = fma52lo(r3, va[3], multiplier);
    r4 = fma52lo(r4, va[4], multiplier);

    r5 = fma52hi(r5, va[4], multiplier);
    r1 = fma52hi(r1, va[0], multiplier);
    r2 = fma52hi(r2, va[1], multiplier);
    r3 = fma52hi(r3, va[2], multiplier);
    r4 = fma52hi(r4, va[3], multiplier);

    // reduce
    REDUCE_ROUND(0, 1, 5);
    r0 = fma52lo(r0, _mm512_srli_epi64(r4, 47), MOD_2_255);

    // trim top r4 bits that were already reduced above
    r4 = and64(r4, MASK_47);

    storeu64(&vr[0], r0);
    storeu64(&vr[1], r1);
    storeu64(&vr[2], r2);
    storeu64(&vr[3], r3);
    storeu64(&vr[4], r4);
}