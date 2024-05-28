#define PRIME52_L0 0x0002631a5cf5d3ed
#define PRIME52_L1 0x000dea2f79cd6581
#define PRIME52_L2 0x000000000014def9
#define PRIME52_L3 0x0000000000000000
#define PRIME52_L4 0x0000100000000000

#define PRIME52_RR0 0x0009d265e952d13b
#define PRIME52_RR1 0x000d63c715bea69f
#define PRIME52_RR2 0x0005be65cb687604
#define PRIME52_RR3 0x0003dceec73d217f
#define PRIME52_RR4 0x000009411b7c309a

// __ALIGN64 static const uint64_t VPRIME52L0[8] = 
//     { PRIME52_L0, PRIME52_L0, PRIME52_L0, PRIME52_L0, 
//       PRIME52_L0, PRIME52_L0, PRIME52_L0, PRIME52_L0 };

// __ALIGN64 static const uint64_t VPRIME52L0[8] = 
//     { PRIME52_L0, PRIME52_L0, PRIME52_L0, PRIME52_L0, 
//       PRIME52_L0, PRIME52_L0, PRIME52_L0, PRIME52_L0 };

// __ALIGN64 static const uint64_t VPRIME52L0[8] = 
//     { PRIME52_L0, PRIME52_L0, PRIME52_L0, PRIME52_L0, 
//       PRIME52_L0, PRIME52_L0, PRIME52_L0, PRIME52_L0 };

// __ALIGN64 static const uint64_t VPRIME52L0[8] = 
//     { PRIME52_L0, PRIME52_L0, PRIME52_L0, PRIME52_L0, 
//       PRIME52_L0, PRIME52_L0, PRIME52_L0, PRIME52_L0 };

// __ALIGN64 static const uint64_t VPRIME52L0[8] = 
//     { PRIME52_L0, PRIME52_L0, PRIME52_L0, PRIME52_L0, 
//       PRIME52_L0, PRIME52_L0, PRIME52_L0, PRIME52_L0 };


// __ALIGN64 static const uint64_t VPRIME52L0[8] = 
//     { PRIME52_L0, PRIME52_L0, PRIME52_L0, PRIME52_L0, 
//       PRIME52_L0, PRIME52_L0, PRIME52_L0, PRIME52_L0 };

__ALIGN64 static const U64 L[5] = 
    {{PRIME52_L0, PRIME52_L0, PRIME52_L0, PRIME52_L0, PRIME52_L0, PRIME52_L0, PRIME52_L0, PRIME52_L0}, 
     {PRIME52_L1, PRIME52_L1, PRIME52_L1, PRIME52_L1, PRIME52_L1, PRIME52_L1, PRIME52_L1, PRIME52_L1}, 
     {PRIME52_L2, PRIME52_L2, PRIME52_L2, PRIME52_L2, PRIME52_L2, PRIME52_L2, PRIME52_L2, PRIME52_L2}, 
     {PRIME52_L3, PRIME52_L3, PRIME52_L3, PRIME52_L3, PRIME52_L3, PRIME52_L3, PRIME52_L3, PRIME52_L3}, 
     {PRIME52_L4, PRIME52_L4, PRIME52_L4, PRIME52_L4, PRIME52_L4, PRIME52_L4, PRIME52_L4, PRIME52_L4}};

__ALIGN64 static const U64 RR[5] = 
    {{PRIME52_RR0, PRIME52_RR0, PRIME52_RR0, PRIME52_RR0, PRIME52_RR0, PRIME52_RR0, PRIME52_RR0, PRIME52_RR0}, 
     {PRIME52_RR1, PRIME52_RR1, PRIME52_RR1, PRIME52_RR1, PRIME52_RR1, PRIME52_RR1, PRIME52_RR1, PRIME52_RR1}, 
     {PRIME52_RR2, PRIME52_RR2, PRIME52_RR2, PRIME52_RR2, PRIME52_RR2, PRIME52_RR2, PRIME52_RR2, PRIME52_RR2}, 
     {PRIME52_RR3, PRIME52_RR3, PRIME52_RR3, PRIME52_RR3, PRIME52_RR3, PRIME52_RR3, PRIME52_RR3, PRIME52_RR3}, 
     {PRIME52_RR4, PRIME52_RR4, PRIME52_RR4, PRIME52_RR4, PRIME52_RR4, PRIME52_RR4, PRIME52_RR4, PRIME52_RR4}};

__ALIGN64 static const U64 LFACTOR = {0x51da312547e1b, 0x51da312547e1b, 0x51da312547e1b, 0x51da312547e1b, 0x51da312547e1b, 0x51da312547e1b, 0x51da312547e1b, 0x51da312547e1b};

void scalarToVector(U64 pub52_mb8[5], const uint8_t bn[8][32]) 
{
    // __ALIGN64 U64 pub52_mb8[5];
    ifma_BNU_to_mb8((uint64_t (*)[8])pub52_mb8, (const uint8_t (*)[32])bn);
}

// const __m512i LFACTOR = _mm512_set1_epi64(0x51da312547e1b);

__INLINE void sc25519_montgomery_reduce_mul(U64 out[], const U64 a[], const U64 b[]) {
    // __m512i LFACTOR = _mm512_set1_epi64(0x51da312547e1b);

    // __m512i L[5] = {_mm512_set1_epi64(0x0002631a5cf5d3ed), 
    //                   _mm512_set1_epi64(0x000dea2f79cd6581), 
    //                   _mm512_set1_epi64(0x000000000014def9), 
    //                   _mm512_set1_epi64(0x0000000000000000), 
    //                   _mm512_set1_epi64(0x0000100000000000)};

    U64 r0, r1, r2, r3, r4, r5, r6, r7, r8, r9;
    U64 u0, u1, u2, u3, u4;
    U64 v0, v1, v2, v3, v4;

    __m512i carry = get_zero64();

    U64 *va = (U64*) a;
    U64 *vb = (U64*) b;
    U64 *vr = (U64*) out;

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

//r0;u0
    u0 = mul52lo(r0, LFACTOR);
    carry = fma52hi(_mm512_srli_epi64(fma52lo(r0, u0, L[0]), 52), u0, L[0]);
//r1;u1
    r1 = fma52lo(r1, u0, L[1]);

    r1 = add64(carry, r1);
    u1 = mul52lo(r1, LFACTOR);
    carry = fma52hi(_mm512_srli_epi64(fma52lo(r1, u1, L[0]), 52), u1, L[0]);
//r2;u2
    r2 = fma52hi(r2, u0, L[1]);
    
    r2 = fma52lo(r2, u0, L[2]);
    r2 =                        fma52lo(r2, u1, L[1]);

    r2 = add64(carry, r2);
    u2 = mul52lo(r2, LFACTOR);
    carry = fma52hi(_mm512_srli_epi64(fma52lo(r2, u2, L[0]), 52), u2, L[0]);
//r3;u3
    r3 = fma52hi(r3, u0, L[2]);
    r3 =                        fma52hi(r3, u1, L[1]);
    
    // r3 = fma52lo(r3, u0, L[3]);
    r3 =                        fma52lo(r3, u1, L[2]);
    r3 =                                                fma52lo(r3, u2, L[1]);

    r3 = add64(carry, r3);
    u3 = mul52lo(r3, LFACTOR);
    carry = fma52hi(_mm512_srli_epi64(fma52lo(r3, u3, L[0]), 52), u3, L[0]);
//r4;u4
    // r4 = fma52hi(r4, u0, L[3]);
    r4 =                        fma52hi(r4, u1, L[2]);
    r4 =                                                fma52hi(r4, u2, L[1]);

    r4 = fma52lo(r4, u0, L[4]);
    // r4 =                        fma52lo(r4, u1, L[3]);
    r4 =                                                fma52lo(r4, u2, L[2]);
    r4 =                                                                        fma52lo(r4, u3, L[1]);
    
    r4 = add64(carry, r4);
    u4 = mul52lo(r4, LFACTOR);
    carry = fma52hi(_mm512_srli_epi64(fma52lo(r4, u4, L[0]), 52), u4, L[0]);
//r5;v0
    r5 = fma52hi(r5, u0, L[4]);
    // r5 =                        fma52hi(r5, u1, L[3]);
    r5 =                                                fma52hi(r5, u2, L[2]);
    r5 =                                                                        fma52hi(r5, u3, L[1]);

    r5 =                        fma52lo(r5, u1, L[4]);
    // r5 =                                                fma52lo(r5, u2, L[3]);
    r5 =                                                                        fma52lo(r5, u3, L[2]);
    r5 =                                                                                                fma52lo(r5, u4, L[1]);

    r5 = add64(carry, r5);
    v0 = and64_const(r5, (1ULL << 52) - 1);
    carry = _mm512_srli_epi64(r5, 52);
//r6;v1
    r6 =                        fma52hi(r6, u1, L[4]);
    // r6 =                                                fma52hi(r6, u2, L[3]);
    r6 =                                                                        fma52hi(r6, u3, L[2]);
    r6 =                                                                                                fma52hi(r6, u4, L[1]);
    
    r6 =                                                fma52lo(r6, u2, L[4]);
    // r6 =                                                                        fma52lo(r6, u3, L[3]);
    r6 =                                                                                                fma52lo(r6, u4, L[2]);

    r6 = add64(carry, r6);
    v1 = and64_const(r6, (1ULL << 52) - 1);
    carry = _mm512_srli_epi64(r6, 52);

//r7;v2
    r7 =                                                fma52hi(r7, u2, L[4]);
    // r7 =                                                                        fma52hi(r7, u3, L[3]);
    r7 =                                                                                                fma52hi(r7, u4, L[2]);

    r7 =                                                                        fma52lo(r7, u3, L[4]);
    // r7 =                                                                                                fma52lo(r7, u4, L[3]);

    r7 = add64(carry, r7);
    v2 = and64_const(r7, (1ULL << 52) - 1);
    carry = _mm512_srli_epi64(r7, 52);

//r8;v3
    r8 =                                                                        fma52hi(r8, u3, L[4]);
    // r8 =                                                                                                fma52hi(r7, u4, L[3]);

    r8 =                                                                                                fma52lo(r8, u4, L[4]);

    r8 = add64(carry, r8);
    v3 = and64_const(r8, (1ULL << 52) - 1);
    carry = _mm512_srli_epi64(r8, 52);

//r9;v4
    r9 =                                                                                                fma52hi(r9, u4, L[4]);

    r9 = add64(carry, r9);
    v4 = and64_const(r9, (1ULL << 52) - 1);

    storeu64(&vr[0], v0);
    storeu64(&vr[1], v1);
    storeu64(&vr[2], v2);
    storeu64(&vr[3], v3);
    storeu64(&vr[4], v4);
}

__INLINE void sc25519_sub(U64 out[], const U64 a[], const U64 b[]) {

    // __m512i L[5] = {_mm512_set1_epi64(0x0002631a5cf5d3ed), 
    //                   _mm512_set1_epi64(0x000dea2f79cd6581), 
    //                   _mm512_set1_epi64(0x000000000014def9), 
    //                   _mm512_set1_epi64(0x0000000000000000), 
    //                   _mm512_set1_epi64(0x0000100000000000)};

    U64 *ca = (U64*) a;
    U64 *cb = (U64*) b;
    U64 *cr = (U64*) out;

    U64 borrow = get_zero64();
    uint32_t i;
    for(i = 0; i < 5; ++i)
    {
        borrow = sub64(ca[i], add64(cb[i], _mm512_srli_epi64(borrow, 63)));
        cr[i] = and64_const(borrow, (1ULL << 52) - 1);
    }

    U64 underflow_mask = sub64(_mm512_xor_epi64(_mm512_srli_epi64(borrow, 63), _mm512_set1_epi64(1)), _mm512_set1_epi64(1));

    borrow = get_zero64();
    for(i = 0; i < 5; ++i)
    {
        borrow = add64(add64(cr[i], _mm512_srli_epi64(borrow, 52)), and64(L[i], underflow_mask));
        cr[i] = and64_const(borrow, (1ULL << 52) - 1);
    }
}

__INLINE void sc25519_mul(U64 out[], const U64 a[], const U64 b[]) {
    U64 ab[5];
    // __m512i RR[5] = {_mm512_set1_epi64(0x0009d265e952d13b), 
    //                    _mm512_set1_epi64(0x000d63c715bea69f), 
    //                    _mm512_set1_epi64(0x0005be65cb687604), 
    //                    _mm512_set1_epi64(0x0003dceec73d217f), 
    //                    _mm512_set1_epi64(0x000009411b7c309a)};
    sc25519_montgomery_reduce_mul(ab, a, b);

    sc25519_montgomery_reduce_mul(out, ab, RR);
    // __m512i L[5] = {_mm512_set1_epi64(0x0002631a5cf5d3ed), 
    //                   _mm512_set1_epi64(0x000dea2f79cd6581), 
    //                   _mm512_set1_epi64(0x000000000014def9), 
    //                   _mm512_set1_epi64(0x0000000000000000), 
    //                   _mm512_set1_epi64(0x0000100000000000)};
    sc25519_sub(out, out, L);
}