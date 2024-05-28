// #include "fieldarithavx.h"

void printU64(U64 out[5])
    {
        __m128i re;
        for (int i = 0; i < 5; i++) {

        re = _mm512_extracti64x2_epi64(out[i], 0);
        printf("%0llx ", re[0]);
        printf("%0llx ", re[1]);
        re = _mm512_extracti64x2_epi64(out[i], 1);
        printf("%0llx ", re[0]);
        printf("%0llx ", re[1]);
        re = _mm512_extracti64x2_epi64(out[i], 2);
        printf("%0llx ", re[0]);
        printf("%0llx ", re[1]);
        re = _mm512_extracti64x2_epi64(out[i], 3);
        printf("%0llx ", re[0]);
        printf("%0llx ", re[1]);
        printf("\n");
        }
    }

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

        __mb_mask b = _mm512_cmp_epi64_mask(e, get_zero64(), _MM_CMPINT_LT);

        swap = mask_xor (swap, b);
        fe52mb8_cswap(x2, x3, swap);
        fe52mb8_cswap(z2, z3, swap);
        swap = b;
        fe52_sub(tmp0, x3, z3);
        fe52_sub(tmp1, x2, z2);
        fe52_add(x2, x2, z2); 
        fe52_add(z2, x3, z3);

        // #ifdef USE_DUAL_MUL_SQR
            // ed25519_mul_dual(z3, z2, x2, tmp0, z2, tmp1);
        // #else
            fe52_mul(z3, x2, tmp0);
            fe52_mul(z2, z2, tmp1);
        // #endif

        // #ifdef USE_DUAL_MUL_SQR
            // ed25519_sqr_dual(tmp0, tmp1, tmp1, x2);
        // #else
            fe52_sqr(tmp0, tmp1);
            fe52_sqr(tmp1, x2);
        // #endif

        fe52_add(x3, z3, z2);
        fe52_sub(z2, z3, z2);
        fe52_mul(x2, tmp1, tmp0);
        fe52_sub(tmp1, tmp1, tmp0);
        fe52_sqr(z2, z2);
        fe52_mul121666(z3, tmp1);
        fe52_sqr(x3, x3);
        fe52_add(tmp0, tmp0, z3);

        // #ifdef USE_DUAL_MUL_SQR
            // ed25519_mul_dual(z3, z2, x1, z2, tmp1, tmp0);
        // #else
            fe52_mul(z3, x1, z2);
            fe52_mul(z2, tmp1, tmp0);
        // #endif

        e = slli64(e, 1);
    }

    fe52mb8_cswap(x2, x3, swap);
    fe52mb8_cswap(z2, z3, swap);

    fe52mb8_inv_mod25519(z2, z2);
    fe52_mul(out, x2, z2);
}

// static void
// curve25519_scalarmult_motgomery_fs(curve25519_key mypublic, const curve25519_key n, const curve25519_key point) {
// 	curve25519_key t;
// 	bignum25519 x2 = {1}, z2 = {0}, z3 = {1}, x3, x1;
// 	bignum25519 t0, t1, t2, t3, e0, e1, e2, e3;
//     int32_t     pos;
//     uint32_t    swap;
//     uint32_t    bit;
// 	curve25519_expand(x1, point);
// 	// print51(x1);
// 	curve25519_copy(x3, x1);
// 	int32_t i;
// 	for (i = 0; i < 32; i++) {
//         t[i] = n[i];
//     }

// 	swap = 0;
// 	int32_t zeroNum;
//     for (pos = 254; pos >= 0; --pos) {
//         bit = t[pos / 8] >> (pos & 7);
//         bit &= 1;
//         swap ^= bit;
//         curve25519_swap_conditional(x2, x3, swap);
//         curve25519_swap_conditional(z2, z3, swap);
//         swap = bit;
// 		curve25519_copy(t0, x2);
// 		curve25519_copy(t1, z2);
// 		curve25519_copy(t2, x3);
// 		curve25519_copy(t3, z3);
//         curve25519_add(x2, x2, z2);
//         curve25519_sub(t0, t0, t1);
// 		curve25519_add(x3, x3, z3);
// 		curve25519_sub(t2, t2, t3);

// 		curve25519_copy(e0, x2);
// 		curve25519_copy(e1, t0);
//         curve25519_square(z2, x2);
//         curve25519_square(t1, t0);
//         curve25519_mul(z3, x3, e1);
//         curve25519_mul(t3, t2, e0);

// 		curve25519_copy(e2, z2);
// 		curve25519_copy(e3, t1);
// 		curve25519_copy(x3, z3);
// 		curve25519_copy(t2, t3);

//         curve25519_mul(x2, z2, t1);
//         curve25519_sub(e3, e2, e3);
//         curve25519_add(t3, z3, t3);
//         curve25519_sub(t2, x3, t2);

// 		curve25519_scalar_product(t1, e3, 121666);
// 		curve25519_square(x3, t3);
//         curve25519_square(e1, t2);

// 		curve25519_add(e2, t1, e2);
//         curve25519_mul(z3, e1, x1);
//         curve25519_mul(z2, e2, e3);
//     }
//     curve25519_swap_conditional(x2, x3, swap);
//     curve25519_swap_conditional(z2, z3, swap);

//     curve25519_recip(z2, z2);
//     curve25519_mul(x2, x2, z2);

//     curve25519_contract(mypublic, x2);
// }

static void x25519_scalar_mul_(U64 out[], U64 scalar[], U64 point[])
{
    __ALIGN64 U64 x1[5], x2[5], x3[5];
    __ALIGN64 U64        z2[5], z3[5];
    __ALIGN64 U64 t0[5], t1[5], t2[5], t3[5], e0[5], e1[5], e2[5], e3[5];

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

        __mb_mask b = _mm512_cmp_epi64_mask(e, get_zero64(), _MM_CMPINT_LT);

        swap = mask_xor (swap, b);
        fe52mb8_cswap(x2, x3, swap);
        fe52mb8_cswap(z2, z3, swap);
        swap = b;
        fe52_copy_mb(t0, x2);
        fe52_copy_mb(t1, z2);
        fe52_copy_mb(t2, x3);
        fe52_copy_mb(t3, z3);
        fe52_add(x2, x2, z2);
        fe52_sub(t0, t0, t1);
		fe52_add(x3, x3, z3);
		fe52_sub(t2, t2, t3);

		fe52_copy_mb(e0, x2);
		fe52_copy_mb(e1, t0);
        fe52_sqr(z2, x2);
        fe52_sqr(t1, t0);
        fe52_mul(z3, x3, e1);
        fe52_mul(t3, t2, e0);

		fe52_copy_mb(e2, z2);
		fe52_copy_mb(e3, t1);
		fe52_copy_mb(x3, z3);
		fe52_copy_mb(t2, t3);

        fe52_mul(x2, z2, t1);
        fe52_sub(e3, e2, e3);
        fe52_add(t3, z3, t3);
        fe52_sub(t2, x3, t2);

		fe52_mul121666(t1, e3);
		fe52_sqr(x3, t3);
        fe52_sqr(e1, t2);

		fe52_add(e2, t1, e2);
        fe52_mul(z3, e1, x1);
        fe52_mul(z2, e2, e3);

        e = slli64(e, 1);
    }

    fe52mb8_cswap(x2, x3, swap);
    fe52mb8_cswap(z2, z3, swap);

    fe52mb8_inv_mod25519(z2, z2);
    fe52_mul(out, x2, z2);
}

__INLINE void x25519_scalar_mul_dual(U64 out[], U64 scalar[], U64 point[])
{
    printf("IF montgomery init\n");
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
    U64 vmask = _mm512_slli_epi64(xor64(e, srli64(e, 1)), 1);
    __mb_mask swap = _mm512_cmp_epi64_mask(vmask, get_zero64(), _MM_CMPINT_LT);

    int bitpos;
    for (bitpos = 0; bitpos >= 0; bitpos--) {
        if (63 == (bitpos % 64)) {
            U64 t = e;
            e = loadu64(&scalar[bitpos/64]);
            vmask = xor64(_mm512_slli_epi64(t, 63), xor64(e, srli64(e, 1)));
            swap  = _mm512_cmp_epi64_mask(vmask, get_zero64(), _MM_CMPINT_LT);
        }
        printf("%0x\n", swap);

        fe52mb8_cswap(x2, x3, swap);
        fe52mb8_cswap(z2, z3, swap);

#if (defined(linux) && ((SIMD_LEN)==512))
        // Avoid reordering optimization by compiler
        U64 Z = get_zero64();
        __asm__ ("vpsllq $1, %0, %0 \n"
             "vpcmpq $1, %2, %0, %1\n"
             : "+x" (vmask), "=k" (swap): "x" (Z) : );
#else
        vmask = _mm512_slli_epi64(vmask, 1);
        swap  = _mm512_cmp_epi64_mask(vmask, get_zero64(), _MM_CMPINT_LT);
#endif

        fe52_sub(tmp0, x3, z3);
        fe52_sub(tmp1, x2, z2);
        fe52_add(x2, x2, z2); 
        fe52_add(z2, x3, z3);

        ed25519_mul_dual_wonorm(z3, z2, x2,tmp0, z2,tmp1);

        ed25519_sqr_dual(tmp0, tmp1, tmp1, x2);

        fe52_add(x3, z3, z2);
        fe52_sub(z2, z3, z2);

        fe52mb8_mul_mod25519_wonorm(x2, tmp1, tmp0);

        fe52_sub(tmp1, tmp1, tmp0);

        ed25519_sqr_dual(z2, x3, z2, x3);

        fe52mb8_mul121666_mod25519_wonorm(z3, tmp1);
        fe52_add(tmp0, tmp0, z3);
        printf("x1 before: \n");
    printU64(x1);
    printf("z2 before: \n");
    printU64(z2);

        ed25519_mul_dual_wonorm(z3, z2, x1,z2, tmp1,tmp0);
    }

    printf("x2 before: \n");
    printU64(x2);
    printf("z2 before: \n");
    printU64(z2);
    printf("x3 before: \n");
    printU64(x3);
    printf("z3 before: \n");
    printU64(z3);

    U64 ta;
    for(int ii = 0; ii < 5; ++ii)
    {
        ta = _mm512_mask_mov_epi64((x2[ii]), (__mmask8)(0xFF), (x3[ii])); 
        (x3[ii]) = _mm512_mask_mov_epi64((x3[ii]), (__mmask8)(0xFF), (x2[ii])); 
        (x2[ii]) = ta;
    }
    for(int ii = 0; ii < 5; ++ii)
    {
        ta = _mm512_mask_mov_epi64((z2[ii]), (__mmask8)(0xFF), (z3[ii])); 
        (z3[ii]) = _mm512_mask_mov_epi64((z3[ii]), (__mmask8)(0xFF), (z2[ii])); 
        (z2[ii]) = ta;
    }

        printf("x2 after: \n");
    printU64(x2);
    printf("z2 after: \n");
    printU64(z2);
    printf("x3 after: \n");
    printU64(x3);
    printf("z3 after: \n");
    printU64(z3);

    // normalize z2 and x2 before inversion
    {
      U64 r0 = z2[0];
      U64 r1 = z2[1];
      U64 r2 = z2[2];
      U64 r3 = z2[3];
      U64 r4 = z2[4];
      NORM_(r, 0,1)
      NORM_(r, 1,2)
      NORM_(r, 2,3)
      NORM_(r, 3,4)
     storeu64(&z2[0], r0);
     storeu64(&z2[1], r1);
     storeu64(&z2[2], r2);
     storeu64(&z2[3], r3);
     storeu64(&z2[4], r4);

      r0 = x2[0];
      r1 = x2[1];
      r2 = x2[2];
      r3 = x2[3];
      r4 = x2[4];
      NORM_(r, 0,1)
      NORM_(r, 1,2)
      NORM_(r, 2,3)
      NORM_(r, 3,4)
     storeu64(&x2[0], r0);
     storeu64(&x2[1], r1);
     storeu64(&x2[2], r2);
     storeu64(&x2[3], r3);
     storeu64(&x2[4], r4);
    }
    // printAVX(x2, 5);
    // printf("\n");
    // printAVX(z2, 5);
        printf("x2\n");
    printAVX(x2, 5);
    printf("z2\n");
    printAVX(z2, 5);
    fe52mb8_inv_mod25519(z2, z2);
    fe52_mul(x2, x2, z2);
    fe52mb8_red_p25519(out, x2);
    printf("out\n");
    printAVX(out, 5);
}

DLL_PUBLIC
// mbx_status MB_FUNC_NAME(mbx_x25519_)(uint8_t pa_shared_key[8][32],
//                        const uint8_t pa_private_key[8][32],
//                        const uint8_t pa_public_key[8][32])
mbx_status MB_FUNC_NAME(mbx_x25519_)(uint8_t* const pa_shared_key[8],
                       const uint8_t* const pa_private_key[8],
                       const uint8_t* const pa_public_key[8])
{
    mbx_status status = 0;
    int buf_no;

    /* test input pointers */
    if(NULL==pa_shared_key || NULL==pa_private_key || NULL==pa_public_key) {
        status = MBX_SET_STS_ALL(MBX_STATUS_NULL_PARAM_ERR);
        return status;
    }

    /* check pointers and values */
    // for(buf_no=0; buf_no<8; buf_no++) {
    //     printf("status before: %0x\n", status);
    //     printf("pa_private_key[%d]: %0x ", buf_no, *pa_private_key[buf_no]);
    //     printf("pa_public_key[%d]: %0x\n", buf_no, *pa_public_key[buf_no]);
    //     uint64_t* shared = (uint64_t*) pa_shared_key[buf_no];
    //     const uint64_t* own_private = (const uint64_t*) pa_private_key[buf_no];
    //     const uint64_t* party_public = (const uint64_t*) pa_public_key[buf_no];
    //     printf("status ing: %0x\n", status);
    //         printf("own_private[%d]: %0lx ", buf_no, *(own_private+1));
    //         printf("party_public[%d]: %0lx\n", buf_no, *party_public);
    //     /* if any of pointer NULL set error status */
    //     if(NULL==own_private || NULL==party_public) {
    //         status = MBX_SET_STS(status, buf_no, MBX_STATUS_NULL_PARAM_ERR);
            
    //         continue;
    //     }
    // }
    // printf("status after: %0x\n", status);

    // for(buf_no=0; buf_no<8; buf_no++) {
    //     uint64_t* shared = (uint64_t*) pa_shared_key[buf_no];
    //     const uint64_t* own_private = (const uint64_t*) pa_private_key[buf_no];
    //     const uint64_t* party_public = (const uint64_t*) pa_public_key[buf_no];

    //     /* if any of pointer NULL set error status */
    //     if(NULL==shared || NULL==own_private || NULL==party_public) {
    //         status = MBX_SET_STS(status, buf_no, MBX_STATUS_NULL_PARAM_ERR);
    //         continue;
    //     }
    // }

    for(buf_no=0; buf_no<8; buf_no++) {
        uint64_t (*shared)[4] = (uint64_t (*)[4])pa_shared_key;
        const uint64_t (*own_private)[4] = (uint64_t (*)[4])pa_private_key;
        const uint64_t (*party_public)[4] = (uint64_t (*)[4])pa_public_key;

        // printf("%0x ", (own_private + buf_no));
        // printf("%0x ", *(own_private + buf_no));
        // printf("%0x ", **(own_private + buf_no));
        // printf("\n");

        /* if any of pointer NULL set error status */
        if(NULL==(shared + buf_no) || NULL==(own_private + buf_no) || NULL==(party_public + buf_no)) {
            status = MBX_SET_STS(status, buf_no, MBX_STATUS_NULL_PARAM_ERR);
            printf("status ing: %0x\n", status);
            continue;
        }
    }

    /* continue processing if there are correct parameters */
    if( MBX_IS_ANY_OK_STS(status) ) {
        __ALIGN64 U64 private64_mb8[4];
        __ALIGN64 U64 pub52_mb8[5];
        __ALIGN64 U64 shared52_mb8[5];

        /* get scalars and convert to MB8 */
        ifma_BNU_transpose_copy((uint64_t (*)[8])private64_mb8, (const uint8_t (*)[32])pa_private_key);
        /* decode keys into scalars according to RFC7748 */
        // private64_mb8[0] = and64_const(private64_mb8[0], 0xfffffffffffffff8);
        // private64_mb8[3] = and64_const(private64_mb8[3], 0x7fffffffffffffff);
        // private64_mb8[3] = or64(private64_mb8[3], set64(0x4000000000000000));

        /* get peer's public keys and convert to MB8 */
        ifma_BNU_to_mb8((uint64_t (*)[8])pub52_mb8, (const uint8_t (*)[32])pa_public_key);

        /* RFC7748: (x25519) ... MUST mask the most significant bit in the final byte.
           This is done to preserve compatibility with point formats .. (compact format??)
        */
        //  pub52_mb8[4] = and64(pub52_mb8[4], loadu64(VPRIME25519_HI));

        /* point multiplication */
        // x25519_scalar_mul_dual(shared52_mb8, private64_mb8, pub52_mb8);
        x25519_scalar_mul(shared52_mb8, private64_mb8, pub52_mb8);

        /* test shared secret before return; all-zero output results when the input is a point of small order. */
        __ALIGN64 U64 shared52_sum = shared52_mb8[0];
        shared52_sum = or64(shared52_sum, shared52_mb8[1]);
        shared52_sum = or64(shared52_sum, shared52_mb8[2]);
        shared52_sum = or64(shared52_sum, shared52_mb8[3]);
        shared52_sum = or64(shared52_sum, shared52_mb8[4]);
        uint8_t stt_mask = cmpeq64_mask(shared52_sum, get_zero64());
        status |= MBX_SET_STS_BY_MASK(status, stt_mask, MBX_STATUS_LOW_ORDER_ERR);

        /* convert result back */
        ifma_mb8_to_BNU((uint64_t (*)[4])pa_shared_key, (const uint64_t (*)[8])shared52_mb8);

        /* clear computed shared keys and it sum */
        MB_FUNC_NAME(zero_)((uint64_t (*)[8])shared52_mb8, sizeof(shared52_mb8)/sizeof(U64));
        MB_FUNC_NAME(zero_)((uint64_t (*)[8])&shared52_sum, 1);

        /* clear copy of the secret keys */
        MB_FUNC_NAME(zero_)((uint64_t (*)[8])private64_mb8, sizeof(private64_mb8)/sizeof(U64));
    }

    return status;
}

static void 
curve25519_scalarmult_montgomery(uint8_t out[8][32], const uint8_t key[32], uint8_t point[8][32]) {
    // for(int k = 0; k < 32; ++k)
	// {
	// 	printf("0x%02x, ",key[k]);
	// }
	// printf("\n");
    /* Convert key of array mod to AVX512 mod first */
    uint64_t key64[4];
    uint32_t i, j;
    for(i = 0; i < 4; ++i)
    {
        key64[i] = 0;
        for(j = 0; j < 8; ++j)
        {
            key64[i] |= ((uint64_t)(key[i * 8 + j])) << (j * 8);
        }
    }

    __ALIGN64 U64 key64_mb8[4] = {_mm512_set1_epi64(key64[0]), 
                                  _mm512_set1_epi64(key64[1]), 
                                  _mm512_set1_epi64(key64[2]), 
                                  _mm512_set1_epi64(key64[3])};

    // printf("scalar: \n");
    // printAVX(key64_mb8, 4);

    /* Main code of montgomery ladder */
    mbx_status status = 0;
    uint32_t buf_no;

    /* Test input pointers */
    if(NULL == key || NULL == key) {
        status = MBX_SET_STS_ALL(MBX_STATUS_NULL_PARAM_ERR);
        // return status;
    }

    for(buf_no=0; buf_no<8; buf_no++) {
        uint64_t (*shared)[4] = (uint64_t (*)[4])out;
        const uint64_t (*party_public)[4] = (uint64_t (*)[4])point;

        /* if any of pointer NULL set error status */
        if(NULL==(shared + buf_no) || NULL==(party_public + buf_no)) {
            status = MBX_SET_STS(status, buf_no, MBX_STATUS_NULL_PARAM_ERR);
            // printf("status ing: %0x\n", status);
            continue;
        }
    }

    /* Start to do montgomery ladder */
    if( MBX_IS_ANY_OK_STS(status) ) {
        __ALIGN64 U64 point52_mb8[5];
        __ALIGN64 U64 out52_mb8[5];

        /* get peer's public keys and convert to MB8 */
        ifma_BNU_to_mb8((uint64_t (*)[8])point52_mb8, (const uint8_t (*)[32])point);

        /* point multiplication */
        // x25519_scalar_mul_dual(out52_mb8, key64_mb8, point52_mb8);
        x25519_scalar_mul_(out52_mb8, key64_mb8, point52_mb8);

        __m512i P[5] = {_mm512_set1_epi64(PRIME25519_LO), 
                      _mm512_set1_epi64(PRIME25519_MID), 
                      _mm512_set1_epi64(PRIME25519_MID), 
                      _mm512_set1_epi64(PRIME25519_MID), 
                      _mm512_set1_epi64(PRIME25519_HI)};

        U64 borrow = get_zero64();
        U64 cr[5];
        for(i = 0; i < 5; ++i)
        {
            borrow = sub64(out52_mb8[i], add64(P[i], _mm512_srli_epi64(borrow, 63)));
            cr[i] = and64_const(borrow, (1ULL << 52) - 1);
        }
    
        U64 underflow_mask = sub64(_mm512_xor_epi64(_mm512_srli_epi64(borrow, 63), _mm512_set1_epi64(1)), _mm512_set1_epi64(1));
    
        borrow = get_zero64();
        for(i = 0; i < 5; ++i)
        {
            borrow = add64(add64(cr[i], _mm512_srli_epi64(borrow, 52)), and64(P[i], underflow_mask));
            out52_mb8[i] = and64_const(borrow, (1ULL << 52) - 1);
        }

        /* test shared secret before return; all-zero output results when the input is a point of small order. */
        __ALIGN64 U64 shared52_sum = cr[0];
        shared52_sum = or64(shared52_sum, out52_mb8[1]);
        shared52_sum = or64(shared52_sum, out52_mb8[2]);
        shared52_sum = or64(shared52_sum, out52_mb8[3]);
        shared52_sum = or64(shared52_sum, out52_mb8[4]);
        uint8_t stt_mask = cmpeq64_mask(shared52_sum, get_zero64());
        status |= MBX_SET_STS_BY_MASK(status, stt_mask, MBX_STATUS_LOW_ORDER_ERR);

        /* convert result back */
        ifma_mb8_to_BNU((uint64_t (*)[4])out, (const uint64_t (*)[8])out52_mb8);

        /* clear computed shared keys and it sum */
        MB_FUNC_NAME(zero_)((uint64_t (*)[8])out52_mb8, sizeof(out52_mb8)/sizeof(U64));
        MB_FUNC_NAME(zero_)((uint64_t (*)[8])&shared52_sum, 1);

        /* clear copy of the secret keys */
        MB_FUNC_NAME(zero_)((uint64_t (*)[8])key64_mb8, sizeof(key64_mb8)/sizeof(U64));
    }

    // return status;
}