#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#if defined(COMPILER_MSVC)
	typedef signed char int8_t;
	typedef unsigned char uint8_t;
	typedef signed short int16_t;
	typedef unsigned short uint16_t;
	typedef signed int int32_t;
	typedef unsigned int uint32_t;
	typedef signed __int64 int64_t;
	typedef unsigned __int64 uint64_t;
#else
	#include <stdint.h>
#endif

#ifdef __GNUC__
#define ASM(a) __asm__(a);
#else
#define ASM(a)
#endif

#define __NOINLINE __attribute__((noinline))
#define __INLINE static __inline__
#define __ALIGN64 __attribute__((aligned(64)))
#define DIGIT_MASK ((uint64_t)0xFFFFFFFFFFFFF)
#define BYTES_REV (1)
#define RADIX_CVT (2)
#define DIGIT_SIZE (52)
#define NUMBER_OF_DIGITS(bitsize,digsize) (((bitsize) + (digsize)-1)/(digsize))
#define PROC_LEN (52)
#define DLL_PUBLIC __attribute__ ((visibility ("default")))
typedef uint32_t mbx_status;
typedef __m512i U64;

__INLINE __mmask8 MB_MASK(int L) {
   return (L > 0) ? (__mmask8)0xFF : (__mmask8)0;
}
__INLINE __mmask64 SB_MASK1(int L, int REV)
{
   if (L <= 0)
      return (__mmask64)0x0;
   if (L > PROC_LEN)
      L = PROC_LEN;
   if (REV)
      return (__mmask64)(0xFFFFFFFFFFFFFFFFULL << ((int)sizeof(__m512i) - L));
   return (__mmask64)(0xFFFFFFFFFFFFFFFFULL >> ((int)sizeof(__m512i) - L));
}

#define TRANSPOSE_8xI64x8(X0_, X1_ ,X2_ ,X3_ ,X4_ ,X5_ ,X6_ ,X7_) {\
            __m512i X01L = _mm512_unpacklo_epi64(X0_, X1_); \
            __m512i X23L = _mm512_unpacklo_epi64(X2_, X3_); \
            __m512i X45L = _mm512_unpacklo_epi64(X4_, X5_); \
            __m512i X67L = _mm512_unpacklo_epi64(X6_, X7_); \
            \
            __m512i X01H = _mm512_unpackhi_epi64(X0_, X1_); \
            __m512i X23H = _mm512_unpackhi_epi64(X2_, X3_); \
            __m512i X45H = _mm512_unpackhi_epi64(X4_, X5_); \
            __m512i X67H = _mm512_unpackhi_epi64(X6_, X7_); \
            \
            __m512i X4567L, X0123L, X4567H, X0123H; \
            X4567L = _mm512_shuffle_i64x2(X45L, X67L, 0b01000100 ); \
            X0_ = _mm512_mask_shuffle_i64x2(X01L, 0b11111100, X23L, X4567L, 0b10000000 ); \
            X2_ = _mm512_mask_shuffle_i64x2(X23L, 0b11110011, X01L, X4567L, 0b11010001 ); \
            \
            X0123L = _mm512_shuffle_i64x2(X01L, X23L, 0b11101110 ); \
            X4_ = _mm512_mask_shuffle_i64x2(X45L, 0b11001111, X0123L, X67L, 0b10001000 ); \
            X6_ = _mm512_mask_shuffle_i64x2(X67L, 0b00111111, X0123L, X45L, 0b10111101 ); \
            \
            X4567H = _mm512_shuffle_i64x2(X45H, X67H, 0b01000100 ); \
            X1_ = _mm512_mask_shuffle_i64x2(X01H, 0b11111100, X23H, X4567H, 0b10000000 ); \
            X3_ = _mm512_mask_shuffle_i64x2(X23H, 0b11110011, X01H, X4567H, 0b11010001 ); \
            \
            X0123H = _mm512_shuffle_i64x2(X01H, X23H, 0b11101110 ); \
            X5_ = _mm512_mask_shuffle_i64x2(X45H, 0b11001111, X0123H, X67H, 0b10001000 ); \
            X7_ = _mm512_mask_shuffle_i64x2(X67H, 0b00111111, X0123H, X45H, 0b10111101 ); \
        }

#if defined(_MSC_VER) && (_MSC_VER < 1920)
  // Disable optimization for VS2017 due to AVX512 masking bug
  #define DISABLE_OPTIMIZATION __pragma(optimize( "", off ))
#else
  #define DISABLE_OPTIMIZATION
#endif

DISABLE_OPTIMIZATION
uint8_t ifma_BNU_transpose_copy(uint64_t out_mb8[4][8], const uint8_t bn[8][32])
{
    // Check input parameters
    // assert(bitLen > 0);

    __mmask8 kbn[8];
    int i;
    for (i = 0; i < 8; ++i)
        kbn[i] = (NULL == bn[i]) ? (__mmask8)0 : (__mmask8)0xFF;

    int len = 4;
    int n = 0;
    // for (n = 0; len > 0; n += 8, out_mb8 += 8) {
        __mmask8 kread = (len >= 8) ? 0xFF : (__mmask8)((1 << len) - 1);

        __m512i X0 = _mm512_maskz_loadu_epi64(kread & kbn[0], bn[0] + n);
        __m512i X1 = _mm512_maskz_loadu_epi64(kread & kbn[1], bn[1] + n);
        __m512i X2 = _mm512_maskz_loadu_epi64(kread & kbn[2], bn[2] + n);
        __m512i X3 = _mm512_maskz_loadu_epi64(kread & kbn[3], bn[3] + n);
        __m512i X4 = _mm512_maskz_loadu_epi64(kread & kbn[4], bn[4] + n);
        __m512i X5 = _mm512_maskz_loadu_epi64(kread & kbn[5], bn[5] + n);
        __m512i X6 = _mm512_maskz_loadu_epi64(kread & kbn[6], bn[6] + n);
        __m512i X7 = _mm512_maskz_loadu_epi64(kread & kbn[7], bn[7] + n);

        TRANSPOSE_8xI64x8(X0, X1, X2, X3, X4, X5, X6, X7);

        _mm512_mask_storeu_epi64(&out_mb8[0], MB_MASK(len--), X0);
        _mm512_mask_storeu_epi64(&out_mb8[1], MB_MASK(len--), X1);
        _mm512_mask_storeu_epi64(&out_mb8[2], MB_MASK(len--), X2);
        _mm512_mask_storeu_epi64(&out_mb8[3], MB_MASK(len--), X3);
        // _mm512_mask_storeu_epi64(&out_mb8[4], MB_MASK(len--), X4);
        // _mm512_mask_storeu_epi64(&out_mb8[5], MB_MASK(len--), X5);
        // _mm512_mask_storeu_epi64(&out_mb8[6], MB_MASK(len--), X6);
        // _mm512_mask_storeu_epi64(&out_mb8[7], MB_MASK(len--), X7);
    // }

    return _mm512_cmpneq_epi64_mask(_mm512_loadu_si512((__m512i*)bn), _mm512_setzero_si512());
    // return 0;
}

/*
// transpose 8 SB into MB including (reverse bytes and) radix 2^64 => 2^52 conversion
//
// covers:
//    - 8 BIGNUM     -> mb8
//    - 8 BNU        -> mb8
//    - 8 hex strings -> mb8
*/
DISABLE_OPTIMIZATION
__INLINE void transform_8sb_to_mb8(U64 out_mb8[], int bitLen, const uint8_t inp[8][32], int inpLen[8], int flag) {
   // inverse bytes (reverse=1)
   const __m512i bswap_mask = _mm512_set_epi64(
                     0x0001020304050607, 0x08090a0b0c0d0e0f,
                     0x1011121314151617, 0x18191a1b1c1d1e1f,
                     0x2021222324252627, 0x28292a2b2c2d2e2f,
                     0x3031323334353637, 0x38393a3b3c3d3e3f);
   // repeat words
   const __m512i idx16 = _mm512_set_epi64(
                     0x0019001800170016,
                     0x0016001500140013,
                     0x0013001200110010,
                     0x0010000f000e000d,
                     0x000c000b000a0009,
                     0x0009000800070006,
                     0x0006000500040003,
                     0x0003000200010000);
   // shift right
   const __m512i shiftR = _mm512_set_epi64(
      12, 8, 4, 0, 12, 8, 4, 0);
   // radix 2^52 mask of digits
   __m512i digMask = _mm512_set1_epi64(DIGIT_MASK);

   int bytesRev = flag & BYTES_REV; /* reverse flag */
   int radixCvt = flag & RADIX_CVT; /* radix (64->52) conversion assumed*/

   int inpBytes = NUMBER_OF_DIGITS(bitLen, 8); /* bytes */
   int outDigits = NUMBER_OF_DIGITS(bitLen, DIGIT_SIZE); /* digits */

   int i;
   for (i = 0; inpBytes > 0; i += PROC_LEN, inpBytes -= PROC_LEN, out_mb8 += 8) {
      int sbidx = bytesRev ? inpBytes - (int)sizeof(__m512i) : i;

      __m512i X0 = _mm512_maskz_loadu_epi8(SB_MASK1(inpLen[0] - i, bytesRev), (__m512i*)&inp[0][sbidx]);
      __m512i X1 = _mm512_maskz_loadu_epi8(SB_MASK1(inpLen[1] - i, bytesRev), (__m512i*)&inp[1][sbidx]);
      __m512i X2 = _mm512_maskz_loadu_epi8(SB_MASK1(inpLen[2] - i, bytesRev), (__m512i*)&inp[2][sbidx]);
      __m512i X3 = _mm512_maskz_loadu_epi8(SB_MASK1(inpLen[3] - i, bytesRev), (__m512i*)&inp[3][sbidx]);
      __m512i X4 = _mm512_maskz_loadu_epi8(SB_MASK1(inpLen[4] - i, bytesRev), (__m512i*)&inp[4][sbidx]);
      __m512i X5 = _mm512_maskz_loadu_epi8(SB_MASK1(inpLen[5] - i, bytesRev), (__m512i*)&inp[5][sbidx]);
      __m512i X6 = _mm512_maskz_loadu_epi8(SB_MASK1(inpLen[6] - i, bytesRev), (__m512i*)&inp[6][sbidx]);
      __m512i X7 = _mm512_maskz_loadu_epi8(SB_MASK1(inpLen[7] - i, bytesRev), (__m512i*)&inp[7][sbidx]);

      if (bytesRev) {
         X0 = _mm512_permutexvar_epi8(bswap_mask, X0);
         X1 = _mm512_permutexvar_epi8(bswap_mask, X1);
         X2 = _mm512_permutexvar_epi8(bswap_mask, X2);
         X3 = _mm512_permutexvar_epi8(bswap_mask, X3);
         X4 = _mm512_permutexvar_epi8(bswap_mask, X4);
         X5 = _mm512_permutexvar_epi8(bswap_mask, X5);
         X6 = _mm512_permutexvar_epi8(bswap_mask, X6);
         X7 = _mm512_permutexvar_epi8(bswap_mask, X7);
      }

      if (radixCvt) {
         X0 = _mm512_permutexvar_epi16(idx16, X0);
         X0 = _mm512_srlv_epi64(X0, shiftR);
         X0 = _mm512_and_si512(X0, digMask); /* probably exceeded instruction */

         X1 = _mm512_permutexvar_epi16(idx16, X1);
         X1 = _mm512_srlv_epi64(X1, shiftR);
         X1 = _mm512_and_si512(X1, digMask);

         X2 = _mm512_permutexvar_epi16(idx16, X2);
         X2 = _mm512_srlv_epi64(X2, shiftR);
         X2 = _mm512_and_si512(X2, digMask);

         X3 = _mm512_permutexvar_epi16(idx16, X3);
         X3 = _mm512_srlv_epi64(X3, shiftR);
         X3 = _mm512_and_si512(X3, digMask);

         X4 = _mm512_permutexvar_epi16(idx16, X4);
         X4 = _mm512_srlv_epi64(X4, shiftR);
         X4 = _mm512_and_si512(X4, digMask);

         X5 = _mm512_permutexvar_epi16(idx16, X5);
         X5 = _mm512_srlv_epi64(X5, shiftR);
         X5 = _mm512_and_si512(X5, digMask);

         X6 = _mm512_permutexvar_epi16(idx16, X6);
         X6 = _mm512_srlv_epi64(X6, shiftR);
         X6 = _mm512_and_si512(X6, digMask);

         X7 = _mm512_permutexvar_epi16(idx16, X7);
         X7 = _mm512_srlv_epi64(X7, shiftR);
         X7 = _mm512_and_si512(X7, digMask);
      }

      // transpose 8 digits at a time
      TRANSPOSE_8xI64x8(X0, X1, X2, X3, X4, X5, X6, X7);

      // store transposed digits
      _mm512_mask_storeu_epi64(&out_mb8[0], MB_MASK(outDigits--), X0);
      _mm512_mask_storeu_epi64(&out_mb8[1], MB_MASK(outDigits--), X1);
      _mm512_mask_storeu_epi64(&out_mb8[2], MB_MASK(outDigits--), X2);
      _mm512_mask_storeu_epi64(&out_mb8[3], MB_MASK(outDigits--), X3);
      _mm512_mask_storeu_epi64(&out_mb8[4], MB_MASK(outDigits--), X4);
      _mm512_mask_storeu_epi64(&out_mb8[5], MB_MASK(outDigits--), X5);
      _mm512_mask_storeu_epi64(&out_mb8[6], MB_MASK(outDigits--), X6);
      _mm512_mask_storeu_epi64(&out_mb8[7], MB_MASK(outDigits--), X7);
   }
}

uint8_t ifma_BNU_to_mb8(uint64_t out_mb8[5][8], const uint8_t bn[8][32])
{
   // Check input parameters
//    assert(bitLen > 0);

   int byteLens[8];
   int byteLen = 32;
   int i;
   for (i = 0; i < 8; ++i)
       byteLens[i] = (NULL != bn[i]) ? byteLen : 0;

   transform_8sb_to_mb8((U64*)out_mb8, 256, bn, byteLens, 2);

   return _mm512_cmpneq_epi64_mask(_mm512_loadu_si512((__m512i*)bn), _mm512_setzero_si512());
}

/*
// transpose MB into 8 SB including (reverse bytes and) radix 2^52 => 2^64 conversion
//
// covers:
//    - mb8 -> 8 BNU
//    - mb8 -> 8 hex strings
*/
DISABLE_OPTIMIZATION
__INLINE void transform_mb8_to_8sb(uint64_t out[8][4], int outLen[8], const U64 inp_mb8[5], int bitLen, int flag)
{
   // inverse bytes (reverse=1)
   const __m512i bswap_mask = _mm512_set_epi64(
                     0x0001020304050607, 0x08090a0b0c0d0e0f,
                     0x1011121314151617, 0x18191a1b1c1d1e1f,
                     0x2021222324252627, 0x28292a2b2c2d2e2f,
                     0x3031323334353637, 0x38393a3b3c3d3e3f);

   const __m512i shiftL = _mm512_set_epi64(4, 0, 4, 0, 4, 0, 4, 0);

   const __m512i permutation1 = _mm512_set_epi64(0x3f3f3f3f3f3f3f3f, // {63,63,63,63,63,63,63,63}
                                                0x3f3f3f3f3e3d3c3b,  // {63,63,63,63,62,61,60,59}
                                                0x3737363534333231,  // {55,55,54,53,52,51,50,49}
                                                0x302e2d2c2b2a2928,  // {48,46,45,44,43,42,41,40}
                                                0x1f1f1f1f1f1f1e1d,  // {31,31,31,31,31,31,30,29}
                                                0x1717171716151413,  // {23,23,23,23,22,21,20,19}
                                                0x0f0f0f0e0d0c0b0a,  // {15,15,15,14,13,12,11,10}
                                                0x0706050403020100); // { 7, 6, 5, 4, 3, 2, 1, 0}

   const __m512i permutation2 = _mm512_set_epi64(0x3f3f3f3f3f3f3f3f, // {63,63,63,63,63,63,63,63}
                                                0x3f3f3f3f3f3f3f3f,  // {63,63,63,63,63,63,63,63}
                                                0x3a39383737373737,  // {58,57,56,55,55,55,55,55}
                                                0x2727272727272726,  // {39,39,39,39,39,39,39,38}
                                                0x2524232221201f1f,  // {37,36,35,34,33,32,31,31}
                                                0x1c1b1a1918171717,  // {28,27,26,25,24,23,23,23}
                                                0x1211100f0f0f0f0f,  // {18,17,16,15,15,15,15,15}
                                                0x0908070707070707); // { 9, 8, 7, 7, 7, 7, 7, 7}
   int bytesRev = flag & BYTES_REV; /* reverse flag */
   int radixCvt = flag & RADIX_CVT; /* radix (52->64) conversion assumed */

   int inpDigits = NUMBER_OF_DIGITS(bitLen, DIGIT_SIZE); /* digits */
   int outBytes = NUMBER_OF_DIGITS(bitLen, 8); /* bytes */

   int i;
   for (i = 0; outBytes > 0; i += PROC_LEN, outBytes -= PROC_LEN, inp_mb8 += 8) {
      int sbidx = bytesRev ? outBytes - (int)sizeof(__m512i) : i;

      __m512i X0 = _mm512_maskz_loadu_epi64(MB_MASK(inpDigits--), &inp_mb8[0]);
      __m512i X1 = _mm512_maskz_loadu_epi64(MB_MASK(inpDigits--), &inp_mb8[1]);
      __m512i X2 = _mm512_maskz_loadu_epi64(MB_MASK(inpDigits--), &inp_mb8[2]);
      __m512i X3 = _mm512_maskz_loadu_epi64(MB_MASK(inpDigits--), &inp_mb8[3]);
      __m512i X4 = _mm512_maskz_loadu_epi64(MB_MASK(inpDigits--), &inp_mb8[4]);
    //   __m512i X5 = _mm512_maskz_loadu_epi64(MB_MASK(inpDigits--), &inp_mb8[5]);
    //   __m512i X6 = _mm512_maskz_loadu_epi64(MB_MASK(inpDigits--), &inp_mb8[6]);
    //   __m512i X7 = _mm512_maskz_loadu_epi64(MB_MASK(inpDigits--), &inp_mb8[7]);
      __m512i X5 = _mm512_set1_epi64(0);
      __m512i X6 = _mm512_set1_epi64(0);
      __m512i X7 = _mm512_set1_epi64(0);

      // transpose 8 digits at a time
      TRANSPOSE_8xI64x8(X0, X1, X2, X3, X4, X5, X6, X7);

      if (radixCvt) {
         __m512i T;
         X0 = _mm512_sllv_epi64(X0, shiftL);
         T = _mm512_permutexvar_epi8(permutation1, X0);
         X0 = _mm512_permutexvar_epi8(permutation2, X0);
         X0 = _mm512_or_si512(X0, T);

         X1 = _mm512_sllv_epi64(X1, shiftL);
         T = _mm512_permutexvar_epi8(permutation1, X1);
         X1 = _mm512_permutexvar_epi8(permutation2, X1);
         X1 = _mm512_or_si512(X1, T);

         X2 = _mm512_sllv_epi64(X2, shiftL);
         T = _mm512_permutexvar_epi8(permutation1, X2);
         X2 = _mm512_permutexvar_epi8(permutation2, X2);
         X2 = _mm512_or_si512(X2, T);

         X3 = _mm512_sllv_epi64(X3, shiftL);
         T = _mm512_permutexvar_epi8(permutation1, X3);
         X3 = _mm512_permutexvar_epi8(permutation2, X3);
         X3 = _mm512_or_si512(X3, T);

         X4 = _mm512_sllv_epi64(X4, shiftL);
         T = _mm512_permutexvar_epi8(permutation1, X4);
         X4 = _mm512_permutexvar_epi8(permutation2, X4);
         X4 = _mm512_or_si512(X4, T);

         X5 = _mm512_sllv_epi64(X5, shiftL);
         T = _mm512_permutexvar_epi8(permutation1, X5);
         X5 = _mm512_permutexvar_epi8(permutation2, X5);
         X5 = _mm512_or_si512(X5, T);

         X6 = _mm512_sllv_epi64(X6, shiftL);
         T = _mm512_permutexvar_epi8(permutation1, X6);
         X6 = _mm512_permutexvar_epi8(permutation2, X6);
         X6 = _mm512_or_si512(X6, T);

         X7 = _mm512_sllv_epi64(X7, shiftL);
         T = _mm512_permutexvar_epi8(permutation1, X7);
         X7 = _mm512_permutexvar_epi8(permutation2, X7);
         X7 = _mm512_or_si512(X7, T);
      }

      if (bytesRev) {
         X0 = _mm512_permutexvar_epi8(bswap_mask, X0);
         X1 = _mm512_permutexvar_epi8(bswap_mask, X1);
         X2 = _mm512_permutexvar_epi8(bswap_mask, X2);
         X3 = _mm512_permutexvar_epi8(bswap_mask, X3);
         X4 = _mm512_permutexvar_epi8(bswap_mask, X4);
         X5 = _mm512_permutexvar_epi8(bswap_mask, X5);
         X6 = _mm512_permutexvar_epi8(bswap_mask, X6);
         X7 = _mm512_permutexvar_epi8(bswap_mask, X7);
      }

      // store transposed digits
      _mm512_mask_storeu_epi8(out[0] + sbidx, SB_MASK1(outLen[0] - i, bytesRev), X0);
      _mm512_mask_storeu_epi8(out[1] + sbidx, SB_MASK1(outLen[1] - i, bytesRev), X1);
      _mm512_mask_storeu_epi8(out[2] + sbidx, SB_MASK1(outLen[2] - i, bytesRev), X2);
      _mm512_mask_storeu_epi8(out[3] + sbidx, SB_MASK1(outLen[3] - i, bytesRev), X3);
      _mm512_mask_storeu_epi8(out[4] + sbidx, SB_MASK1(outLen[4] - i, bytesRev), X4);
      _mm512_mask_storeu_epi8(out[5] + sbidx, SB_MASK1(outLen[5] - i, bytesRev), X5);
      _mm512_mask_storeu_epi8(out[6] + sbidx, SB_MASK1(outLen[6] - i, bytesRev), X6);
      _mm512_mask_storeu_epi8(out[7] + sbidx, SB_MASK1(outLen[7] - i, bytesRev), X7);
   }

}

uint8_t ifma_mb8_to_BNU(uint64_t out_bn[8][4], const uint64_t inp_mb8[5][8])
{
    // Check input parameters
    // assert(bitLen > 0);

    // int bnu_bitlen = 256; // gres: output length is multiple 64
    int byteLens[8];
    int i;
    for (i = 0; i < 8; ++i)
        //gres: byteLens[i] = (NULL != out_bn[i]) ? NUMBER_OF_DIGITS(bitLen, 8) : 0;
         byteLens[i] = (NULL != out_bn[i]) ? 32 : 0;

    transform_mb8_to_8sb(out_bn, byteLens, (U64*)inp_mb8, 256, RADIX_CVT);

   return _mm512_cmpneq_epi64_mask(_mm512_loadu_si512((__m512i*)out_bn), _mm512_setzero_si512());
}

typedef __mmask8 __mb_mask;

        #define SIMD_LEN   512
        #define SIMD_BYTES (SIMD_LEN/8)
        #define MB_WIDTH   (SIMD_LEN/64)

        __INLINE U64 loadu64(const void *p) {
            return _mm512_loadu_si512((U64*)p);
        }

        __INLINE U64 loadstream64(const void *p) {
            return _mm512_stream_load_si512 ((U64*)p);
        }

        __INLINE void storeu64(const void *p, U64 v) {
            _mm512_storeu_si512((U64*)p, v);
        }

        #define mask_mov64 _mm512_mask_mov_epi64
        #define set64      _mm512_set1_epi64

        __INLINE U64 fma52lo(U64 a, U64 b, U64 c) {
            return _mm512_madd52lo_epu64(a, b, c);
        }

        __INLINE U64 fma52hi(U64 a, U64 b, U64 c) {
            return _mm512_madd52hi_epu64(a, b, c);
        }

        __INLINE U64 mul52lo(U64 b, U64 c) {
            return _mm512_madd52lo_epu64(_mm512_setzero_si512(), b, c);
        }

        #ifdef __GNUC__
            // memory ops intrinsics - force load from original buffer
            #define _mm512_madd52lo_epu64_(r, a, b, c, o) {\
                r=a; \
                __asm__ ( "vpmadd52luq " #o "(%2), %1, %0" : "+x" (r): "x" (b), "r" (c) ); \
            }

            #define _mm512_madd52hi_epu64_(r, a, b, c, o) {\
                r=a; \
                __asm__ ( "vpmadd52huq " #o "(%2), %1, %0" : "+x" (r): "x" (b), "r" (c) ); \
            }

            __INLINE U64 select64(__mb_mask k, U64 v, U64 *d) {
                __asm__("vmovdqu64 %2, %%zmm0 \n"
                        "vpblendmq %%zmm0, %0, %0 %{%1%} \n"
                : "+v"(v)
                : "Yk"(k), "m"(*d)
                : "zmm0");
                return v;
            }
            
        #else
            // Use IFMA instrinsics for all other compilers
            #define _mm512_madd52lo_epu64_(r, a, b, c, o) {\
                r=fma52lo(a, b, _mm512_loadu_si512((U64*)(((char*)c)+o))); \
            }

            #define _mm512_madd52hi_epu64_(r, a, b, c, o) {\
                r=fma52hi(a, b, _mm512_loadu_si512((U64*)(((char*)c)+o))); \
            }

            #pragma optimize("", off)
            __INLINE U64 select64(__mb_mask k, U64 v, U64 *d) {
                return _mm512_mask_blend_epi64(k, v, _mm512_load_si512(d));
            }
            
            #pragma optimize("", on)
        #endif

        #define fma52lo_mem(r, a, b, c, o) _mm512_madd52lo_epu64_(r, a, b, c, o) // gres
        #define fma52hi_mem(r, a, b, c, o) _mm512_madd52hi_epu64_(r, a, b, c, o) // gres

        __INLINE U64 add64(U64 a, U64 b) {
            return _mm512_add_epi64(a, b);
        }

        __INLINE U64 sub64(U64 a, U64 b) {
            return _mm512_sub_epi64(a, b);
        }

        __INLINE U64 get_zero64() {
            return _mm512_setzero_si512();
        }

        __INLINE void set_zero64(U64 *a) {
            *a = _mm512_xor_si512(*a, *a);
        }

        __INLINE U64 set1(unsigned long long a) {
            return _mm512_set1_epi64((long long)a);
        }

        __INLINE U64 srli64(U64 a, int s) {
            return _mm512_srli_epi64(a, s);
        }

        #define srai64 _mm512_srai_epi64
        #define slli64 _mm512_slli_epi64

        __INLINE U64 and64_const(U64 a, unsigned long long mask) {
            return _mm512_and_epi64(a, _mm512_set1_epi64((long long)mask));
        }

        __INLINE U64 and64(U64 a, U64 mask) {
            return _mm512_and_epi64(a, mask);
        }

        #define or64         _mm512_or_epi64
        #define xor64        _mm512_xor_epi64
        #define cmp64_mask   _mm512_cmp_epi64_mask
        #define cmpeq16_mask _mm512_cmpeq_epi16_mask
        #define cmpeq64_mask _mm512_cmpeq_epi64_mask

        // Mask operations
        #define mask_blend64 _mm512_mask_blend_epi64
        #define mask_add64   _mm512_mask_add_epi64
        #define mask_sub64   _mm512_mask_sub_epi64
        #define maskz_sub64  _mm512_maskz_sub_epi64

        __INLINE __mb_mask is_zero(U64* p, int len) {
            U64 Z = p[0];
            for(int i = 1; i < len; i++) {
                Z = or64(Z, p[i]);
            }

            return cmpeq64_mask(Z, get_zero64());
        }

        #if defined(_MSC_VER) && !defined(__INTEL_COMPILER) && !defined(__INTEL_LLVM_COMPILER) // for MSVC
            #define mask_xor(m1,m2) (__mb_mask)(_mm512_kxor((m1),(m2)))
        #else
            #define mask_xor _kxor_mask8
        #endif
        
        #define get_mask(a)       (a)
        #define get_mask_value(a) (a)

#define MB_FUNC_NAME(name) name ## mb8

#define MBX_STATUS_OK                 (0)
#define MBX_STATUS_MISMATCH_PARAM_ERR (1)
#define MBX_STATUS_NULL_PARAM_ERR     (2)
#define MBX_STATUS_LOW_ORDER_ERR      (4)
#define MBX_STATUS_SIGNATURE_ERR      (8)

__INLINE mbx_status MBX_SET_STS(mbx_status status, int numb, mbx_status sttVal)
{
   numb &= 7; /* 0 <= numb < 8 */
   status &= (mbx_status)(~(0xF << (numb*4)));
   return status |= (sttVal & 0xF) << (numb*4);
}

__INLINE mbx_status MBX_GET_STS(mbx_status status, int numb)
{
   return (status >>(numb*4)) & 0xF;
}
__INLINE mbx_status MBX_SET_STS_ALL(mbx_status stsVal)
{
   return (stsVal<<4*7) | (stsVal<<4*6) | (stsVal<<4*5) | (stsVal<<4*4)  | (stsVal<<4*3) | (stsVal<<4*2) | (stsVal<<4*1) | stsVal;
}

__INLINE mbx_status MBX_SET_STS_BY_MASK(mbx_status status, uint8_t mask, mbx_status sttVal)
{
   int numb;

   for(numb=0; numb<8; numb++) {
      mbx_status buf_stt = (0 - ((mask>>numb) &1)) & sttVal;
      status = MBX_SET_STS(status, numb, buf_stt);
   }
   return status;
}

__INLINE int MBX_IS_ANY_OK_STS(mbx_status status)
{
   int ret = MBX_STATUS_OK==MBX_GET_STS(status, 0)
          || MBX_STATUS_OK==MBX_GET_STS(status, 1)
          || MBX_STATUS_OK==MBX_GET_STS(status, 2)
          || MBX_STATUS_OK==MBX_GET_STS(status, 3)
          || MBX_STATUS_OK==MBX_GET_STS(status, 4)
          || MBX_STATUS_OK==MBX_GET_STS(status, 5)
          || MBX_STATUS_OK==MBX_GET_STS(status, 6)
          || MBX_STATUS_OK==MBX_GET_STS(status, 7);
   return ret;
}

__NOINLINE
void zero_mb8(uint64_t (*out)[8], int len)
{
#if defined(__GNUC__)
   // Avoid dead code elimination for GNU compilers
   ASM("");
#endif
   __m512i T = _mm512_setzero_si512();
   int i;
   for(i=0; i<len; i++)
      _mm512_storeu_si512(out[i], T);
}

/* repetitions */
#define  REP2_DECL(a)   a, a
#define  REP4_DECL(a)   REP2_DECL(a), REP2_DECL(a)
#define  REP8_DECL(a)   REP4_DECL(a), REP4_DECL(a)

void printAVX(__m512i input[], uint32_t k)
{
    __m128i result;
    for (uint32_t i = 0; i < k; i++) {

        result = _mm512_extracti64x2_epi64(input[i], 0);
        printf("%0llx ", result[0]);
        printf("%0llx ", result[1]);
        result = _mm512_extracti64x2_epi64(input[i], 1);
        printf("%0llx ", result[0]);
        printf("%0llx ", result[1]);
        result = _mm512_extracti64x2_epi64(input[i], 2);
        printf("%0llx ", result[0]);
        printf("%0llx ", result[1]);
        result = _mm512_extracti64x2_epi64(input[i], 3);
        printf("%0llx ", result[0]);
        printf("%0llx ", result[1]);
        printf("\n");
    }
}

void printavxone(U64 input)
{
    __m128i result;
        result = _mm512_extracti64x2_epi64(input, 0);
        printf("%0llx ", result[0]);
        printf("%0llx ", result[1]);
        result = _mm512_extracti64x2_epi64(input, 1);
        printf("%0llx ", result[0]);
        printf("%0llx ", result[1]);
        result = _mm512_extracti64x2_epi64(input, 2);
        printf("%0llx ", result[0]);
        printf("%0llx ", result[1]);
        result = _mm512_extracti64x2_epi64(input, 3);
        printf("%0llx ", result[0]);
        printf("%0llx ", result[1]);
        printf("\n");
}