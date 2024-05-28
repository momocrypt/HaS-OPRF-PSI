// #include "ifma_arith_ed25519.h"

// /*
//  * Usage: scalar mult basepoint with comb algorithm
//  * Author: Mengqing Yang
//  * Time: 2023/11/10 17:16 +8
//  * Mail: yangmengqing1\@jd.com
//  */
// void 
// ge25519_scalarmult_base_comb(ge25519 *h, const uint8_t *a)
// {
//     uint8_t         e[256];
//     ge25519_p1p1    r;
//     ge25519         s;
//     ge25519_niels   t;
//     int32_t         i;
//     int32_t         u;
//     int32_t         j;
//     int32_t         k;

//     for (i = 0; i < 32; ++i) {
//         e[8 * i] = (a[i] >> 0) & 1;
//         e[8 * i + 1] = (a[i] >> 1) & 1;
//         e[8 * i + 2] = (a[i] >> 2) & 1;
//         e[8 * i + 3] = (a[i] >> 3) & 1;
//         e[8 * i + 4] = (a[i] >> 4) & 1;
//         e[8 * i + 5] = (a[i] >> 5) & 1;
//         e[8 * i + 6] = (a[i] >> 6) & 1;
//         e[8 * i + 7] = (a[i] >> 7) & 1;
//     }
//     ge25519_p3_0(h);

//         for (j = 0; j < 8; ++j) {
//             u = 0;
//             for (k = 0; k < 8; ++k) {
//                 u += (1 << k) * e[32 * j + 4 * k + 3];
//             }
//             ge25519_comb_base(&t, j, u);
// 			ge25519_nielsadd2(h, &t);
//         }
//         ge25519_double(h, h);
//         for (j = 0; j < 8; ++j) {
//             u = 0;
//             for (k = 0; k < 8; ++k) {
//                 u += (1 << k) * e[32 * j + 4 * k + 2];
//             }
//             ge25519_comb_base(&t, j, u);
//             ge25519_nielsadd2(h, &t);
//         }
// 		ge25519_double(h, h);
//         for (j = 0; j < 8; ++j) {
//             u = 0;
//             for (k = 0; k < 8; ++k) {
//                 u += (1 << k) * e[32 * j + 4 * k + 1];
//             }
//             ge25519_comb_base(&t, j, u);
//             ge25519_nielsadd2(h, &t);
//         }
//         ge25519_double(h, h);
//         for (j = 0; j < 8; ++j) {
//             u = 0;
//             for (k = 0; k < 8; ++k) {
//                 u += (1 << k) * e[32 * j + 4 * k];
//             }
//             ge25519_comb_base(&t, j, u);
//             ge25519_nielsadd2(h, &t);
//         }

// }

/*
// r = p + q
*/
static void ge52_add_precomp(ge52_p1p1_mb *r, const ge52_ext_mb *p, const ge52_precomp_mb *q)
{
   U64 t0[5];

   fe52_add(r->X, p->Y, p->X);      // X3 = Y1+X1
   fe52_sub(r->Y, p->Y, p->X);      // Y3 = Y1-X1
   fe52_mul(r->Z, r->X, q->yaddx);  // Z3 = X3*yplusx2
   fe52_mul(r->Y, r->Y, q->ysubx);  // Y3 = Y3*yminisx2
   fe52_mul(r->T, q->t2d, p->T);    // T3 = T1*xy2d2
   fe52_add(t0, p->Z, p->Z);        // t0 = Z1+Z1
   fe52_sub(r->X, r->Z, r->Y);      // X3 = Z3-Y3 = X3*yplusx2 - Y3*yminisx2 = (Y1+X1)*yplusx2 - (Y1-X1)*yminisx2
   fe52_add(r->Y, r->Z, r->Y);      // Y3 = Z3+Y3 = X3*yplusx2 + Y3*yminisx2 = (Y1+X1)*yplusx2 + (Y1-X1)*yminisx2
   fe52_add(r->Z, t0, r->T);        // Z3 = 2*Z1 + T1*xy2d2
   fe52_sub(r->T, t0, r->T);        // T3 = 2*Z1 - T1*xy2d2
}

/* r = 2 * p */
static void ge_dbl(ge52_p1p1_mb *r, const ge52_homo_mb* p)
{
   fe52_mb t0;

   fe52_sqr(r->X, p->X);
   fe52_sqr(r->Z, p->Y);
   fe52_sqr(r->T, p->Z);
   fe52_add(r->T, r->T, r->T);
   fe52_add(r->Y, p->X, p->Y);
   fe52_sqr(t0, r->Y);
   fe52_add(r->Y, r->Z, r->X);
   fe52_sub(r->Z, r->Z, r->X);
   fe52_sub(r->X, t0, r->Y);
   fe52_sub(r->T, r->T, r->Z);
}

static void
ge52_double(ge52_ext_mb *r, const ge52_ext_mb *p) {
	ge52_p1p1_mb t;
    ge52_homo_mb q;
	ge52_ext_to_homo_mb(&q, p);
    ge_dbl(&t, &q);
	ge52_p1p1_to_ext_mb(r, &t);
}

// #include <stdint.h>
// #include <stdio.h>
// #include <stdlib.h>
// #include <stddef.h>
// #include <string.h>
// typedef __uint64_t uint64_t;


// static void
// curve25519_copy_mb(bignum25519 out, const bignum25519 in) {
// 	out[0] = in[0];
// 	out[1] = in[1];
// 	out[2] = in[2];
// 	out[3] = in[3];
// 	out[4] = in[4];
// }


static void ge25519_comb_base_mb(ge52_precomp_mb* lo8, ge52_precomp_mb* hi8,
                                const uint32_t j, U64 idx)
{
    /* select p0, p1 wrt idx0, idx1 indexes */
    // static const ge25519_niels tabl[8][256] = {
    //     #include "curve25519-yx-basecombtable52.h"
    // };

    uint32_t index[16];
    __m128i mm;

    mm = _mm512_extracti32x4_epi32(idx, 0);
    index[0] = _mm_extract_epi32(mm, 0);
    // mm = _mm512_extracti32x4_epi32(idx, 0);
    index[1] = _mm_extract_epi32(mm, 1);
    // mm = _mm512_extracti32x4_epi32(idx, 0);
    index[2] = _mm_extract_epi32(mm, 2);
    // mm = _mm512_extracti32x4_epi32(idx, 0);
    index[3] = _mm_extract_epi32(mm, 3);

    mm = _mm512_extracti32x4_epi32(idx, 1);
    index[4] = _mm_extract_epi32(mm, 0);
    // mm = _mm512_extracti32x4_epi32(idx, 1);
    index[5] = _mm_extract_epi32(mm, 1);
    // mm = _mm512_extracti32x4_epi32(idx, 1);
    index[6] = _mm_extract_epi32(mm, 2);
    // mm = _mm512_extracti32x4_epi32(idx, 1);
    index[7] = _mm_extract_epi32(mm, 3);

    mm = _mm512_extracti32x4_epi32(idx, 2);
    index[8] = _mm_extract_epi32(mm, 0);
    // mm = _mm512_extracti32x4_epi32(idx, 2);
    index[9] = _mm_extract_epi32(mm, 1);
    // mm = _mm512_extracti32x4_epi32(idx, 2);
    index[10] = _mm_extract_epi32(mm, 2);
    // mm = _mm512_extracti32x4_epi32(idx, 2);
    index[11] = _mm_extract_epi32(mm, 3);

    mm = _mm512_extracti32x4_epi32(idx, 3);
    index[12] = _mm_extract_epi32(mm, 0);
    // mm = _mm512_extracti32x4_epi32(idx, 3);
    index[13] = _mm_extract_epi32(mm, 1);
    // mm = _mm512_extracti32x4_epi32(idx, 3);
    index[14] = _mm_extract_epi32(mm, 2);
    // mm = _mm512_extracti32x4_epi32(idx, 3);
    index[15] = _mm_extract_epi32(mm, 3);

    lo8->ysubx[0] = _mm512_set_epi64(tabl[2 * j][index[14]].ysubx[0], tabl[2 * j][index[12]].ysubx[0], 
                                     tabl[2 * j][index[10]].ysubx[0],  tabl[2 * j][index[8]].ysubx[0], 
                                      tabl[2 * j][index[6]].ysubx[0],  tabl[2 * j][index[4]].ysubx[0], 
                                      tabl[2 * j][index[2]].ysubx[0],  tabl[2 * j][index[0]].ysubx[0]);
    lo8->ysubx[1] = _mm512_set_epi64(tabl[2 * j][index[14]].ysubx[1], tabl[2 * j][index[12]].ysubx[1], 
                                     tabl[2 * j][index[10]].ysubx[1],  tabl[2 * j][index[8]].ysubx[1], 
                                      tabl[2 * j][index[6]].ysubx[1],  tabl[2 * j][index[4]].ysubx[1], 
                                      tabl[2 * j][index[2]].ysubx[1],  tabl[2 * j][index[0]].ysubx[1]);
    lo8->ysubx[2] = _mm512_set_epi64(tabl[2 * j][index[14]].ysubx[2], tabl[2 * j][index[12]].ysubx[2], 
                                     tabl[2 * j][index[10]].ysubx[2],  tabl[2 * j][index[8]].ysubx[2], 
                                      tabl[2 * j][index[6]].ysubx[2],  tabl[2 * j][index[4]].ysubx[2], 
                                      tabl[2 * j][index[2]].ysubx[2],  tabl[2 * j][index[0]].ysubx[2]);                                  
    lo8->ysubx[3] = _mm512_set_epi64(tabl[2 * j][index[14]].ysubx[3], tabl[2 * j][index[12]].ysubx[3], 
                                     tabl[2 * j][index[10]].ysubx[3],  tabl[2 * j][index[8]].ysubx[3], 
                                      tabl[2 * j][index[6]].ysubx[3],  tabl[2 * j][index[4]].ysubx[3], 
                                      tabl[2 * j][index[2]].ysubx[3],  tabl[2 * j][index[0]].ysubx[3]);
    lo8->ysubx[4] = _mm512_set_epi64(tabl[2 * j][index[14]].ysubx[4], tabl[2 * j][index[12]].ysubx[4], 
                                     tabl[2 * j][index[10]].ysubx[4],  tabl[2 * j][index[8]].ysubx[4], 
                                      tabl[2 * j][index[6]].ysubx[4],  tabl[2 * j][index[4]].ysubx[4], 
                                      tabl[2 * j][index[2]].ysubx[4],  tabl[2 * j][index[0]].ysubx[4]);
    
    lo8->yaddx[0] = _mm512_set_epi64(tabl[2 * j][index[14]].yaddx[0],  tabl[2 * j][index[12]].yaddx[0], 
                                     tabl[2 * j][index[10]].yaddx[0],   tabl[2 * j][index[8]].yaddx[0], 
                                      tabl[2 * j][index[6]].yaddx[0],   tabl[2 * j][index[4]].yaddx[0], 
                                      tabl[2 * j][index[2]].yaddx[0],   tabl[2 * j][index[0]].yaddx[0]);
    lo8->yaddx[1] = _mm512_set_epi64(tabl[2 * j][index[14]].yaddx[1],  tabl[2 * j][index[12]].yaddx[1], 
                                     tabl[2 * j][index[10]].yaddx[1],   tabl[2 * j][index[8]].yaddx[1], 
                                      tabl[2 * j][index[6]].yaddx[1],   tabl[2 * j][index[4]].yaddx[1], 
                                      tabl[2 * j][index[2]].yaddx[1],   tabl[2 * j][index[0]].yaddx[1]);
    lo8->yaddx[2] = _mm512_set_epi64(tabl[2 * j][index[14]].yaddx[2],  tabl[2 * j][index[12]].yaddx[2], 
                                     tabl[2 * j][index[10]].yaddx[2],   tabl[2 * j][index[8]].yaddx[2], 
                                      tabl[2 * j][index[6]].yaddx[2],   tabl[2 * j][index[4]].yaddx[2], 
                                      tabl[2 * j][index[2]].yaddx[2],   tabl[2 * j][index[0]].yaddx[2]);
    lo8->yaddx[3] = _mm512_set_epi64(tabl[2 * j][index[14]].yaddx[3],  tabl[2 * j][index[12]].yaddx[3], 
                                     tabl[2 * j][index[10]].yaddx[3],   tabl[2 * j][index[8]].yaddx[3], 
                                      tabl[2 * j][index[6]].yaddx[3],   tabl[2 * j][index[4]].yaddx[3], 
                                      tabl[2 * j][index[2]].yaddx[3],   tabl[2 * j][index[0]].yaddx[3]);
    lo8->yaddx[4] = _mm512_set_epi64(tabl[2 * j][index[14]].yaddx[4],  tabl[2 * j][index[12]].yaddx[4], 
                                     tabl[2 * j][index[10]].yaddx[4],   tabl[2 * j][index[8]].yaddx[4], 
                                      tabl[2 * j][index[6]].yaddx[4],   tabl[2 * j][index[4]].yaddx[4], 
                                      tabl[2 * j][index[2]].yaddx[4],   tabl[2 * j][index[0]].yaddx[4]);

    lo8->t2d[0]   = _mm512_set_epi64(tabl[2 * j][index[14]].t2d[0], tabl[2 * j][index[12]].t2d[0], 
                                     tabl[2 * j][index[10]].t2d[0],  tabl[2 * j][index[8]].t2d[0], 
                                      tabl[2 * j][index[6]].t2d[0],  tabl[2 * j][index[4]].t2d[0], 
                                      tabl[2 * j][index[2]].t2d[0],  tabl[2 * j][index[0]].t2d[0]);
    lo8->t2d[1]   = _mm512_set_epi64(tabl[2 * j][index[14]].t2d[1], tabl[2 * j][index[12]].t2d[1], 
                                     tabl[2 * j][index[10]].t2d[1],  tabl[2 * j][index[8]].t2d[1], 
                                      tabl[2 * j][index[6]].t2d[1],  tabl[2 * j][index[4]].t2d[1], 
                                      tabl[2 * j][index[2]].t2d[1],  tabl[2 * j][index[0]].t2d[1]);
    lo8->t2d[2]   = _mm512_set_epi64(tabl[2 * j][index[14]].t2d[2], tabl[2 * j][index[12]].t2d[2], 
                                     tabl[2 * j][index[10]].t2d[2],  tabl[2 * j][index[8]].t2d[2], 
                                      tabl[2 * j][index[6]].t2d[2],  tabl[2 * j][index[4]].t2d[2], 
                                      tabl[2 * j][index[2]].t2d[2],  tabl[2 * j][index[0]].t2d[2]);
    lo8->t2d[3]   = _mm512_set_epi64(tabl[2 * j][index[14]].t2d[3], tabl[2 * j][index[12]].t2d[3], 
                                     tabl[2 * j][index[10]].t2d[3],  tabl[2 * j][index[8]].t2d[3], 
                                      tabl[2 * j][index[6]].t2d[3],  tabl[2 * j][index[4]].t2d[3], 
                                      tabl[2 * j][index[2]].t2d[3],  tabl[2 * j][index[0]].t2d[3]);
    lo8->t2d[4]   = _mm512_set_epi64(tabl[2 * j][index[14]].t2d[4], tabl[2 * j][index[12]].t2d[4], 
                                     tabl[2 * j][index[10]].t2d[4],  tabl[2 * j][index[8]].t2d[4], 
                                      tabl[2 * j][index[6]].t2d[4],  tabl[2 * j][index[4]].t2d[4], 
                                      tabl[2 * j][index[2]].t2d[4],  tabl[2 * j][index[0]].t2d[4]);

    hi8->ysubx[0] = _mm512_set_epi64(tabl[2 * j + 1][index[15]].ysubx[0], tabl[2 * j + 1][index[13]].ysubx[0], 
                                     tabl[2 * j + 1][index[11]].ysubx[0],  tabl[2 * j + 1][index[9]].ysubx[0], 
                                      tabl[2 * j + 1][index[7]].ysubx[0],  tabl[2 * j + 1][index[5]].ysubx[0], 
                                      tabl[2 * j + 1][index[3]].ysubx[0],  tabl[2 * j + 1][index[1]].ysubx[0]);
    hi8->ysubx[1] = _mm512_set_epi64(tabl[2 * j + 1][index[15]].ysubx[1], tabl[2 * j + 1][index[13]].ysubx[1], 
                                     tabl[2 * j + 1][index[11]].ysubx[1],  tabl[2 * j + 1][index[9]].ysubx[1], 
                                      tabl[2 * j + 1][index[7]].ysubx[1],  tabl[2 * j + 1][index[5]].ysubx[1], 
                                      tabl[2 * j + 1][index[3]].ysubx[1],  tabl[2 * j + 1][index[1]].ysubx[1]);
    hi8->ysubx[2] = _mm512_set_epi64(tabl[2 * j + 1][index[15]].ysubx[2], tabl[2 * j + 1][index[13]].ysubx[2], 
                                     tabl[2 * j + 1][index[11]].ysubx[2],  tabl[2 * j + 1][index[9]].ysubx[2], 
                                      tabl[2 * j + 1][index[7]].ysubx[2],  tabl[2 * j + 1][index[5]].ysubx[2], 
                                      tabl[2 * j + 1][index[3]].ysubx[2],  tabl[2 * j + 1][index[1]].ysubx[2]);
    hi8->ysubx[3] = _mm512_set_epi64(tabl[2 * j + 1][index[15]].ysubx[3], tabl[2 * j + 1][index[13]].ysubx[3], 
                                     tabl[2 * j + 1][index[11]].ysubx[3],  tabl[2 * j + 1][index[9]].ysubx[3], 
                                      tabl[2 * j + 1][index[7]].ysubx[3],  tabl[2 * j + 1][index[5]].ysubx[3], 
                                      tabl[2 * j + 1][index[3]].ysubx[3],  tabl[2 * j + 1][index[1]].ysubx[3]);
    hi8->ysubx[4] = _mm512_set_epi64(tabl[2 * j + 1][index[15]].ysubx[4], tabl[2 * j + 1][index[13]].ysubx[4], 
                                     tabl[2 * j + 1][index[11]].ysubx[4],  tabl[2 * j + 1][index[9]].ysubx[4], 
                                      tabl[2 * j + 1][index[7]].ysubx[4],  tabl[2 * j + 1][index[5]].ysubx[4], 
                                      tabl[2 * j + 1][index[3]].ysubx[4],  tabl[2 * j + 1][index[1]].ysubx[4]);

    hi8->yaddx[0] = _mm512_set_epi64(tabl[2 * j + 1][index[15]].yaddx[0], tabl[2 * j + 1][index[13]].yaddx[0], 
                                     tabl[2 * j + 1][index[11]].yaddx[0],  tabl[2 * j + 1][index[9]].yaddx[0], 
                                      tabl[2 * j + 1][index[7]].yaddx[0],  tabl[2 * j + 1][index[5]].yaddx[0], 
                                      tabl[2 * j + 1][index[3]].yaddx[0],  tabl[2 * j + 1][index[1]].yaddx[0]);
    hi8->yaddx[1] = _mm512_set_epi64(tabl[2 * j + 1][index[15]].yaddx[1], tabl[2 * j + 1][index[13]].yaddx[1], 
                                     tabl[2 * j + 1][index[11]].yaddx[1],  tabl[2 * j + 1][index[9]].yaddx[1], 
                                      tabl[2 * j + 1][index[7]].yaddx[1],  tabl[2 * j + 1][index[5]].yaddx[1], 
                                      tabl[2 * j + 1][index[3]].yaddx[1],  tabl[2 * j + 1][index[1]].yaddx[1]);
    hi8->yaddx[2] = _mm512_set_epi64(tabl[2 * j + 1][index[15]].yaddx[2], tabl[2 * j + 1][index[13]].yaddx[2], 
                                     tabl[2 * j + 1][index[11]].yaddx[2],  tabl[2 * j + 1][index[9]].yaddx[2], 
                                      tabl[2 * j + 1][index[7]].yaddx[2],  tabl[2 * j + 1][index[5]].yaddx[2], 
                                      tabl[2 * j + 1][index[3]].yaddx[2],  tabl[2 * j + 1][index[1]].yaddx[2]);
    hi8->yaddx[3] = _mm512_set_epi64(tabl[2 * j + 1][index[15]].yaddx[3], tabl[2 * j + 1][index[13]].yaddx[3], 
                                     tabl[2 * j + 1][index[11]].yaddx[3],  tabl[2 * j + 1][index[9]].yaddx[3], 
                                      tabl[2 * j + 1][index[7]].yaddx[3],  tabl[2 * j + 1][index[5]].yaddx[3], 
                                      tabl[2 * j + 1][index[3]].yaddx[3],  tabl[2 * j + 1][index[1]].yaddx[3]);
    hi8->yaddx[4] = _mm512_set_epi64(tabl[2 * j + 1][index[15]].yaddx[4], tabl[2 * j + 1][index[13]].yaddx[4], 
                                     tabl[2 * j + 1][index[11]].yaddx[4],  tabl[2 * j + 1][index[9]].yaddx[4], 
                                      tabl[2 * j + 1][index[7]].yaddx[4],  tabl[2 * j + 1][index[5]].yaddx[4], 
                                      tabl[2 * j + 1][index[3]].yaddx[4],  tabl[2 * j + 1][index[1]].yaddx[4]);

    hi8->t2d[0]   = _mm512_set_epi64(tabl[2 * j + 1][index[15]].t2d[0], tabl[2 * j + 1][index[13]].t2d[0], 
                                     tabl[2 * j + 1][index[11]].t2d[0],  tabl[2 * j + 1][index[9]].t2d[0], 
                                      tabl[2 * j + 1][index[7]].t2d[0],  tabl[2 * j + 1][index[5]].t2d[0], 
                                      tabl[2 * j + 1][index[3]].t2d[0],  tabl[2 * j + 1][index[1]].t2d[0]);
    hi8->t2d[1]   = _mm512_set_epi64(tabl[2 * j + 1][index[15]].t2d[1], tabl[2 * j + 1][index[13]].t2d[1], 
                                     tabl[2 * j + 1][index[11]].t2d[1],  tabl[2 * j + 1][index[9]].t2d[1], 
                                      tabl[2 * j + 1][index[7]].t2d[1],  tabl[2 * j + 1][index[5]].t2d[1], 
                                      tabl[2 * j + 1][index[3]].t2d[1],  tabl[2 * j + 1][index[1]].t2d[1]);
    hi8->t2d[2]   = _mm512_set_epi64(tabl[2 * j + 1][index[15]].t2d[2], tabl[2 * j + 1][index[13]].t2d[2], 
                                     tabl[2 * j + 1][index[11]].t2d[2],  tabl[2 * j + 1][index[9]].t2d[2], 
                                      tabl[2 * j + 1][index[7]].t2d[2],  tabl[2 * j + 1][index[5]].t2d[2], 
                                      tabl[2 * j + 1][index[3]].t2d[2],  tabl[2 * j + 1][index[1]].t2d[2]);
    hi8->t2d[3]   = _mm512_set_epi64(tabl[2 * j + 1][index[15]].t2d[3], tabl[2 * j + 1][index[13]].t2d[3], 
                                     tabl[2 * j + 1][index[11]].t2d[3],  tabl[2 * j + 1][index[9]].t2d[3], 
                                      tabl[2 * j + 1][index[7]].t2d[3],  tabl[2 * j + 1][index[5]].t2d[3], 
                                      tabl[2 * j + 1][index[3]].t2d[3],  tabl[2 * j + 1][index[1]].t2d[3]);
    hi8->t2d[4]   = _mm512_set_epi64(tabl[2 * j + 1][index[15]].t2d[4], tabl[2 * j + 1][index[13]].t2d[4], 
                                     tabl[2 * j + 1][index[11]].t2d[4],  tabl[2 * j + 1][index[9]].t2d[4], 
                                      tabl[2 * j + 1][index[7]].t2d[4],  tabl[2 * j + 1][index[5]].t2d[4], 
                                      tabl[2 * j + 1][index[3]].t2d[4],  tabl[2 * j + 1][index[1]].t2d[4]);

    // for(uint32_t i = 0; i < 5; ++i)
    // {
    //     lo8->ysubx[i] = _mm512_set_epi64(tabl[2 * j][index[14]].ysubx[i], tabl[2 * j][index[12]].ysubx[i], 
    //                                     tabl[2 * j][index[10]].ysubx[i], tabl[2 * j][index[8]].ysubx[i], 
    //                                     tabl[2 * j][index[6]].ysubx[i], tabl[2 * j][index[4]].ysubx[i], 
    //                                     tabl[2 * j][index[2]].ysubx[i], tabl[2 * j][index[0]].ysubx[i]);
    //     lo8->yaddx[i] = _mm512_set_epi64(tabl[2 * j][index[14]].yaddx[i], tabl[2 * j][index[12]].yaddx[i], 
    //                                     tabl[2 * j][index[10]].yaddx[i], tabl[2 * j][index[8]].yaddx[i], 
    //                                     tabl[2 * j][index[6]].yaddx[i], tabl[2 * j][index[4]].yaddx[i], 
    //                                     tabl[2 * j][index[2]].yaddx[i], tabl[2 * j][index[0]].yaddx[i]);
    //     lo8->t2d[i]   = _mm512_set_epi64(tabl[2 * j][index[14]].t2d[i], tabl[2 * j][index[12]].t2d[i], 
    //                                     tabl[2 * j][index[10]].t2d[i], tabl[2 * j][index[8]].t2d[i], 
    //                                     tabl[2 * j][index[6]].t2d[i], tabl[2 * j][index[4]].t2d[i], 
    //                                     tabl[2 * j][index[2]].t2d[i], tabl[2 * j][index[0]].t2d[i]);

    //     hi8->ysubx[i] = _mm512_set_epi64(tabl[2 * j + 1][index[15]].ysubx[i], tabl[2 * j + 1][index[13]].ysubx[i], 
    //                                     tabl[2 * j + 1][index[11]].ysubx[i], tabl[2 * j + 1][index[9]].ysubx[i], 
    //                                     tabl[2 * j + 1][index[7]].ysubx[i], tabl[2 * j + 1][index[5]].ysubx[i], 
    //                                     tabl[2 * j + 1][index[3]].ysubx[i], tabl[2 * j + 1][index[1]].ysubx[i]);
    //     hi8->yaddx[i] = _mm512_set_epi64(tabl[2 * j + 1][index[15]].yaddx[i], tabl[2 * j + 1][index[13]].yaddx[i], 
    //                                     tabl[2 * j + 1][index[11]].yaddx[i], tabl[2 * j + 1][index[9]].yaddx[i], 
    //                                     tabl[2 * j + 1][index[7]].yaddx[i], tabl[2 * j + 1][index[5]].yaddx[i], 
    //                                     tabl[2 * j + 1][index[3]].yaddx[i], tabl[2 * j + 1][index[1]].yaddx[i]);
    //     hi8->t2d[i]   = _mm512_set_epi64(tabl[2 * j + 1][index[15]].t2d[i], tabl[2 * j + 1][index[13]].t2d[i], 
    //                                     tabl[2 * j + 1][index[11]].t2d[i], tabl[2 * j + 1][index[9]].t2d[i], 
    //                                     tabl[2 * j + 1][index[7]].t2d[i], tabl[2 * j + 1][index[5]].t2d[i], 
    //                                     tabl[2 * j + 1][index[3]].t2d[i], tabl[2 * j + 1][index[1]].t2d[i]);
    // }

    // for(uint32_t i = 0; i < 5; ++i)
    // {
    //     lo8->ysubx[i] = _mm512_set_epi64(tabl[2 * j][i].ysubx[index[14]], tabl[2 * j][i].ysubx[index[12]], 
    //                                     tabl[2 * j][i].ysubx[index[10]], tabl[2 * j][i].ysubx[index[8]], 
    //                                     tabl[2 * j][i].ysubx[index[6]], tabl[2 * j][i].ysubx[index[4]], 
    //                                     tabl[2 * j][i].ysubx[index[2]], tabl[2 * j][i].ysubx[index[0]]);
    //     lo8->yaddx[i] = _mm512_set_epi64(tabl[2 * j][i].yaddx[index[14]], tabl[2 * j][i].yaddx[index[12]], 
    //                                     tabl[2 * j][i].yaddx[index[10]], tabl[2 * j][i].yaddx[index[8]], 
    //                                     tabl[2 * j][i].yaddx[index[6]], tabl[2 * j][i].yaddx[index[4]], 
    //                                     tabl[2 * j][i].yaddx[index[2]], tabl[2 * j][i].yaddx[index[0]]);
    //     lo8->t2d[i]   = _mm512_set_epi64(tabl[2 * j][i].t2d[index[14]], tabl[2 * j][i].t2d[index[12]], 
    //                                     tabl[2 * j][i].t2d[index[10]], tabl[2 * j][i].t2d[index[8]], 
    //                                     tabl[2 * j][i].t2d[index[6]], tabl[2 * j][i].t2d[index[4]], 
    //                                     tabl[2 * j][i].t2d[index[2]], tabl[2 * j][i].t2d[index[0]]);

    //     hi8->ysubx[i] = _mm512_set_epi64(tabl[2 * j + 1][i].ysubx[index[15]], tabl[2 * j + 1][i].ysubx[index[13]], 
    //                                     tabl[2 * j + 1][i].ysubx[index[11]], tabl[2 * j + 1][i].ysubx[index[9]], 
    //                                     tabl[2 * j + 1][i].ysubx[index[7]], tabl[2 * j + 1][i].ysubx[index[5]], 
    //                                     tabl[2 * j + 1][i].ysubx[index[3]], tabl[2 * j + 1][i].ysubx[index[1]]);
    //     hi8->yaddx[i] = _mm512_set_epi64(tabl[2 * j + 1][i].yaddx[index[15]], tabl[2 * j + 1][i].yaddx[index[13]], 
    //                                     tabl[2 * j + 1][i].yaddx[index[11]], tabl[2 * j + 1][i].yaddx[index[9]], 
    //                                     tabl[2 * j + 1][i].yaddx[index[7]], tabl[2 * j + 1][i].yaddx[index[5]], 
    //                                     tabl[2 * j + 1][i].yaddx[index[3]], tabl[2 * j + 1][i].yaddx[index[1]]);
    //     hi8->t2d[i]   = _mm512_set_epi64(tabl[2 * j + 1][i].t2d[index[15]], tabl[2 * j + 1][i].t2d[index[13]], 
    //                                     tabl[2 * j + 1][i].t2d[index[11]], tabl[2 * j + 1][i].t2d[index[9]], 
    //                                     tabl[2 * j + 1][i].t2d[index[7]], tabl[2 * j + 1][i].t2d[index[5]], 
    //                                     tabl[2 * j + 1][i].t2d[index[3]], tabl[2 * j + 1][i].t2d[index[1]]);
    // }

}



/*
 * Usage: scalar mult basepoint with comb algorithm
 * Author: Mengqing Yang
 * Time: 2023/11/10 17:16 +8
 * Mail: yangmengqing1\@jd.com
 */
void 
ge25519_scalarmult_base_comb_mb(ge52_ext_mb *h, const U64 scalar[4])
{
    U64             S;
    ge52_precomp_mb t1, t2;
    ge52_p1p1_mb p1p1;
    ge52_homo_mb r;
    int32_t         j;

    neutral_ge52_ext_mb(h);

        for (j = 0; j < 4; ++j) {
            S = _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 3), _mm512_set1_epi32(1));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 6), _mm512_set1_epi32(2)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 9), _mm512_set1_epi32(4)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 12), _mm512_set1_epi32(8)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 15), _mm512_set1_epi32(16)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 18), _mm512_set1_epi32(32)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 21), _mm512_set1_epi32(64)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 24), _mm512_set1_epi32(128)));

            ge25519_comb_base_mb(&t1, &t2, j, S);
            ge52_add_precomp(&p1p1, h, &t1);
            ge52_p1p1_to_ext_mb(h, &p1p1);
            ge52_add_precomp(&p1p1, h, &t2);
            ge52_p1p1_to_ext_mb(h, &p1p1);
        }

        ge52_ext_to_homo_mb(&r, h);
        ge52_double(h, h);

        for (j = 0; j < 4; ++j) {
            S = _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 2), _mm512_set1_epi32(1));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 5), _mm512_set1_epi32(2)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 8), _mm512_set1_epi32(4)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 11), _mm512_set1_epi32(8)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 14), _mm512_set1_epi32(16)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 17), _mm512_set1_epi32(32)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 20), _mm512_set1_epi32(64)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 23), _mm512_set1_epi32(128)));

            ge25519_comb_base_mb(&t1, &t2, j, S);
            ge52_add_precomp(&p1p1, h, &t1);
            ge52_p1p1_to_ext_mb(h, &p1p1);
            ge52_add_precomp(&p1p1, h, &t2);
            ge52_p1p1_to_ext_mb(h, &p1p1);
        }

        ge52_ext_to_homo_mb(&r, h);
        ge52_double(h, h);

        for (j = 0; j < 4; ++j) {
            S = _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 1), _mm512_set1_epi32(1));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 4), _mm512_set1_epi32(2)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 7), _mm512_set1_epi32(4)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 10), _mm512_set1_epi32(8)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 13), _mm512_set1_epi32(16)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 16), _mm512_set1_epi32(32)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 19), _mm512_set1_epi32(64)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 22), _mm512_set1_epi32(128)));

            ge25519_comb_base_mb(&t1, &t2, j, S);
            ge52_add_precomp(&p1p1, h, &t1);
            ge52_p1p1_to_ext_mb(h, &p1p1);
            ge52_add_precomp(&p1p1, h, &t2);
            ge52_p1p1_to_ext_mb(h, &p1p1);
        }

        ge52_ext_to_homo_mb(&r, h);
        ge52_double(h, h);

        for (j = 0; j < 4; ++j) {
            S = _mm512_and_epi64(scalar[j], _mm512_set1_epi32(1));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 3), _mm512_set1_epi32(2)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 6), _mm512_set1_epi32(4)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 9), _mm512_set1_epi32(8)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 12), _mm512_set1_epi32(16)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 15), _mm512_set1_epi32(32)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 18), _mm512_set1_epi32(64)));
            S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[j], 21), _mm512_set1_epi32(128)));

            ge25519_comb_base_mb(&t1, &t2, j, S);
            ge52_add_precomp(&p1p1, h, &t1);
            ge52_p1p1_to_ext_mb(h, &p1p1);
            ge52_add_precomp(&p1p1, h, &t2);
            ge52_p1p1_to_ext_mb(h, &p1p1);
        }

}

/*
 * Usage: scalar mult basepoint with comb algorithm
 * Author: Mengqing Yang
 * Time: 2024/02/05 10:46 +8
 * Mail: yangmengqing1\@jd.com
 */
void 
ge25519_scalarmult_base_comb_mb_expand(ge52_ext_mb *h, const U64 scalar[4])
{
    U64             S;
    ge52_precomp_mb t1, t2;
    ge52_p1p1_mb p1p1;
    ge52_homo_mb r;

    neutral_ge52_ext_mb(h);

        S = _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 3), _mm512_set1_epi32(1));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 6), _mm512_set1_epi32(2)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 9), _mm512_set1_epi32(4)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 12), _mm512_set1_epi32(8)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 15), _mm512_set1_epi32(16)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 18), _mm512_set1_epi32(32)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 21), _mm512_set1_epi32(64)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 24), _mm512_set1_epi32(128)));
        ge25519_comb_base_mb(&t1, &t2, 0, S);
        ge52_add_precomp(&p1p1, h, &t1);
        ge52_p1p1_to_ext_mb(h, &p1p1);
        ge52_add_precomp(&p1p1, h, &t2);
        ge52_p1p1_to_ext_mb(h, &p1p1);

        S = _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 3), _mm512_set1_epi32(1));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 6), _mm512_set1_epi32(2)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 9), _mm512_set1_epi32(4)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 12), _mm512_set1_epi32(8)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 15), _mm512_set1_epi32(16)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 18), _mm512_set1_epi32(32)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 21), _mm512_set1_epi32(64)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 24), _mm512_set1_epi32(128)));
        ge25519_comb_base_mb(&t1, &t2, 1, S);
        ge52_add_precomp(&p1p1, h, &t1);
        ge52_p1p1_to_ext_mb(h, &p1p1);
        ge52_add_precomp(&p1p1, h, &t2);
        ge52_p1p1_to_ext_mb(h, &p1p1);

        S = _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 3), _mm512_set1_epi32(1));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 6), _mm512_set1_epi32(2)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 9), _mm512_set1_epi32(4)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 12), _mm512_set1_epi32(8)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 15), _mm512_set1_epi32(16)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 18), _mm512_set1_epi32(32)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 21), _mm512_set1_epi32(64)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 24), _mm512_set1_epi32(128)));
        ge25519_comb_base_mb(&t1, &t2, 2, S);
        ge52_add_precomp(&p1p1, h, &t1);
        ge52_p1p1_to_ext_mb(h, &p1p1);
        ge52_add_precomp(&p1p1, h, &t2);
        ge52_p1p1_to_ext_mb(h, &p1p1);

        S = _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 3), _mm512_set1_epi32(1));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 6), _mm512_set1_epi32(2)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 9), _mm512_set1_epi32(4)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 12), _mm512_set1_epi32(8)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 15), _mm512_set1_epi32(16)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 18), _mm512_set1_epi32(32)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 21), _mm512_set1_epi32(64)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 24), _mm512_set1_epi32(128)));
        ge25519_comb_base_mb(&t1, &t2, 3, S);
        ge52_add_precomp(&p1p1, h, &t1);
        ge52_p1p1_to_ext_mb(h, &p1p1);
        ge52_add_precomp(&p1p1, h, &t2);
        ge52_p1p1_to_ext_mb(h, &p1p1);

        ge52_ext_to_homo_mb(&r, h);
        ge52_double(h, h);

        S = _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 2), _mm512_set1_epi32(1));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 5), _mm512_set1_epi32(2)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 8), _mm512_set1_epi32(4)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 11), _mm512_set1_epi32(8)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 14), _mm512_set1_epi32(16)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 17), _mm512_set1_epi32(32)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 20), _mm512_set1_epi32(64)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 23), _mm512_set1_epi32(128)));
        ge25519_comb_base_mb(&t1, &t2, 0, S);
        ge52_add_precomp(&p1p1, h, &t1);
        ge52_p1p1_to_ext_mb(h, &p1p1);
        ge52_add_precomp(&p1p1, h, &t2);
        ge52_p1p1_to_ext_mb(h, &p1p1);

        S = _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 2), _mm512_set1_epi32(1));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 5), _mm512_set1_epi32(2)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 8), _mm512_set1_epi32(4)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 11), _mm512_set1_epi32(8)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 14), _mm512_set1_epi32(16)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 17), _mm512_set1_epi32(32)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 20), _mm512_set1_epi32(64)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 23), _mm512_set1_epi32(128)));
        ge25519_comb_base_mb(&t1, &t2, 1, S);
        ge52_add_precomp(&p1p1, h, &t1);
        ge52_p1p1_to_ext_mb(h, &p1p1);
        ge52_add_precomp(&p1p1, h, &t2);
        ge52_p1p1_to_ext_mb(h, &p1p1);

        S = _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 2), _mm512_set1_epi32(1));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 5), _mm512_set1_epi32(2)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 8), _mm512_set1_epi32(4)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 11), _mm512_set1_epi32(8)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 14), _mm512_set1_epi32(16)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 17), _mm512_set1_epi32(32)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 20), _mm512_set1_epi32(64)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 23), _mm512_set1_epi32(128)));
        ge25519_comb_base_mb(&t1, &t2, 2, S);
        ge52_add_precomp(&p1p1, h, &t1);
        ge52_p1p1_to_ext_mb(h, &p1p1);
        ge52_add_precomp(&p1p1, h, &t2);
        ge52_p1p1_to_ext_mb(h, &p1p1);

        S = _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 2), _mm512_set1_epi32(1));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 5), _mm512_set1_epi32(2)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 8), _mm512_set1_epi32(4)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 11), _mm512_set1_epi32(8)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 14), _mm512_set1_epi32(16)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 17), _mm512_set1_epi32(32)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 20), _mm512_set1_epi32(64)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 23), _mm512_set1_epi32(128)));
        ge25519_comb_base_mb(&t1, &t2, 3, S);
        ge52_add_precomp(&p1p1, h, &t1);
        ge52_p1p1_to_ext_mb(h, &p1p1);
        ge52_add_precomp(&p1p1, h, &t2);
        ge52_p1p1_to_ext_mb(h, &p1p1);

        ge52_ext_to_homo_mb(&r, h);
        ge52_double(h, h);

        S = _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 1), _mm512_set1_epi32(1));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 4), _mm512_set1_epi32(2)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 7), _mm512_set1_epi32(4)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 10), _mm512_set1_epi32(8)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 13), _mm512_set1_epi32(16)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 16), _mm512_set1_epi32(32)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 19), _mm512_set1_epi32(64)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 22), _mm512_set1_epi32(128)));
        ge25519_comb_base_mb(&t1, &t2, 0, S);
        ge52_add_precomp(&p1p1, h, &t1);
        ge52_p1p1_to_ext_mb(h, &p1p1);
        ge52_add_precomp(&p1p1, h, &t2);
        ge52_p1p1_to_ext_mb(h, &p1p1);

        S = _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 1), _mm512_set1_epi32(1));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 4), _mm512_set1_epi32(2)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 7), _mm512_set1_epi32(4)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 10), _mm512_set1_epi32(8)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 13), _mm512_set1_epi32(16)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 16), _mm512_set1_epi32(32)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 19), _mm512_set1_epi32(64)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 22), _mm512_set1_epi32(128)));
        ge25519_comb_base_mb(&t1, &t2, 1, S);
        ge52_add_precomp(&p1p1, h, &t1);
        ge52_p1p1_to_ext_mb(h, &p1p1);
        ge52_add_precomp(&p1p1, h, &t2);
        ge52_p1p1_to_ext_mb(h, &p1p1);

        S = _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 1), _mm512_set1_epi32(1));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 4), _mm512_set1_epi32(2)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 7), _mm512_set1_epi32(4)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 10), _mm512_set1_epi32(8)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 13), _mm512_set1_epi32(16)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 16), _mm512_set1_epi32(32)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 19), _mm512_set1_epi32(64)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 22), _mm512_set1_epi32(128)));
        ge25519_comb_base_mb(&t1, &t2, 2, S);
        ge52_add_precomp(&p1p1, h, &t1);
        ge52_p1p1_to_ext_mb(h, &p1p1);
        ge52_add_precomp(&p1p1, h, &t2);
        ge52_p1p1_to_ext_mb(h, &p1p1);

        S = _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 1), _mm512_set1_epi32(1));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 4), _mm512_set1_epi32(2)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 7), _mm512_set1_epi32(4)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 10), _mm512_set1_epi32(8)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 13), _mm512_set1_epi32(16)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 16), _mm512_set1_epi32(32)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 19), _mm512_set1_epi32(64)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 22), _mm512_set1_epi32(128)));
        ge25519_comb_base_mb(&t1, &t2, 3, S);
        ge52_add_precomp(&p1p1, h, &t1);
        ge52_p1p1_to_ext_mb(h, &p1p1);
        ge52_add_precomp(&p1p1, h, &t2);
        ge52_p1p1_to_ext_mb(h, &p1p1);

        ge52_ext_to_homo_mb(&r, h);
        ge52_double(h, h);

        S = _mm512_and_epi64(scalar[0], _mm512_set1_epi32(1));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 3), _mm512_set1_epi32(2)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 6), _mm512_set1_epi32(4)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 9), _mm512_set1_epi32(8)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 12), _mm512_set1_epi32(16)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 15), _mm512_set1_epi32(32)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 18), _mm512_set1_epi32(64)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[0], 21), _mm512_set1_epi32(128)));
        ge25519_comb_base_mb(&t1, &t2, 0, S);
        ge52_add_precomp(&p1p1, h, &t1);
        ge52_p1p1_to_ext_mb(h, &p1p1);
        ge52_add_precomp(&p1p1, h, &t2);
        ge52_p1p1_to_ext_mb(h, &p1p1);

        S = _mm512_and_epi64(scalar[1], _mm512_set1_epi32(1));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 3), _mm512_set1_epi32(2)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 6), _mm512_set1_epi32(4)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 9), _mm512_set1_epi32(8)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 12), _mm512_set1_epi32(16)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 15), _mm512_set1_epi32(32)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 18), _mm512_set1_epi32(64)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[1], 21), _mm512_set1_epi32(128)));
        ge25519_comb_base_mb(&t1, &t2, 1, S);
        ge52_add_precomp(&p1p1, h, &t1);
        ge52_p1p1_to_ext_mb(h, &p1p1);
        ge52_add_precomp(&p1p1, h, &t2);
        ge52_p1p1_to_ext_mb(h, &p1p1);

        S = _mm512_and_epi64(scalar[2], _mm512_set1_epi32(1));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 3), _mm512_set1_epi32(2)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 6), _mm512_set1_epi32(4)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 9), _mm512_set1_epi32(8)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 12), _mm512_set1_epi32(16)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 15), _mm512_set1_epi32(32)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 18), _mm512_set1_epi32(64)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[2], 21), _mm512_set1_epi32(128)));
        ge25519_comb_base_mb(&t1, &t2, 2, S);
        ge52_add_precomp(&p1p1, h, &t1);
        ge52_p1p1_to_ext_mb(h, &p1p1);
        ge52_add_precomp(&p1p1, h, &t2);
        ge52_p1p1_to_ext_mb(h, &p1p1);

        S = _mm512_and_epi64(scalar[3], _mm512_set1_epi32(1));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 3), _mm512_set1_epi32(2)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 6), _mm512_set1_epi32(4)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 9), _mm512_set1_epi32(8)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 12), _mm512_set1_epi32(16)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 15), _mm512_set1_epi32(32)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 18), _mm512_set1_epi32(64)));
        S = _mm512_or_epi32(S, _mm512_and_epi64(_mm512_srli_epi32(scalar[3], 21), _mm512_set1_epi32(128)));
        ge25519_comb_base_mb(&t1, &t2, 3, S);
        ge52_add_precomp(&p1p1, h, &t1);
        ge52_p1p1_to_ext_mb(h, &p1p1);
        ge52_add_precomp(&p1p1, h, &t2);
        ge52_p1p1_to_ext_mb(h, &p1p1);

}


void 
curve25519_scalarmult_basepoint_comb_mb(uint8_t* pk[8], const uint8_t* const sk[8])
{
    ge52_ext_mb p;
    fe52_mb yplusz, zminusy;
    __ALIGN64 U64 scalar[4];
    ifma_BNU_transpose_copy((uint64_t (*)[8])scalar, (const uint8_t (*)[32])sk);

    ge25519_scalarmult_base_comb_mb(&p, scalar);

    fe52_mb_add_mod25519(yplusz, p.Y, p.Z);
    fe52_mb_sub_mod25519(zminusy, p.Z, p.Y);
    fe52_mb_inv_mod25519(zminusy, zminusy);
    fe52_mb_mul_mod25519(yplusz, yplusz, zminusy);
    ifma_mb8_to_BNU((uint64_t (*)[4])pk, (const uint64_t (*)[8])yplusz);
}

// #include "curve25519-yx-AVX512-modL.h"
static void 
ed25519_basemult_comb(uint8_t out[8][32], const uint8_t key[32], uint8_t data[8][32]) {
    /* Convert key of array mod to AVX512 mod first */
    uint64_t key64[4];
    uint64_t key52[5];
    uint32_t i, j;
    for(i = 0; i < 4; ++i)
    {
        key64[i] = 0;
        for(j = 0; j < 8; ++j)
        {
            key64[i] |= ((uint64_t)(key[i * 8 + j])) << (j * 8);
        }
    }
    uint64_t mask = (1ULL << 52) - 1;
    uint64_t top_mask = (1ULL << 48) - 1;
    key52[0] =   key64[0]                            & mask;
    key52[1] = ((key64[0] >> 52) | (key64[1] << 12)) & mask;
    key52[2] = ((key64[1] >> 40) | (key64[2] << 24)) & mask;
    key52[3] = ((key64[2] >> 28) | (key64[3] << 36)) & mask;
    key52[4] =  (key64[3] >> 16)                     & top_mask;
    __ALIGN64 U64 key52_mb8[5] = {_mm512_set1_epi64(key52[0]), 
                                  _mm512_set1_epi64(key52[1]), 
                                  _mm512_set1_epi64(key52[2]), 
                                  _mm512_set1_epi64(key52[3]), 
                                  _mm512_set1_epi64(key52[4])};

    __ALIGN64 U64 data52_mb8[5];
    __ALIGN64 U64 scalar52_mb8[5];
    __ALIGN64 U64 scalar64_mb8[4];
    ifma_BNU_to_mb8((uint64_t (*)[8])data52_mb8, (const uint8_t (*)[32])data);
   
    sc25519_mul(scalar52_mb8, key52_mb8, data52_mb8);

    scalar64_mb8[0] = _mm512_or_epi64(                  scalar52_mb8[0],      _mm512_slli_epi64(scalar52_mb8[1], 52));
    scalar64_mb8[1] = _mm512_or_epi64(_mm512_srli_epi64(scalar52_mb8[1], 12), _mm512_slli_epi64(scalar52_mb8[2], 40));
    scalar64_mb8[2] = _mm512_or_epi64(_mm512_srli_epi64(scalar52_mb8[2], 24), _mm512_slli_epi64(scalar52_mb8[3], 28));
    scalar64_mb8[3] = _mm512_or_epi64(_mm512_srli_epi64(scalar52_mb8[3], 36), _mm512_slli_epi64(scalar52_mb8[4], 16));

    ge52_ext_mb p;
    fe52_mb yplusz, zminusy;

    ge25519_scalarmult_base_comb_mb(&p, scalar64_mb8);

    fe52_mb_add_mod25519(yplusz, p.Y, p.Z);
    fe52_mb_sub_mod25519(zminusy, p.Z, p.Y);
    fe52_mb_inv_mod25519(zminusy, zminusy);
    fe52_mb_mul_mod25519(yplusz, yplusz, zminusy);

    __m512i P[5] = {_mm512_set1_epi64(PRIME25519_LO), 
                      _mm512_set1_epi64(PRIME25519_MID), 
                      _mm512_set1_epi64(PRIME25519_MID), 
                      _mm512_set1_epi64(PRIME25519_MID), 
                      _mm512_set1_epi64(PRIME25519_HI)};

    U64 borrow = get_zero64();
    fe52_mb cr;
    for(i = 0; i < 5; ++i)
    {
        borrow = sub64(yplusz[i], add64(P[i], _mm512_srli_epi64(borrow, 63)));
        cr[i] = and64_const(borrow, (1ULL << 52) - 1);
    }

    U64 underflow_mask = sub64(_mm512_xor_epi64(_mm512_srli_epi64(borrow, 63), _mm512_set1_epi64(1)), _mm512_set1_epi64(1));

    borrow = get_zero64();
    for(i = 0; i < 5; ++i)
    {
        borrow = add64(add64(cr[i], _mm512_srli_epi64(borrow, 52)), and64(P[i], underflow_mask));
        cr[i] = and64_const(borrow, (1ULL << 52) - 1);
    }

    ifma_mb8_to_BNU((uint64_t (*)[4])out, (const uint64_t (*)[8])cr);

    /* clear computed intermediate keys */
    MB_FUNC_NAME(zero_)((uint64_t (*)[8])key52_mb8, sizeof(key52_mb8)/sizeof(U64));
    MB_FUNC_NAME(zero_)((uint64_t (*)[8])scalar52_mb8, sizeof(scalar52_mb8)/sizeof(U64));
    MB_FUNC_NAME(zero_)((uint64_t (*)[8])scalar64_mb8, sizeof(scalar64_mb8)/sizeof(U64));
    for(i = 0; i < 4; ++i)
    {
        key64[i] = 0;
    }
    for(i = 0; i < 5; ++i)
    {
        key52[i] = 0;
    }
}

// __INLINE void sc25519_sub(U64 out[], const U64 a[], const U64 b[]) {

//     __m512i p[5] = {_mm512_set1_epi64(PRIME25519_LO), 
//                       _mm512_set1_epi64(PRIME25519_MID), 
//                       _mm512_set1_epi64(PRIME25519_MID), 
//                       _mm512_set1_epi64(PRIME25519_MID), 
//                       _mm512_set1_epi64(PRIME25519_HI)};

//     U64 *ca = (U64*) a;
//     U64 *cb = (U64*) b;
//     U64 *cr = (U64*) out;

//     U64 borrow = get_zero64();
//     uint32_t i;
//     for(i = 0; i < 5; ++i)
//     {
//         borrow = sub64(ca[i], add64(cb[i], _mm512_srli_epi64(borrow, 63)));
//         cr[i] = and64_const(borrow, (1ULL << 52) - 1);
//     }

//     U64 underflow_mask = sub64(_mm512_xor_epi64(_mm512_srli_epi64(borrow, 63), _mm512_set1_epi64(1)), _mm512_set1_epi64(1));

//     borrow = get_zero64();
//     for(i = 0; i < 5; ++i)
//     {
//         borrow = add64(add64(cr[i], _mm512_srli_epi64(borrow, 52)), and64(L[i], underflow_mask));
//         cr[i] = and64_const(borrow, (1ULL << 52) - 1);
//     }
// }