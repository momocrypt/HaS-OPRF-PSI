/***********************************************************************************
*               Hash as Scalar: A new OPRF-PSI Protocal

* Abstract: Details of each phase for HaS-OPRF-PSI protocol based on FourQ 
************************************************************************************/   

#include "FourQ_HaS_OPRF_PSI.h"
#include <stdio.h>

// For C++
#ifdef __cplusplus
extern "C" {
#endif

void random_scalar_generator(uint64_t* a)
{ // Generating a pseudo-random scalar value in [0, 2^256-1] 
  // NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.
    unsigned char* string = (unsigned char*)&a[0];
    unsigned int i;

    for (i = 0; i < (sizeof(uint64_t)*NWORDS64_ORDER); i++) {
        string[i] = (unsigned char)rand();             
    }
}

// Server setup: generate alpha, r; compute alpha_inv, R. Send R to Client.
void server_setup(digit_t *alpha, digit_t *alpha_inv, digit_t *r, point_affine *R)
{
    //server: sk alpha generate, alpha & r were previously declared in the prior program but cannot be assigned values.
    random_scalar_generator(alpha);
    random_scalar_generator(r);
    //server: sk alpha inversion, alpha_inv = alpha^{-1}
    uint64_t alpha_cache1[4];
    uint64_t alpha_cache2[4];
    to_Montgomery(alpha, alpha_cache1);
    Montgomery_inversion_mod_order(alpha_cache1, alpha_cache2);
    from_Montgomery(alpha_cache2, alpha_inv);
    //server: setup to compute R = g^{r·alpha} = g^{alpha_cache1}
    Scalar_mul(alpha, r, alpha_cache1);
    ecc_mul_comb(alpha_cache1, R);
}

// Client setup: generate beta.
void client_setup(digit_t *beta)
{
    //server: sk beta generate, beta was previously declared in the prior program but cannot be assigned values.
    random_scalar_generator(beta);
}

/*
 * Usage: Server Full Evaluation
 * Author: Mengqing Yang
 * Time: 2023/03/22 13:31 +8
 * Mail: yangmengqing\@pku.edu.cn
 * Input: server sk inversion $alpha_inv$, server random number $r$, original data $Value$.
 * Output: server Full Evaluation output $encrypt1$.
 * Y = g^{alpha_inv·H1(x) + r·H2(x)}
 */
void server_evaluation(unsigned char *Value, digit_t *alpha_inv, digit_t *r, unsigned char *encrypt1)
{
    point_t Y;
    unsigned char FirstHashValue[64], SecondHashValue[32];
    uint64_t data1[4], data2[4];

    //sha512:
    CryptoHashFunction(Value, 32, FirstHashValue);
    data1[0] = (uint64_t)FirstHashValue[0];
    data1[1] = (uint64_t)FirstHashValue[8];
    data1[2] = (uint64_t)FirstHashValue[16];
    data1[3] = (uint64_t)FirstHashValue[24];
    //blake2b:
    blake2b(SecondHashValue, 32, Value, 32, NULL, 0);
    data2[0] = (uint64_t)SecondHashValue[0];
    data2[1] = (uint64_t)SecondHashValue[8];
    data2[2] = (uint64_t)SecondHashValue[16];
    data2[3] = (uint64_t)SecondHashValue[24];
    //server: compute scalar
    uint64_t scalar1[4], scalar2[4], scalar[4];
    //server: scalar multiplications & addition
    //scalar = alpha_inv·H1(x) + r·H2(x)
    Scalar_mul(alpha_inv, data1, scalar1);
    Scalar_mul(r, data2, scalar2);
    add_mod_order(scalar1, scalar2, scalar);
    //server: base point scalar multiplication using mLSB-Comb
    // Y = g^{scalar}
    ecc_mul_comb(scalar, Y);
    encode(Y, encrypt1);
}

/*
 * Usage: Client Blinding
 * Author: Mengqing Yang
 * Time: 2023/03/22 14:12 +8
 * Mail: yangmengqing\@pku.edu.cn
 * Input: client sk $beta$, point $R$ received from server, original data $Value$.
 * Output: client Blinding output $encrypt2$.
 * 
 * The $B1 = g^{beta·H1(y)}·R^{beta·H2(y)}$ calculation temporarily uses FourQ native functions, 
 * because after workflow optimization, the Blind calculation time has no impact on the overall process. 
 */
void client_blinding(unsigned char *Value, digit_t *beta, point_affine *R, unsigned char *encrypt2)
{
    point_t B1;
    unsigned char FirstHashValue[64], SecondHashValue[32];
    uint64_t data1[4], data2[4];

    //sha512:
    CryptoHashFunction(Value, 32, FirstHashValue);
    data1[0] = (uint64_t)FirstHashValue[0];
    data1[1] = (uint64_t)FirstHashValue[8];
    data1[2] = (uint64_t)FirstHashValue[16];
    data1[3] = (uint64_t)FirstHashValue[24];
    //blake2b:
    blake2b(SecondHashValue, 32, Value, 32, NULL, 0);
    data2[0] = (uint64_t)SecondHashValue[0];
    data2[1] = (uint64_t)SecondHashValue[8];
    data2[2] = (uint64_t)SecondHashValue[16];
    data2[3] = (uint64_t)SecondHashValue[24];
    
    //client: scalar multiplication
    uint64_t scalar1[4], scalar2[4], scalar[4];
    //scalar1 = beta·H1(y)
    Scalar_mul(beta, data1, scalar1);
    //scalar2 = beta·H2(y)
    Scalar_mul(beta, data2, scalar2);
    // B1 = g^{beta·H1(y)}·R^{beta·H2(y)}
    ecc_mul_double(scalar1, R, scalar2, B1);
    encode(B1, encrypt2);
}

// Client sk inversion: compute beta_inv when Client finishing Blinding.
void client_inv(digit_t *beta, digit_t *beta_inv)
{
    uint64_t beta_cache1[4];
    uint64_t beta_cache2[4];
    //client: sk beta inversion, beta_inv = beta^{-1}
    to_Montgomery(beta, beta_cache1);
    Montgomery_inversion_mod_order(beta_cache1, beta_cache2);
    from_Montgomery(beta_cache2, beta_inv);
}

/*
 * Usage: Server BlindEvaluation
 * Author: Mengqing Yang
 * Time: 2023/03/22 15:22 +8
 * Mail: yangmengqing\@pku.edu.cn
 * Input: server sk inversion $alpha_inv$, string $encrypt2$ received from client, which is the y-coordinate encoded from point $B1$.
 * Output: server BlindEvaluation output $encrypt2$.
 * B2 = B1^{alpha_inv}
 */
void server_blindevaluation(digit_t *alpha_inv, unsigned char *encrypt2)
{
    //server: scalar multiplication
    point_t B1, B2;
    // B2 = B1^{alpha_inv}
    decode(encrypt2, B1);
    ecc_mul(B1, alpha_inv, B2, false);
    encode(B2, encrypt2);
}

/*
 * Usage: Client Finalizaiton
 * Author: Mengqing Yang
 * Time: 2023/03/22 16:39 +8
 * Mail: yangmengqing\@pku.edu.cn
 * Input: client sk inversion $beta_inv$, string $encrypt2$ received from server, which is the y-coordinate encoded from point $B2$.
 * Output: client Finalizaiton output $encrypt2$.
 * N = B2^{beta_inv}
 */
void client_finalization(digit_t *beta_inv, unsigned char *encrypt2)
{
    //client: scalar multiplication
    point_t B2, N;
    // N = B2^{beta_inv}
    decode(encrypt2, B2);
    ecc_mul(B2, beta_inv, N, false);
    encode(N, encrypt2);
}

#ifdef __cplusplus
}
#endif
