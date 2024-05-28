/***********************************************************************************
*               Hash as Scalar: A new OPRF-PSI Protocal

* Abstract: API of each phase for HaS-OPRF-PSI protocol based on FourQ 
************************************************************************************/   

#include "FourQ_api.h"
#include "FourQ_params.h"
#include "../random/random.h"
#include "../sha512/sha512.h"
#include "../blake2b/blake2.h"
#include <stdio.h>

// For C++
#ifdef __cplusplus
extern "C" {
#endif

void random_scalar_generator(uint64_t* a);

// Server setup: generate $alpha$, $r$; compute alpha_inv, $R$. Send $R$ to Client.
void server_setup(digit_t *alpha, digit_t *alpha_inv, digit_t *r, point_affine *R);

// Client setup: generate $beta$.
void client_setup(digit_t *beta);

// Server Full Evaluation: $Y = g^{alpha_inv·H1(x) + r·H2(x)}$. Send $encrypt1 = y(Y)$ to Server.
void server_evaluation(unsigned char *Value, digit_t *alpha_inv, digit_t *r, unsigned char *encrypt1);

// Client Blinding: B1 = g^{beta·H1(y)}·R^{beta·H2(y)}. Send $encrypt2 = y(B1)$ to Server.
void client_blinding(unsigned char *Value, digit_t *beta, point_affine *R, unsigned char *encrypt2);

// Client sk inversion: compute beta_inv when Client finishing Blinding.
void client_inv(digit_t *beta, digit_t *beta_inv);

// Server BlindEvaluation: B2 = B1^{alpha_inv}. Send $encrypt2 = y(B2)$ to Client.
void server_blindevaluation(digit_t *alpha_inv, unsigned char *encrypt2);

// Client Finalizaiton: N = B2^{beta_inv}. Hold $encrypt2 = y(N)$.
void client_finalization(digit_t *beta_inv, unsigned char *encrypt2);

#ifdef __cplusplus
}
#endif
