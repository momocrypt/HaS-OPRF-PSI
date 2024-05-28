/***********************************************************************************
* FourQlib: a high-performance crypto library based on the elliptic curve FourQ
*
*    Copyright (c) Microsoft Corporation. All rights reserved.
*
* Abstract: testing code for cryptographic functions based on FourQ 
************************************************************************************/   

#include "../FourQ_api.h"
#include "../FourQ_params.h"
#include "../../random/random.h"
#include "../../sha512/sha512.h"
#include "../../blake2b/blake2.h"
#include "test_extras.h"
#include <stdio.h>


// Benchmark and test parameters  
#if defined(GENERIC_IMPLEMENTATION)
    #define BENCH_LOOPS       100       // Number of iterations per bench
    #define TEST_LOOPS        100       // Number of iterations per test
#else 
    #define BENCH_LOOPS       10000
    #define TEST_LOOPS        1000
#endif

void printPoint(point_t BB) 
{
    printf("Px = %0lx, %0lx, %0lx, %0lx\n", BB->x[0][0], BB->x[0][1], BB->x[1][0], BB->x[1][1]);
    printf("Py = %0lx, %0lx, %0lx, %0lx\n", BB->y[0][0], BB->y[0][1], BB->y[1][0], BB->y[1][1]);
}

ECCRYPTO_STATUS SchnorrQ_test()
{ // Test the SchnorrQ digital signature scheme
    int n, passed;       
    void *msg = NULL; 
    unsigned int len, valid = false;
    unsigned char SecretKey[32], PublicKey[32], Signature[64];
    ECCRYPTO_STATUS Status = ECCRYPTO_SUCCESS;

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Testing the SchnorrQ signature scheme: \n\n"); 

    passed = 1;
    for (n = 0; n < TEST_LOOPS; n++)
    {    
        // Signature key generation
        Status = SchnorrQ_FullKeyGeneration(SecretKey, PublicKey);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }  

        // Signature computation
        msg = "a";  
        len = 1;
        Status = SchnorrQ_Sign(SecretKey, PublicKey, msg, len, Signature);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }    

        // Valid signature test
        Status = SchnorrQ_Verify(PublicKey, msg, len, Signature, &valid);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }    
        if (valid == false) {
            passed = 0;
            break;
        }

        // Invalid signature test (flipping one bit of the message)
        msg = "b";  
        Status = SchnorrQ_Verify(PublicKey, msg, len, Signature, &valid);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }    
        if (valid == true) {
            passed = 0;
            break;
        }
    } 
    if (passed==1) printf("  Signature tests.................................................................. PASSED");
    else { printf("  Signature tests... FAILED"); printf("\n"); Status = ECCRYPTO_ERROR_SIGNATURE_VERIFICATION; }
    printf("\n");
    
    return Status;
}


ECCRYPTO_STATUS SchnorrQ_run()
{ // Benchmark the SchnorrQ digital signature scheme 
    int n;
    unsigned long long cycles, cycles1, cycles2;   
    void *msg = NULL;
    unsigned int len = 0, valid = false;
    unsigned char SecretKey[32], PublicKey[32], Signature[64];
    ECCRYPTO_STATUS Status = ECCRYPTO_SUCCESS;

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Benchmarking the SchnorrQ signature scheme: \n\n"); 
    
    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles(); 
        Status = SchnorrQ_FullKeyGeneration(SecretKey, PublicKey);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }    
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  SchnorrQ's key generation runs in ............................................... %8lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");
    
    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles(); 
        Status = SchnorrQ_Sign(SecretKey, PublicKey, msg, len, Signature);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }    
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  SchnorrQ's signing runs in ...................................................... %8lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");
    
    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles(); 
        Status = SchnorrQ_Verify(PublicKey, msg, len, Signature, &valid);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }    
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  SchnorrQ's verification runs in ................................................. %8lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");
    
    return Status;
}


ECCRYPTO_STATUS compressedkex_test()
{ // Test ECDH key exchange based on FourQ
    int n, passed;
    unsigned int i;
    unsigned char SecretKeyA[32], PublicKeyA[32], SecretAgreementA[32];
    unsigned char SecretKeyB[32], PublicKeyB[32], SecretAgreementB[32];
    ECCRYPTO_STATUS Status = ECCRYPTO_SUCCESS;

    printf("\n--------------------------------------------------------------------------------------------------------\n\n");
    printf("Testing DH key exchange using compressed, 32-byte public keys: \n\n");

    passed = 1;
    for (n = 0; n < TEST_LOOPS; n++)
    {
        // Alice's keypair generation
        Status = CompressedKeyGeneration(SecretKeyA, PublicKeyA);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }
        // Bob's keypair generation
        Status = CompressedKeyGeneration(SecretKeyB, PublicKeyB);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }

        // Alice's shared secret computation
        Status = CompressedSecretAgreement(SecretKeyA, PublicKeyB, SecretAgreementA);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }
        // Bob's shared secret computation
        Status = CompressedSecretAgreement(SecretKeyB, PublicKeyA, SecretAgreementB);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }

        for (i = 0; i < 32; i++) {
            if (SecretAgreementA[i] != SecretAgreementB[i]) {
                passed = 0;
                break;
            }
        }
    }
    if (passed==1) printf("  DH key exchange tests............................................................ PASSED");
    else { printf("  DH key exchange tests... FAILED"); printf("\n"); Status = ECCRYPTO_ERROR_SHARED_KEY; }
    printf("\n");

    return Status;
}


ECCRYPTO_STATUS compressedkex_run()
{ // Benchmark ECDH key exchange based on FourQ 
    int n;
    unsigned long long cycles, cycles1, cycles2;
    unsigned char SecretKeyA[32], PublicKeyA[32], SecretAgreementA[32];
    unsigned char SecretKeyB[32], PublicKeyB[32];
    ECCRYPTO_STATUS Status = ECCRYPTO_SUCCESS;

    printf("\n--------------------------------------------------------------------------------------------------------\n\n");
    printf("Benchmarking DH key exchange using compressed, 32-byte public keys: \n\n");

    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        Status = CompressedKeyGeneration(SecretKeyA, PublicKeyA);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }
        cycles2 = cpucycles();
        cycles = cycles + (cycles2 - cycles1);
    }
    printf("  Keypair generation runs in ...................................................... %8lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");
    
    Status = CompressedKeyGeneration(SecretKeyB, PublicKeyB);
    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        Status = CompressedSecretAgreement(SecretKeyA, PublicKeyB, SecretAgreementA);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }
        cycles2 = cpucycles();
        cycles = cycles + (cycles2 - cycles1);
    }
    printf("  Secret agreement runs in ........................................................ %8lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    return Status;
}


ECCRYPTO_STATUS kex_test()
{ // Test ECDH key exchange based on FourQ
    int n, passed;
    unsigned int i;
    unsigned char SecretKeyA[32], PublicKeyA[64], SecretAgreementA[32];
    unsigned char SecretKeyB[32], PublicKeyB[64], SecretAgreementB[32];
    ECCRYPTO_STATUS Status = ECCRYPTO_SUCCESS;

    printf("\n--------------------------------------------------------------------------------------------------------\n\n");
    printf("Testing DH key exchange using uncompressed, 64-byte public keys: \n\n");

    passed = 1;
    for (n = 0; n < TEST_LOOPS; n++)
    {
        // Alice's keypair generation
        Status = KeyGeneration(SecretKeyA, PublicKeyA);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }
        // Bob's keypair generation
        Status = KeyGeneration(SecretKeyB, PublicKeyB);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }

        // Alice's shared secret computation
        Status = SecretAgreement(SecretKeyA, PublicKeyB, SecretAgreementA);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }
        // Bob's shared secret computation
        Status = SecretAgreement(SecretKeyB, PublicKeyA, SecretAgreementB);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }

        for (i = 0; i < 32; i++) {
            if (SecretAgreementA[i] != SecretAgreementB[i]) {
                passed = 0;
                break;
            }
        }
    }
    if (passed==1) printf("  DH key exchange tests............................................................ PASSED");
    else { printf("  DH key exchange tests... FAILED"); printf("\n"); Status = ECCRYPTO_ERROR_SHARED_KEY; }
    printf("\n");

    return Status;
}


ECCRYPTO_STATUS kex_run()
{ // Benchmark ECDH key exchange based on FourQ 
    int n;
    unsigned long long cycles, cycles1, cycles2;
    unsigned char SecretKeyA[32], PublicKeyA[64], SecretAgreementA[32];
    unsigned char SecretKeyB[32], PublicKeyB[64];
    ECCRYPTO_STATUS Status = ECCRYPTO_SUCCESS;

    printf("\n--------------------------------------------------------------------------------------------------------\n\n");
    printf("Benchmarking DH key exchange using uncompressed, 64-byte public keys: \n\n");

    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        Status = KeyGeneration(SecretKeyA, PublicKeyA);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }
        cycles2 = cpucycles();
        cycles = cycles + (cycles2 - cycles1);
    }
    printf("  Keypair generation runs in ...................................................... %8lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    Status = KeyGeneration(SecretKeyB, PublicKeyB);
    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        Status = SecretAgreement(SecretKeyA, PublicKeyB, SecretAgreementA);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }
        cycles2 = cpucycles();
        cycles = cycles + (cycles2 - cycles1);
    }
    printf("  Secret agreement runs in ........................................................ %8lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    return Status;
}


ECCRYPTO_STATUS hash2curve_test()
{ // Test hashing to FourQ
    int n, passed;
    point_t P, Q;
    point_extproj_t R;
    unsigned char Value[32], HashedValue[64];
    f2elm_t* f2elmt = (f2elm_t*)&HashedValue[0];
    ECCRYPTO_STATUS Status = ECCRYPTO_SUCCESS;

    printf("\n--------------------------------------------------------------------------------------------------------\n\n");
    printf("Testing hashing to FourQ: \n\n");

    passed = 1;
    for (n = 0; n < TEST_LOOPS; n++)
    {
        RandomBytesFunction(Value, 32);
        CryptoHashFunction(Value, 32, HashedValue);
        mod1271(((felm_t*)f2elmt)[0]);
        mod1271(((felm_t*)f2elmt)[1]);

        // Hash GF(p^2) element to curve
        Status = HashToCurve((felm_t*)f2elmt, P);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }
        hash2curve_unsafe((felm_t*)f2elmt, Q);  // Non-constant-time version for testing        
        if (fp2compare64((uint64_t*)P->x,(uint64_t*)Q->x)!=0 || fp2compare64((uint64_t*)P->y,(uint64_t*)Q->y)!=0) { passed=0; break; }

        // Check if point is on the curve
        point_setup(P, R);
        if (!ecc_point_validate(R)) { passed=0; break; }
    }
    if (passed==1) printf("  Hash to FourQ tests.............................................................. PASSED");
    else { printf("  Hash to FourQ tests... FAILED"); printf("\n"); Status = ECCRYPTO_ERROR_HASH_TO_CURVE; }
    printf("\n");

    return Status;
}


ECCRYPTO_STATUS hash2curve_run()
{ // Benchmark hashing to FourQ 
    int n;
    unsigned long long cycles, cycles1, cycles2;
    point_t P;
    unsigned char Value[32], HashedValue[64];
    f2elm_t* f2elmt = (f2elm_t*)&HashedValue[0];
    ECCRYPTO_STATUS Status = ECCRYPTO_SUCCESS;

    printf("\n--------------------------------------------------------------------------------------------------------\n\n");
    printf("Benchmarking hashing to FourQ: \n\n");

    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        RandomBytesFunction(Value, 32);
        CryptoHashFunction(Value, 32, HashedValue);
        mod1271(((felm_t*)f2elmt)[0]);
        mod1271(((felm_t*)f2elmt)[1]);

        cycles1 = cpucycles();
        Status = HashToCurve((felm_t*)f2elmt, P);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }
        cycles2 = cpucycles();
        cycles = cycles + (cycles2 - cycles1);
    }
    printf("  Hashing to FourQ runs in ....................................................... %8lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    return Status;
}

void print256(unsigned char testout[32])
{
	for(int k = 0; k < 32; ++k)
	{
		printf("%02x ",testout[k]);
	}
	printf("\n");
}

ECCRYPTO_STATUS oprfpsi_test()
{ // Benchmark hashing to FourQ 
    int n;
    point_t P, Q, R;
    unsigned char encrypt1[32];
    unsigned char encrypt2[32];
    unsigned char Value[32], HashedValue[64];
    f2elm_t* f2elmt = (f2elm_t*)&HashedValue[0];
    ECCRYPTO_STATUS Status = ECCRYPTO_SUCCESS;
    
    int passed = 1;

    printf("\n--------------------------------------------------------------------------------------------------------\n\n");
    printf("Testing 2Hash-OPRF-PSI using FourQ: \n\n");

    for (n = 0; n < BENCH_LOOPS; n++)
    {
        //client: data generate
        RandomBytesFunction(Value, 32);
        CryptoHashFunction(Value, 32, HashedValue);
        mod1271(((felm_t*)f2elmt)[0]);
        mod1271(((felm_t*)f2elmt)[1]);

        Status = HashToCurve((felm_t*)f2elmt, P);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }
        //server: scalar generate
        uint64_t alpha[4];
        random_scalar_test(alpha);
        // uint64_t alpha[4] = {1, 0, 0, 0};
        //server: scalar multiplication
        ecc_mul(P, alpha, Q, true);
        encode(Q, encrypt1);

        //client: scalar generate
        uint64_t beta[4];
        uint64_t beta_inv[4];
        uint64_t beta_cache1[4];
        uint64_t beta_cache2[4];
        random_scalar_test(beta);
        //client: scalar inversion
        to_Montgomery(beta, beta_cache1);
        Montgomery_inversion_mod_order(beta_cache1, beta_cache2);
        from_Montgomery(beta_cache2, beta_inv);
        //client: scalar multiplication
        ecc_mul(P, beta, Q, true);
        encode(Q, encrypt2);

        //server: scalar multiplication
        decode(encrypt2, Q);
        ecc_mul(Q, alpha, R, false);
        encode(R, encrypt2);

        //client: scalar multiplication
        decode(encrypt2, R);
        ecc_mul(R, beta_inv, Q, false);
        encode(Q, encrypt2);

        if (fp2compare64((uint64_t*)encrypt1,(uint64_t*)encrypt2)!=0) { passed=0; break; }

    }
    if (passed==1) printf("  2Hash-OPRF-PSI tests.............................................................. PASSED");
    else { printf("  2Hash-OPRF-PSI tests... FAILED"); printf("\n"); Status = ECCRYPTO_ERROR_HASH_TO_CURVE; }
    printf("\n");
}

ECCRYPTO_STATUS oprfpsi_run()
{ // Benchmark hashing to FourQ 
    int n;
    unsigned long long cyclesA, cyclesB, cyclesC, cyclesD, cycles1, cycles2;
    point_t P, Q, R;
    unsigned char encrypt1[32];
    unsigned char encrypt2[32];
    unsigned char Value[32], HashedValue[64];
    f2elm_t* f2elmt = (f2elm_t*)&HashedValue[0];
    ECCRYPTO_STATUS Status = ECCRYPTO_SUCCESS;
    
    int passed = 1;

    printf("\n--------------------------------------------------------------------------------------------------------\n\n");
    printf("Benchmarking of 2Hash-OPRF-PSI using FourQ: \n\n");

    cyclesA = 0, cyclesB = 0, cyclesC = 0, cyclesD = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        //client: data generate
        RandomBytesFunction(Value, 32);

        CryptoHashFunction(Value, 32, HashedValue);

        mod1271(((felm_t*)f2elmt)[0]);
        mod1271(((felm_t*)f2elmt)[1]);
        Status = HashToCurve((felm_t*)f2elmt, P);
        if (Status != ECCRYPTO_SUCCESS) {
            return Status;
        }
        //server: scalar generate
        uint64_t alpha[4];
        random_scalar_test(alpha);
        // uint64_t alpha[4] = {1, 0, 0, 0};
        //server: scalar multiplication
        ecc_mul(P, alpha, Q, true);
        encode(Q, encrypt1);

        cycles2 = cpucycles();
        cyclesA = cyclesA + (cycles2 - cycles1);

        cycles1 = cpucycles();
        //client: scalar generate
        uint64_t beta[4];
        uint64_t beta_inv[4];
        uint64_t beta_cache1[4];
        uint64_t beta_cache2[4];
        random_scalar_test(beta);
        //client: scalar inversion
        to_Montgomery(beta, beta_cache1);
        Montgomery_inversion_mod_order(beta_cache1, beta_cache2);
        from_Montgomery(beta_cache2, beta_inv);
        //client: scalar multiplication
        ecc_mul(P, beta, Q, true);
        encode(Q, encrypt2);
        cycles2 = cpucycles();
        cyclesB = cyclesB + (cycles2 - cycles1);

        cycles1 = cpucycles();
        //server: scalar multiplication
        decode(encrypt2, Q);
        ecc_mul(Q, alpha, R, false);
        encode(R, encrypt2);
        cycles2 = cpucycles();
        cyclesC = cyclesC + (cycles2 - cycles1);

        cycles1 = cpucycles();
        //client: scalar multiplication
        decode(encrypt2, R);
        ecc_mul(R, beta_inv, Q, false);
        encode(Q, encrypt2);
        cycles2 = cpucycles();
        cyclesD = cyclesD + (cycles2 - cycles1);

    }
    printf("  2Hash-OPRF-PSI runs in .......Ev: %6lld,  Blind: %6lld,  BlindEv: %6lld,  Finalize: %6lld ", cyclesA/BENCH_LOOPS, cyclesB/BENCH_LOOPS, cyclesC/BENCH_LOOPS, cyclesD/BENCH_LOOPS); print_unit;
    printf("\n");

    return Status;
}

ECCRYPTO_STATUS haspsi_test()
{ // Benchmark hashing to FourQ 
    int n;
    point_t R, B1, B2, Y, N, P, Q;
    unsigned char encrypt1[32];
    unsigned char encrypt2[32];
    unsigned char Value[32], FirstHashValue[64], SecondHashValue[32];
    uint64_t data1[4], data2[4];
    ECCRYPTO_STATUS Status = ECCRYPTO_SUCCESS;
    
    int passed = 1;

    printf("\n--------------------------------------------------------------------------------------------------------\n\n");
    printf("Testing HaS-OPRF-PSI using FourQ: \n\n");

    //server: key generate
    uint64_t alpha[4], r[4];
    random_scalar_test(alpha);
    random_scalar_test(r);
    //server: key inversion
    uint64_t alpha_inv[4];
    uint64_t alpha_cache1[4];
    uint64_t alpha_cache2[4];
    to_Montgomery(alpha, alpha_cache1);
    Montgomery_inversion_mod_order(alpha_cache1, alpha_cache2);
    from_Montgomery(alpha_cache2, alpha_inv);
    //server: setup to compute R = g^{r·alpha} = g^{alpha_cache1}
    Scalar_mul(alpha, r, alpha_cache1);
    ecc_mul_comb(alpha_cache1, R);

    for (n = 0; n < BENCH_LOOPS; n++)
    // for (n = 0; n < 10; n++)
    {
        //data generate
        RandomBytesFunction(Value, 32);

        //compute 2 hashes
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
        // scalar = alpha_inv·H1(x) + r·H2(x)
        Scalar_mul(alpha_inv, data1, scalar1);
        Scalar_mul(r, data2, scalar2);
        add_mod_order(scalar1, scalar2, scalar);
        
        //server: base point scalar multiplication using mLSB-Comb
        // Y = g^{scalar}
        ecc_mul_comb(scalar, Y);
        encode(Y, encrypt1);

        //client: client key beta generate
        uint64_t beta[4];
        uint64_t beta_inv[4];
        uint64_t beta_cache1[4];
        uint64_t beta_cache2[4];
        random_scalar_test(beta);
        //client: scalar inversion
        to_Montgomery(beta, beta_cache1);
        Montgomery_inversion_mod_order(beta_cache1, beta_cache2);
        from_Montgomery(beta_cache2, beta_inv);
        //client: scalar multiplication
        Scalar_mul(beta, data1, scalar1);
        // ecc_mul_comb(scalar2, P);
        Scalar_mul(beta, data2, scalar2);
        // ecc_mul(R, scalar2, B1, false);
        // eccadd_affine(P, B1);
        ecc_mul_double(scalar1, R, scalar2, B1);
        encode(B1, encrypt2);

        //server: scalar multiplication
        decode(encrypt2, B1);
        ecc_mul(B1, alpha_inv, B2, false);
        encode(B2, encrypt2);

        //client: scalar multiplication
        decode(encrypt2, B2);
        ecc_mul(B2, beta_inv, N, false);
        encode(N, encrypt2);

        // printf("Num. %d: encrypt1: \n", n);
        // print256(encrypt1);
        // printf("Num. %d: encrypt2: \n", n);
        // print256(encrypt2);

        if (fp2compare64((uint64_t*)encrypt1,(uint64_t*)encrypt2)!=0) { passed=0; break; }

    }
    if (passed==1) printf("  HaS-OPRF-PSI tests.............................................................. PASSED");
    else { printf("  HaS-OPRF-PSI tests... FAILED"); printf("\n"); Status = ECCRYPTO_ERROR_HASH_TO_CURVE; }
    printf("\n");
}

ECCRYPTO_STATUS haspsi_run()
{ // Benchmark hashing to FourQ 
    int n;
    unsigned long long cyclesA, cyclesB, cyclesC, cyclesD, cycles1, cycles2;
    point_t R, B1, B2, Y, N, P, Q;
    unsigned char encrypt1[32];
    unsigned char encrypt2[32];
    unsigned char Value[32], FirstHashValue[64], SecondHashValue[32];
    uint64_t data1[4], data2[4];
    ECCRYPTO_STATUS Status = ECCRYPTO_SUCCESS;

    printf("\n--------------------------------------------------------------------------------------------------------\n\n");
    printf("Benchmarking of HaS-OPRF-PSI using FourQ: \n\n");

    cyclesA = 0, cyclesB = 0, cyclesC = 0, cyclesD = 0;    

    //server: sk alpha generate
    uint64_t alpha[4], r[4];
    random_scalar_test(alpha);
    random_scalar_test(r);
    //server: sk alpha inversion, alpha_inv = alpha^{-1}
    uint64_t alpha_inv[4];
    uint64_t alpha_cache1[4];
    uint64_t alpha_cache2[4];
    to_Montgomery(alpha, alpha_cache1);
    Montgomery_inversion_mod_order(alpha_cache1, alpha_cache2);
    from_Montgomery(alpha_cache2, alpha_inv);
    //server: setup to compute R = g^{r·alpha} = g^{alpha_cache1}
    Scalar_mul(alpha, r, alpha_cache1);
    ecc_mul_comb(alpha_cache1, R);

    //client: sk beta generate
    uint64_t beta[4];
    uint64_t beta_inv[4];
    uint64_t beta_cache1[4];
    uint64_t beta_cache2[4];
    random_scalar_test(beta);
    //client: sk beta inversion, beta_inv = beta^{-1}
    to_Montgomery(beta, beta_cache1);
    Montgomery_inversion_mod_order(beta_cache1, beta_cache2);
    from_Montgomery(beta_cache2, beta_inv);

    for (n = 0; n < BENCH_LOOPS; n++)
    {

        RandomBytesFunction(Value, 32);

        cycles1 = cpucycles();
        //compute 2 hashes
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
        cycles2 = cpucycles();
        cyclesA = cyclesA + (cycles2 - cycles1);

        cycles1 = cpucycles();
        //client: scalar multiplication
        //scalar1 = beta·H1(y)
        Scalar_mul(beta, data1, scalar1);
        //scalar2 = beta·H2(y)
        Scalar_mul(beta, data2, scalar2);
        // B1 = g^{beta·H1(y)}·R^{beta·H2(y)}
        ecc_mul_double(scalar1, R, scalar2, B1);
        encode(B1, encrypt2);
        cycles2 = cpucycles();
        cyclesB = cyclesB + (cycles2 - cycles1);

        cycles1 = cpucycles();
        //server: scalar multiplication
        // B2 = B1^{alpha_inv}
        decode(encrypt2, B1);
        ecc_mul(B1, alpha_inv, B2, false);
        encode(B2, encrypt2);
        cycles2 = cpucycles();
        cyclesC = cyclesC + (cycles2 - cycles1);

        cycles1 = cpucycles();
        //client: scalar multiplication
        // N = B2^{beta_inv}
        decode(encrypt2, B2);
        ecc_mul(B2, beta_inv, N, false);
        encode(N, encrypt2);
        cycles2 = cpucycles();
        cyclesD = cyclesD + (cycles2 - cycles1);
    }
    printf("  HaS-OPRF-PSI runs in .........Ev: %6lld,  Blind: %6lld,  BlindEv: %6lld,  Finalize: %6lld ", cyclesA/BENCH_LOOPS, cyclesB/BENCH_LOOPS, cyclesC/BENCH_LOOPS, cyclesD/BENCH_LOOPS); print_unit;
    printf("\n");

    return Status;
}

int main()
{
    ECCRYPTO_STATUS Status = ECCRYPTO_SUCCESS;
    
    Status = SchnorrQ_test();         // Test SchnorrQ signature scheme
    if (Status != ECCRYPTO_SUCCESS) {
        printf("\n\n   Error detected: %s \n\n", FourQ_get_error_message(Status));
        return false;
    }
    Status = SchnorrQ_run();          // Benchmark SchnorrQ signature scheme
    if (Status != ECCRYPTO_SUCCESS) {
        printf("\n\n   Error detected: %s \n\n", FourQ_get_error_message(Status));
        return false;
    }

    Status = compressedkex_test();    // Test Diffie-Hellman key exchange using compressed public keys
    if (Status != ECCRYPTO_SUCCESS) {
        printf("\n\n   Error detected: %s \n\n", FourQ_get_error_message(Status));
        return false;
    }
    Status = compressedkex_run();     // Benchmark Diffie-Hellman key exchange using compressed public keys
    if (Status != ECCRYPTO_SUCCESS) {
        printf("\n\n   Error detected: %s \n\n", FourQ_get_error_message(Status));
        return false;
    }

    Status = kex_test();              // Test Diffie-Hellman key exchange using uncompressed public keys
    if (Status != ECCRYPTO_SUCCESS) {
        printf("\n\n   Error detected: %s \n\n", FourQ_get_error_message(Status));
        return false;
    }
    Status = kex_run();               // Benchmark Diffie-Hellman key exchange using uncompressed public keys
    if (Status != ECCRYPTO_SUCCESS) {
        printf("\n\n   Error detected: %s \n\n", FourQ_get_error_message(Status));
        return false;
    }
    
    Status = hash2curve_test();       // Test hash to FourQ function
    if (Status != ECCRYPTO_SUCCESS) {
        printf("\n\n   Error detected: %s \n\n", FourQ_get_error_message(Status));
        return false;
    }
    Status = hash2curve_run();        // Benchmark hash to FourQ function
    if (Status != ECCRYPTO_SUCCESS) {
        printf("\n\n   Error detected: %s \n\n", FourQ_get_error_message(Status));
        return false;
    }
    Status = oprfpsi_test();        // Benchmark hash to FourQ function
    Status = oprfpsi_run();        // Benchmark hash to FourQ function
    if (Status != ECCRYPTO_SUCCESS) {
        printf("\n\n   Error detected: %s \n\n", FourQ_get_error_message(Status));
        return false;
    }

    Status = haspsi_test(); 
    Status = haspsi_run();        // Benchmark hash to FourQ function
    if (Status != ECCRYPTO_SUCCESS) {
        printf("\n\n   Error detected: %s \n\n", FourQ_get_error_message(Status));
        return false;
    }
    // testBaseMult();

    return true;
}