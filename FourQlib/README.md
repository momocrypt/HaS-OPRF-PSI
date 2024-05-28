## HasOPRFPSI based on FourQlib

## Contents

* [`FourQ_HaS_OPRF_PSI`](HaS_OPRF_PSI/FourQlib/FourQ_64bit_and_portable/): Given APIs and detailed implimentations in FourQ_HaS_OPRF_PSI.h and FourQ_HaS_OPRF_PSI.c
* [`oprf_psi_tests`](HaS_OPRF_PSI/FourQlib/FourQ_64bit_and_portable/tests): Given comparisons of the classic 2HashOPRFPSI and HaSOPRFPSI in oprf_psi_tests.c, including the correctness and benchmarks tests.

**Quick Start** cd in the dir HaS_OPRF_PSI/FourQlib/FourQ_64bit_and_portable/tests, one can quickly test the classic 2HashOPRFPSI and HaSOPRFPSI on any Intel x86 skylake or icelake microarchitectures by executing:

```sh
$ gcc -O3 -o oprf_psi_tests ../AMD64/fp2_1271_AVX2.S ../eccp2.c ../crypto_util.c ../eccp2_core.c ../eccp2_no_endo.c ../hash_to_curve.c ../kex.c ../schnorrq.c ../../sha512/sha512.c ../../blake2b/blake2b-ref.c ../../random/random.c test_extras.c -D USE_ENDO -D _AVX2_ -mavx2 oprf_psi_tests.c
```

and then run 

```sh
$ ./oprf_psi_tests
```

to verify the correctness and benchmarks.

* To imply `_asm_` optimization, add `-D _ASM_`.
* `-D USE_ENDO` is used for the endomorphism optimization in FourQ curve.
* `-D _AVX2_ -mavx2` is used for the AVX2 optimization (while only in table lookup) in FourQ curve.