#include "../FourQ_api.h"
#include "../FourQ_params.h"
#include "../../random/random.h"
#include "../../sha512/sha512.h"
#include "test_extras.h"
#include <stdio.h>

// #include "../FourQ_internal.h"
// // #include "../FourQ_params.h"
// #include "../FourQ_tables.h"
// #if defined(GENERIC_IMPLEMENTATION)
//     #include "../generic/fp.h"
// #elif (TARGET == TARGET_AMD64)
//     #include "../AMD64/fp_x64.h"
// #elif (TARGET == TARGET_ARM64)
//     #include "../ARM64/fp_arm64.h"
// #endif

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

void testBaseMult()
{
    uint64_t k[4] = {1, 0, 0, 0};
    point_t Q, BB;
    unsigned int j, w = W_FIXEDBASE, v = V_FIXEDBASE, d = D_FIXEDBASE, e = E_FIXEDBASE;
    unsigned int digit = 0, digits[NBITS_ORDER_PLUS_ONE+(W_FIXEDBASE*V_FIXEDBASE)-1] = {0}; 
    digit_t temp[NWORDS_ORDER];
    point_extproj_t R;
    point_precomp_t S;
    int i, ii;

	modulo_order(k, temp);                                      // temp = k mod (order) 
	conversion_to_odd(temp, temp);                              // Converting scalar to odd using the prime subgroup order
	mLSB_set_recode((uint64_t*)temp, digits);                   // Scalar recoding

    // Extracting initial digit 
    digit = digits[w*d-1];
    for (i = (int)((w-1)*d-1); i >= (int)(2*d-1); i = i-d)           
    {
        digit = 2*digit + digits[i];
    }
    // Initialize R = (x+y,y-x,2dt) with a point from the table
	table_lookup_fixed_base(((point_precomp_t*)&FIXED_BASE_TABLE)+(v-1)*(1 << (w-1)), S, digit, digits[d-1]);
    R5_to_R1(S, R);                                             // Converting to representation (X:Y:1:Ta:Tb)

    eccnorm(R, BB);
    printPoint(BB);

    for (j = 0; j < (v-1); j++)
    {
        digit = digits[w*d-(j+1)*e-1];
        for (i = (int)((w-1)*d-(j+1)*e-1); i >= (int)(2*d-(j+1)*e-1); i = i-d)           
        {
            digit = 2*digit + digits[i];
        }
        // Extract point in (x+y,y-x,2dt) representation
        table_lookup_fixed_base(((point_precomp_t*)&FIXED_BASE_TABLE)+(v-j-2)*(1 << (w-1)), S, digit, digits[d-(j+1)*e-1]);
        eccmadd(S, R);                                          // R = R+S using representations (X,Y,Z,Ta,Tb) <- (X,Y,Z,Ta,Tb) + (x+y,y-x,2dt) 
        eccnorm(R, BB);
        printPoint(BB);
    }

    for (ii = (e-2); ii >= 0; ii--)
    {
        eccdouble(R);                                           // R = 2*R using representations (X,Y,Z,Ta,Tb) <- 2*(X,Y,Z)
        for (j = 0; j < v; j++)
        {
            digit = digits[w*d-j*e+ii-e];
            for (i = (int)((w-1)*d-j*e+ii-e); i >= (int)(2*d-j*e+ii-e); i = i-d)           
            {
                digit = 2*digit + digits[i];
            }
            // Extract point in (x+y,y-x,2dt) representation
            table_lookup_fixed_base(((point_precomp_t*)&FIXED_BASE_TABLE)+(v-j-1)*(1 << (w-1)), S, digit, digits[d-j*e+ii-e]);
            eccmadd(S, R);                                      // R = R+S using representations (X,Y,Z,Ta,Tb) <- (X,Y,Z,Ta,Tb) + (x+y,y-x,2dt)
        }        
    }     
    eccnorm(R, Q);

    printPoint(Q);
}

int main()
{
    testBaseMult();

    return true;
}