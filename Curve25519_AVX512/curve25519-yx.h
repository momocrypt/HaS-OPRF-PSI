#if !defined(ED25519_SUFFIX)
#define ED25519_SUFFIX 
#endif

#define ED25519_FN3(fn,suffix) fn##suffix
#define ED25519_FN2(fn,suffix) ED25519_FN3(fn,suffix)
#define ED25519_FN(fn)         ED25519_FN2(fn,ED25519_SUFFIX)

#include <assert.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/syscall.h>
#include <linux/perf_event.h>
#include <string.h>

static long long cpucycles(void)
{
  static int rdpmcworks = 1;
  long long result;

  while (rdpmcworks) {
    static int fdperf = -1;
    static struct perf_event_mmap_page *buf = 0;
#ifdef THREADING
    unsigned int seq;
    long long index;
    long long offset;
#endif

    if (fdperf == -1) {
      static struct perf_event_attr attr;
      memset(&attr,0,sizeof attr);
      attr.type = PERF_TYPE_HARDWARE;
      attr.config = PERF_COUNT_HW_CPU_CYCLES;
      attr.exclude_kernel = 1;
      fdperf = syscall(__NR_perf_event_open, &attr, 0, -1, -1, 0);
      if (fdperf == -1) {
        rdpmcworks = 0;
        break;
      }
      buf = mmap(NULL, sysconf(_SC_PAGESIZE), PROT_READ, MAP_SHARED, fdperf, 0);
      if (buf == 0) {
        rdpmcworks = 0;
        break;
      }
    }

#ifdef THREADING
    do {
      seq = buf->lock;
      asm volatile("" ::: "memory");
      index = buf->index;
      offset = buf->offset;
      asm volatile("rdpmc;shlq $32,%%rdx;orq %%rdx,%%rax"
        : "=a"(result) : "c"(index-1) : "%rdx");
      asm volatile("" ::: "memory");
    } while (buf->lock != seq);

    result += offset;
#else
    asm volatile("rdpmc;shlq $32,%%rdx;orq %%rdx,%%rax"
      : "=a"(result) : "c"(0) : "%rdx");
#endif

    result &= 0xffffffffffff;
    return result;
  }

  asm volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
    : "=a" (result) ::  "%rdx");
  return result;
}

u_int64_t t1, t2, sum;

#include "curve25519-yx-AVX512-def.h"
#include "curve25519-yx-AVX512-modp.h"
#include "curve25519-yx-AVX512-modp2.h"
#include "curve25519-yx-AVX512-pointgroup.h"
#include "curve25519-yx-AVX512-modL.h"
#include "curve25519-yx-AVX512-montgomery.h"
#include "curve25519-yx-basecombtable52.h"
// #include "basecombtableT.h"
#include "curve25519-yx-AVX512-comb.h"

typedef uint8_t curved25519[32];

void
ED25519_FN(clamp) (curved25519 out, const curved25519 in) {
	for (int i = 0; i < 32; i++) out[i] = in[i];
	out[0] &= 248;
	out[31] &= 127;
	out[31] |= 64;
}

void
ED25519_FN(x25519_scalarmult_mb8) (uint8_t out[8][32], const uint8_t secretkey[32], uint8_t point[8][32]) {
    curve25519_scalarmult_montgomery(out, secretkey, point);
}

void
ED25519_FN(ed25519_basemult_mb8) (uint8_t out[8][32], const uint8_t secretkey[32], uint8_t data[8][32]) {
    ed25519_basemult_comb(out, secretkey, data);
}

// void ED25519_FN(sc25519_invert32) (curved25519_key recip, curved25519_key s) {
// 	bignum256modm out, in;
// 	expand_raw256_modm(in, s);
// 	sc25519_invert(out, in);
// 	contract256_modm(recip, out);
// }