#include <stdio.h>
#include <stdint.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdbool.h>

static int powmod(int phi, int g, int exp, int N) {
    int retVal, i;

    retVal = phi % N;
    for (i = 0; i < exp; i++)
        retVal = (retVal * g) % N;
    return(retVal);
}

#ifdef FLTEMU_EMU // emulation header

#define MPFR_WANT_FLOAT128
#define mpfr_float128 __float128
#include <mpfr.h>

typedef union {
    __uint128_t b;
    __float128  v;
} quad;

typedef union {
    uint64_t b;
    double   v;
} dbl;

typedef struct {
    int16_t e;
    int64_t m;
} DOUBLE;

#define DECLARE_VAR(var) int16_t var##_0; int64_t var##_1

#define DECLARE_DATA_SE11M52(var, ref, s, e, m) DECLARE_VAR(var); _fuse_unpack_SE11M52(var##_0, var##_1, s, e, m)

#define LOAD(var, X, idx) var##_0 = (X)[idx].e; var##_1 = (X)[idx].m

#define LOAD_T(var, X, idx) var##_0 = (X##_0)[idx]; var##_1 = (X##_1)[idx]

#define LOAD_C(var, x) var##_0 = (x##_0); var##_1 = (x##_1)

#define STORE(Y, idx, var) (Y)[idx].e = var##_0; (Y)[idx].m = var##_1

#define STORE_T(Y, idx, var) (Y##_0)[idx] = var##_0; (Y##_1)[idx] = var##_1

#define ADD(c, a, b) _fuse_add(c##_0, c##_1, a##_0, a##_1, b##_0, b##_1)

#define SUB(c, a, b) _fuse_sub(c##_0, c##_1, a##_0, a##_1, b##_0, b##_1)

#define MUL(c, a, b) _fuse_mul(c##_0, c##_1, a##_0, a##_1, b##_0, b##_1)

#define DECLARE_ARRAY(var, n) int16_t var##_0[n]; int64_t var##_1[n]

static inline __float128 cospi(int64_t x, int64_t n)
{
    mpfr_t arg, result;
    mpfr_init2(arg, 256);
    mpfr_init2(result, 256);
    mpfr_set_si(arg, x, MPFR_RNDN);
    mpfr_div_si(arg, arg, n, MPFR_RNDN);
    mpfr_cospi(result, arg, MPFR_RNDN);
    __float128 out = mpfr_get_float128(result, MPFR_RNDN);
    mpfr_clear(arg);
    mpfr_clear(result);
    return out;
}

static inline __float128 sinpi(int64_t x, int64_t n)
{
    mpfr_t arg, result;
    mpfr_init2(arg, 256);
    mpfr_init2(result, 256);
    mpfr_set_si(arg, x, MPFR_RNDN);
    mpfr_div_si(arg, arg, n, MPFR_RNDN);
    mpfr_sinpi(result, arg, MPFR_RNDN);
    __float128 out = mpfr_get_float128(result, MPFR_RNDN);
    mpfr_clear(arg);
    mpfr_clear(result);
    return out;
}

static inline __float128 rcp_int(int64_t n)
{
    mpfr_t x;
    mpfr_init2(x, 256);
    mpfr_set_sj(x, n, MPFR_RNDN);
    mpfr_ui_div(x, 1, x, MPFR_RNDN);
    __float128 r = mpfr_get_float128(x, MPFR_RNDN);
    mpfr_clear(x);
    return r;
}

#define QUAD_to_BITS(val) (((quad){ .v = (val) }).b)

#define COSPI(var, x, n) _fuse_unpack128(var##_0, var##_1, QUAD_to_BITS(cospi(x, n)))

#define SINPI(var, x, n) _fuse_unpack128(var##_0, var##_1, QUAD_to_BITS(sinpi(x, n)))

#define ASSIGN_NTH(var, idx, val) var##_0[idx] = val##_0; var##_1[idx] = val##_1

#define ASSIGN(var, val) var##_0 = val##_0; var##_1 = val##_1

#define NEG(b, a) (b##_0) = (a##_0); (b##_1) = (0 - (a##_1))

#define MUL_BY2(b, a) (b##_0) = ((a##_0) + 1); (b##_1) = (a##_1)

#define DIV_BY2(b, a) (b##_0) = ((a##_0) - 1); (b##_1) = (a##_1)

#define DIV_BY4(b, a) (b##_0) = ((a##_0) - 2); (b##_1) = (a##_1)

#define DBL_to_BITS(val) (((dbl){ .v = (val) }).b)

#define RCP(out, in) dbl tmp; _fuse_pack(tmp.b, in.e, in.m); tmp.v = 1.0/tmp.v; _fuse_unpack(out##_0, out##_1, tmp.b) // dummy

#define RCP_INT(b, a) _fuse_unpack128(b##_0, b##_1, QUAD_to_BITS(rcp_int(a)))

#define ASSIGN_NTH_ZERO(var, idx) var##_0[idx] = 0x8000; var##_1[idx] = 0

#define ASSIGN_ZERO(var) var##_0 = 0x8000; var##_1 = 0

#define TCAST_DOUBLE(b, a) _fuse_unpack(b##_0, b##_1, DBL_to_BITS((double)(a)))


///////////////////////////// Arithmetic /////////////////////////////
// parts of VND
#define _sign_int65(_v) ((_Bool)((_v) >> 64))
static inline _Bool _carry_addq(__uint64_t a, __uint64_t b) { return _sign_int65((__uint128_t)a + (__uint128_t)b); }
static inline __uint64_t _addq(__uint64_t a, __uint64_t b) { return a + b; }
#define	_addqc(_c, _cflag, _a, _b) {(_cflag) = _carry_addq((_a), (_b)); (_c) = _addq((_a), (_b));}

#define _fuse_mul(_ce, _cm, _ae, _am, _be, _bm){ \
    /* if (am==0 || bm==0), we hardcode ouput exponent to 0. This can help improve addition's precision (the add will shift out improtant bits if we outptut a non-zero exponent and zero mantissa) */ \
    (_ce) = (((_am) == 0) || ((_bm) == 0)) ? 0x8000 : ((_ae)+(_be)+2); /* output exponent will be (ae + be + 2) to compensate the bias of output mantissa */ \
    (_cm) = ((int64_t)((((__int128_t)(int64_t)(_am)) * ((__int128_t)(int64_t)(_bm))) >> 64)); /* directly take the upper half of product as the output mantissa (this introduces 1/4 bias due to our way of constructing mantissa) */ \
}

#define _fuse_add(_ce, _cm, _ae, _am, _be, _bm){ \
    int16_t _diff0 = (_ae) - (_be); /* ae-be */ \
    int16_t _diff1 = (_be) - (_ae); /* be-ae */ \
    _Bool _flag = (_ae) < (_be); /* flag = (ae<be) */ \
    int16_t _shamt = _flag ? _diff1 : _diff0; /* shift amount (absolute value of diff) */ \
    int64_t _sm = _flag ? (_am) : (_bm); /* choose the mantissa to be shifted (the operand coresponding to the smaller exponent) */ \
    int64_t _km = _flag ? (_bm) : (_am); /* choose the mantissa to be kept (the operand coresponding to the larger exponent) */ \
    _sm = _sm >> _shamt; /* arithmetic right shift the mantissa selected to be shifted (km) */ \
    _sm = (_shamt < 64) ? _sm : 0; /* if the shift amount >= 64, then hardcode the result to be 0 */ \
    _Bool _cmc; int64_t _sum; /* _sum will be the output mantissa */ \
    /* add with carry for lower half: output lower half's sum in _sum; output lower half's carry in _cmc */ \
    /* _cmc = (_Bool)(((__uint128_t)_km + (__uint128_t)_sm) >> 64); */ \
    /* _sum = _km + _sm; */ \
    _addqc(_sum, _cmc, _km, _sm); \
    _Bool _ovf = (((_sm^_sum)&(_km^_sum)) < 0); /* detect sign overflow: overflow happens when the operands have the same sign and the sumb s sign is different */ \
    (_cm) = _ovf ? ((uint64_t)(((((__uint128_t)((uint64_t)_cmc)) << 64) | ((uint64_t)(_sum))) >> 1)) : _sum; /* if overflow, we want _cm = {_cmc,_sum>>1} (_cm will be biased by 1/2 due to the right shift).*/ \
    /* the output exponent will be the larger inputted exponent; if overflow, the output exponent needs to plus 1 to compensate the 1/2 bias we had on the output mantissa*/ \
    int16_t _larger_e = _flag ? (_be) : (_ae); \
    (_ce) = _ovf ? (_larger_e+1) : _larger_e; \
}

#define _fuse_sub(_ce, _cm, _ae, _am, _be, _bm){ \
    int64_t _bm_c = 0 - (_bm); \
    _fuse_add(_ce, _cm, _ae, _am, _be, _bm_c); \
}

#define _fuse_norm(_ce, _cm, _ae, _am){ \
    uint64_t msb = (_am) & 0x8000000000000000; \
    bool neg_mant = msb != 0; \
    int64_t _m = neg_mant ? (0 - (_am)) : (_am); \
    bool _all_0 = ((_m & 0x3fffffffffffffff) == 0); \
    unsigned _shamt = __builtin_clzll(_m) - 1u; \
    int64_t _m_shifted = (_m) << _shamt; \
    int16_t _ae_shifted = (_ae) - (uint16_t)_shamt;\
    (_ce) = _all_0 ? (_ae) : _ae_shifted; \
    (_m) = _all_0 ? (_m) : _m_shifted; \
    (_cm) = neg_mant ? (0 - (_m)) : (_m);\
}

#define _fuse_unpack(_ae, _am, _a){ \
    uint64_t (_sign) = ((_a) >> 63); \
    uint64_t (_mant) = ((_a) << 10); \
    (_mant) = ((_mant) & 0x3fffffffffffffff); \
    uint64_t (_exp) = ((_a) >> 52); \
    (_exp) = ((_exp) & 0x00000000000007ff); \
    (_mant) = ((_exp) != 0) ? ((_mant) | 0x4000000000000000): (_mant); /* if not denorm, add leading 1 */ \
    (_am) = ((_sign) == 1) ? (0 - (_mant)): (_mant); /* transfer sign to mant */ \
    (_ae) = ((int16_t)(_exp) - 1023); \
}

#define _fuse_unpack_SE11M52(_ae, _am, _s, _e, _m){ \
    uint64_t (_exp) = ((_e) & 0x07ff); \
    uint64_t (_mant) = ((_exp) != 0) ? ((_m) | 0x4000000000000000): (_m); /* if not denorm, add leading 1 */ \
    (_am) = ((_s) == -1) ? (0 - (_mant)): (_mant); /* transfer sign to mant */ \
    (_ae) = ((int16_t)(_exp) - 1023); \
}

#define _fuse_pack(_a, _ae, _am){ \
    uint64_t _msb = (_am) & 0x8000000000000000; \
    bool neg_mant = _msb != 0; \
    int64_t _m = neg_mant ? (0 - (_am)) : (_am); \
    bool _all_0 = ((_m & 0x7fffffffffffffff) == 0); \
    unsigned _shamt = __builtin_clzll(_m) - 1u; \
    int64_t _m_shifted = (_m) << _shamt; \
    int16_t _ae_shifted = (_ae) - (uint16_t)_shamt;\
    (_ae) = _all_0 ? (_ae) : _ae_shifted; \
    (_m) = _all_0 ? (_m) : _m_shifted; \
    uint64_t _exp = ((int64_t)(_ae)) + 1023; \
    bool _denorm = ((int64_t)(_exp)) <=  0; \
    (_exp) = _denorm ? 0 : (_exp); \
    uint64_t (_diff) = _denorm ? ((0-(_exp)) + 10) : 10; \
    (_m) = (_m) >> (_diff); \
    (_exp) = (_all_0) ? 0 : (_exp); /* 0 check */  \
    _Bool _ovf = ((_exp) >= 2047); /* ovf check */ \
    (_exp) = _ovf ? 2047 : (_exp); /* ovf check */ \
    (_m) = _ovf ? 0 : (_m); /* ovf check */ \
    (_m) = ((_m) & 0x000fffffffffffff); \
    (_exp) = ((_exp) << 52); \
    (_m) = ((_exp) | (_m)); \
    (_a) = ((_m) | (_msb)); \
}

#define _fuse_unpack128(_ae, _am, _a){ \
    uint64_t _sign = (_a) >> 127; \
    __uint128_t _am128 = (_a) << 14; \
    _am128 = _am128 & (((__uint128_t)1 << 126) - 1); \
    _am128 = _am128 >> 64; \
    __uint128_t _exp = ((_a) >> 112); \
    _exp = _exp & 0x07fff; \
    uint64_t _mant = (_exp != 0) ? (((uint64_t)_am128) | 0x4000000000000000) : ((uint64_t)_am128); /* if not denorm, add leading 1 */ \
    _mant = (_sign != 0) ? (0-_mant) : _mant; /* transfer sign to mant */ \
    (_ae) = (int16_t)(_exp - 16383); \
    (_am) = _mant;  \
}

#else // native FP header

#include <math.h>

#define DOUBLE double

#define DECLARE_VAR(var) double var
#define DECLARE_ARRAY(var, n) double var[n]
#define DECLARE_DATA_SE11M52(var, d, s, e, m) double (var) = (d) 

#define NEG(b, a)   (b) = -(a)
#define ADD(c, a, b) ((c) = (a) + (b))
#define SUB(c, a, b) ((c) = (a) - (b))
#define MUL(c, a, b) ((c) = (a) * (b))

#define MUL_BY2(b, a) (b) = ((a)*2.0)
#define DIV_BY2(b, a) (b) = ((a)/2.0)
#define DIV_BY4(b, a) (b) = ((a)/4.0)
#define RCP(b, a)   (b) = (1.0/(a))

#define LOAD(c, X, idx) ((c) = (X)[idx])
#define STORE(Y, idx, c) ((Y)[idx] = (c))

#define LOAD_T(c, X, idx) ((c) = (X)[idx])
#define STORE_T(Y, idx, c) ((Y)[idx] = (c))

#define LOAD_C(a, b) (a) = (b)

#define ASSIGN(a, b) (a) = (b)
#define ASSIGN_NTH(X, idx, c) ((X)[idx] = (c))
#define ASSIGN_ZERO(a) (a) = (0.0)
#define ASSIGN_NTH_ZERO(X, idx) ((X)[idx] = (0.0))

#define RCP_INT(a, b) (a) = (1.0 / ((double)(b)))
#define TCAST_DOUBLE(a, b) (a) = ((double)(b))

#define PI 3.14159265358979323846
#define COSPI(out, x, n) ((out) = cos(PI * (double)(x) / (double)(n)))
#define SINPI(out, x, n) ((out) = sin(PI * (double)(x) / (double)(n)))

#define FRAC_PI(out, x, n) ((out) = (PI * (double)(x) / (double)(n)))

#endif

