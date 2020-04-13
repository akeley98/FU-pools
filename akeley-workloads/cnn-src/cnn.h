#include <stddef.h>
#include <xmmintrin.h>

const size_t kFeatures = 48;
const size_t kKernel = 5;
const size_t kImSize = 64;
const size_t kInImSize = kImSize + kKernel - 1;
const size_t kOutImSize = kImSize / 2;

// Hard wired
const size_t kITile = 4;
const size_t kHTile = 1;
const size_t kWTile = 8;

static_assert(kImSize % kOutImSize == 0);
static_assert(kOutImSize % kWTile == 0);
static_assert(kOutImSize % kHTile == 0);
static_assert(kFeatures % kITile == 0);

void fill_output_tile(
    const float (*I) [kInImSize][kInImSize],
    const float (*W) [kFeatures][kKernel][kKernel],
    const float *B,
    float (*O) [kOutImSize][kOutImSize],
    size_t i_tile, size_t h_tile, size_t w_tile);

struct float8
{
    __m128 a, b;

    float8 () = default;

    float8 (float s)
    {
        a = _mm_set1_ps(s);
        b = _mm_set1_ps(s);
    }
};

inline float8 max(float8 lhs, float8 rhs)
{
    float8 result;
    result.a = _mm_max_ps(lhs.a, rhs.a);
    result.b = _mm_max_ps(lhs.b, rhs.b);
    return result;
}

inline void vstore8(float8 data, size_t idx, float* ptr)
{
    ptr += idx * 8;
    _mm_storeu_ps(&ptr[0], data.a);
    _mm_storeu_ps(&ptr[4], data.b);
}

struct float16
{
    __m128 a, b, c, d;

    float16 () = default;

    float16 (float s)
    {
        a = _mm_set1_ps(s);
        b = _mm_set1_ps(s);
        c = _mm_set1_ps(s);
        d = _mm_set1_ps(s);
    }

    float16& operator += (float16 o)
    {
        a = _mm_add_ps(a, o.a);
        b = _mm_add_ps(b, o.b);
        c = _mm_add_ps(c, o.c);
        d = _mm_add_ps(d, o.d);
        return *this;
    }

    float8 s02468ace() const
    {
        float8 result;
        result.a = _mm_shuffle_ps(a, b, 0x88);
        result.b = _mm_shuffle_ps(c, d, 0x88);
        return result;
    }

    float8 s13579bdf() const
    {
        float8 result;
        result.a = _mm_shuffle_ps(a, b, 0xDD);
        result.b = _mm_shuffle_ps(c, d, 0xDD);
        return result;
    }
};

inline float16 operator* (float f, float16 v)
{
    __m128 fv = _mm_set1_ps(f);
    float16 result;
    result.a = _mm_mul_ps(fv, v.a);
    result.b = _mm_mul_ps(fv, v.b);
    result.c = _mm_mul_ps(fv, v.c);
    result.d = _mm_mul_ps(fv, v.d);
    return result;
}

inline float16 vload16(size_t idx, const float* ptr)
{
    ptr += idx * 16;
    float16 result;
    result.a = _mm_loadu_ps(&ptr[0]);
    result.b = _mm_loadu_ps(&ptr[4]);
    result.c = _mm_loadu_ps(&ptr[8]);
    result.d = _mm_loadu_ps(&ptr[12]);
    return result;
}

inline float16 max(float16 lhs, float16 rhs)
{
    float16 result;
    result.a = _mm_max_ps(lhs.a, rhs.a);
    result.b = _mm_max_ps(lhs.b, rhs.b);
    result.c = _mm_max_ps(lhs.c, rhs.c);
    result.d = _mm_max_ps(lhs.d, rhs.d);
    return result;
}
