#include <stdint.h>
#include <stdio.h>

#include "cnn.h"

struct rng_state
{
    uint64_t a = 0x19980724cafebabeu;
    uint64_t b = 0x20010106feef1f0fu;
} rng;

static int16_t xorshift128p(rng_state* state)
{
    uint64_t t = state->a;
    uint64_t const s = state->b;
    state->a = s;
    t ^= t << 23; // a
    t ^= t >> 17; // b
    t ^= s ^ (s >> 26); // c
    state->b = t;
    return int16_t((t + s) >> 48);
}

#define RANDOMIZE_FLOAT_ARRAY(arr_name) \
for (size_t _i = 0; _i < (sizeof(arr_name) / sizeof(float)); ++_i) \
{ \
    float* flat_arr = reinterpret_cast<float*>(arr_name); \
    flat_arr[_i] = xorshift128p(&rng) * 1e-5f; \
}

volatile float anti_optimize;

int main()
{
    static float input[kFeatures][kInImSize][kInImSize];
    RANDOMIZE_FLOAT_ARRAY(input);
    static float weights[kFeatures][kFeatures][kKernel][kKernel];
    RANDOMIZE_FLOAT_ARRAY(weights);
    static float bias[kFeatures];
    RANDOMIZE_FLOAT_ARRAY(bias);
    static float output[kFeatures][kOutImSize][kOutImSize];

    for (size_t i_tile = 0; i_tile < kFeatures; i_tile += kITile) {
        printf("i_tile = %i (of %i)\n", int(i_tile), int(kFeatures));
        for (size_t h_tile = 0; h_tile < kOutImSize; h_tile += kHTile) {
            for (size_t w_tile = 0; w_tile < kOutImSize; w_tile += kWTile) {
                fill_output_tile(
                    input, weights, bias, output, i_tile, h_tile, w_tile);
                anti_optimize = output[i_tile][h_tile][w_tile];
            }
        }
    }
    printf("%f\n", anti_optimize);
}
