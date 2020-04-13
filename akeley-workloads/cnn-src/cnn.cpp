// High ILP test case. Translated from the CS 133 OpenCL Convolutional
// Neural Network assignment, which uses aggressive vectorization and
// loop transformations.
//
// Note: Due to time constraints, I haven't really tried to verify
// that this translation is correct. Don't use it for your
// self-driving car project.
#include "cnn.h"

void fill_output_tile(
    const float (*I) [kInImSize][kInImSize],
    const float (*W) [kFeatures][kKernel][kKernel],
    const float *B,
    float (*O) [kOutImSize][kOutImSize],
    size_t i_tile, size_t h_tile, size_t w_tile)
{
    const size_t i0 = i_tile + 0;
    const size_t i1 = i_tile + 1;
    const size_t i2 = i_tile + 2;
    const size_t i3 = i_tile + 3;
    const float bias0 = B[i0];
    const float bias1 = B[i1];
    const float bias2 = B[i2];
    const float bias3 = B[i3];
    float16 sum0_0 = (float16) (bias0);
    float16 sum0_1 = (float16) (bias0);
    float16 sum1_0 = (float16) (bias1);
    float16 sum1_1 = (float16) (bias1);
    float16 sum2_0 = (float16) (bias2);
    float16 sum2_1 = (float16) (bias2);
    float16 sum3_0 = (float16) (bias3);
    float16 sum3_1 = (float16) (bias3);
    for (size_t j = 0; j < kFeatures; ++j) {
        const float (*I_layer)[kInImSize] = I[j];
        const float (*W_layer0)[kKernel] = W[i0][j];
        const float (*W_layer1)[kKernel] = W[i1][j];
        const float (*W_layer2)[kKernel] = W[i2][j];
        const float (*W_layer3)[kKernel] = W[i3][j];
        for (size_t p = 0; p < kKernel; ++p) {
            for (size_t q = 0; q < kKernel; ++q) {
                const float weight0 = W_layer0[p][q];
                const float weight1 = W_layer1[p][q];
                const float weight2 = W_layer2[p][q];
                const float weight3 = W_layer3[p][q];
                const size_t hp0 = h_tile * 2 + 0 + p;
                const size_t hp1 = h_tile * 2 + 1 + p;
                const size_t wq = w_tile * 2 + q;
                const float16 img_row0 = vload16(0, &I_layer[hp0][wq]);
                const float16 img_row1 = vload16(0, &I_layer[hp1][wq]);
                sum0_0 += weight0 * img_row0;
                sum1_0 += weight1 * img_row0;
                sum2_0 += weight2 * img_row0;
                sum3_0 += weight3 * img_row0;
                sum0_1 += weight0 * img_row1;
                sum1_1 += weight1 * img_row1;
                sum2_1 += weight2 * img_row1;
                sum3_1 += weight3 * img_row1;
            }
        }
    }

    size_t h;
    const float8 zero = (float8) (0.0);
    float16 vmax0;
    float16 vmax1;
    float16 vmax2;
    float16 vmax3;
    float8 hmax0, out_row0;
    float8 hmax1, out_row1;
    float8 hmax2, out_row2;
    float8 hmax3, out_row3;
    h = h_tile + 0;
    vmax0 = max(sum0_0, sum0_1);
    vmax1 = max(sum1_0, sum1_1);
    vmax2 = max(sum2_0, sum2_1);
    vmax3 = max(sum3_0, sum3_1);
    hmax0 = max(vmax0.s02468ace(), vmax0.s13579bdf());
    hmax1 = max(vmax1.s02468ace(), vmax1.s13579bdf());
    hmax2 = max(vmax2.s02468ace(), vmax2.s13579bdf());
    hmax3 = max(vmax3.s02468ace(), vmax3.s13579bdf());
    out_row0 = max(zero, hmax0);
    out_row1 = max(zero, hmax1);
    out_row2 = max(zero, hmax2);
    out_row3 = max(zero, hmax3);
    vstore8(out_row0, 0, &O[i0][h][w_tile]);
    vstore8(out_row1, 0, &O[i1][h][w_tile]);
    vstore8(out_row2, 0, &O[i2][h][w_tile]);
    vstore8(out_row3, 0, &O[i3][h][w_tile]);
}
