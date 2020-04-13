// Workload based on "real world" game code. This is some code that I
// wrote as part of a Unity plugin that produces particle-based
// visualizations of tumors. I think it's a reasonable benchmark for
// "typical" code that was written with real-world time constraints,
// isn't all that optimized for performance, and does significant
// dynamic memory allocation.
//
// Within if statements, I had to replace float compares with integer
// compares, because gem5 does not implement the 'fcomi' instruction.
// I still get that warning (not sure from where), but at least I should
// have minimized the effect that missing 'fcomi' has on validity.

#include <assert.h>
#include <stdio.h>
#include <vector>
#include <stddef.h>
#include "bedrock_particles.h"

constexpr size_t particle_system_count = 4;
constexpr size_t step_count = 3;
constexpr size_t interpolate_count = 10;
constexpr float tumor_cells = 3000e6f;
constexpr float ctl_cells = 150e6f;
constexpr float treg_cells = 90e6f;
constexpr float t_helper_cells = 180e6f;

int main()
{
    std::vector<h_particles_snapshot> particle_systems;
    std::vector<visual_particle> vp_vector(4000);

    for (size_t i = 0; i < particle_system_count; ++i) {
        particle_systems.push_back(hpsBedrock_new_particles_snapshot());
    }

    for (size_t step = 0; step < step_count; ++step) {
        for (size_t i = 0; i < particle_system_count; ++i) {
            printf("%i/%i %i/%i\n", int(step), int(step_count),
                int(i), int(particle_system_count));
            Bedrock_simulation_step_fn_inplace(
                particle_systems[i],
                tumor_cells,
                ctl_cells,
                t_helper_cells,
                treg_cells);
            for (size_t interp = 0; interp < interpolate_count; ++interp) {
                return_code rc = Bedrock_interpolate_visual_particles(
                    particle_systems[i],
                    1.0f / interpolate_count * interp,
                    vp_vector.data(),
                    vp_vector.size() * sizeof *vp_vector.data());
                assert(rc >= 0);
            }
        }
    }
}
