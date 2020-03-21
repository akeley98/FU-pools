#ifndef BEDROCK_PARTICLES_H_
#define BEDROCK_PARTICLES_H_
// Bedrock 4.0 particle simulation prototype declarations.
// 
// By David Zhao Akeley (dza724`at`gmail.com) 2019-09-23
// 
// A particle simulation keeps track of a cloud of particles that
// visually represents a tumor and its changes through time. It is
// mainly designed for its aesthetic value: it is not a real
// "simulation" in that it's not involved in actually determining
// tumor behaviour.
// 
//     Particle Simulation Goals:
// 
// * Run multiple simulations concurrently using the Unity Job System.
// 
// * Visually represent tumor, CTL, T helper, and Treg cells as
//   particles representing (on the order of) 1 million cells each.
// 
// * Accept target cell counts in real time (of the 4 types mentioned
//   above) and modify the simulation in a visually plausible way to
//   represent cells added to/lost from the tumor. For example, if the
//   tumor cell count increases, then some tumor cell particles will
//   split in half.
// 
// * Maintain a cloud of particles of mostly spherical shape and
//   constant density (cells/cubic unit). This is done by computing a
//   voxel grid that keeps track of cell densities in space, and by
//   stochastically repositioning particles using a cost function that
//   incentivizes moving inward towards the origin and outwards of
//   overfilled voxels.
//   
// 
//     Particle Simulation Concepts:
// 
       struct visual_particle;
//
// Structure of 7 floats used to position and color one particle of a
// simulation on screen. A particle simulation can be queried every
// frame for an array of visual_particle structs to draw on screen.
// This struct is kept deliberately simple since it is used to pass
// data between the C++ and Unity C# parts of the simulation.
// 
       struct particles_snapshot;
// 
// The particle simulation splits time into discrete simulation steps
// (that are much longer than a single game frame). A
// particles_snapshot holds the state that some particle simulation
// has at a certain simulation step, along with the state it had at
// recent earlier steps. The simulation step function computes the
// particles_snapshot for step N+1 given the stapshot at step N.
// 
// To facilitate multi-threading, at a high level we prefer to program
// functionally and create modified copies of objects rather than
// modify them in place. The name particles_snapshot (as opposed to
// particles_state, etc.) reinforces this view.
// 
       struct sim_particle;
// 
// "Simulation particle" used in particles_snapshot to store a
// particle's position, cell count, and history. Distinct from
// visual_particle, which is a much simpler structure.
// 
// Note that a single sim_particle may produce multiple
// visual_particles, because CTL and T helper particles are viewed
// as attachments on a tumor sim_particle.
// 
//     Visual Particle Interpolation:
// 
// We need to interpolate between particle positions when drawing a
// particle simulation in the time between simulation steps. Because
// particles_snapshot stores recent history, this "interpolation" can
// be done using just one particles_snapshot (1) -- rendering in the
// time between steps N-1 and N can be done by interpolating just the
// snapshot data for step N. (2) Conceptually, visual particle
// interpolation is the process that converts simulation particles
// into visual particles.
// 
// (1) This complicates the design of a simulation particle, but
// simplifies multi-threading by minimizing dependencies.
// 
// (2) In reality, there is some lag when rendering (not least of
// which is because we would have to see the future to have the
// particles_snapshot for step N be visible at step N-1), but this
// just complicates this conceptual discussion.
// 
#include <stdint.h>
       typedef int64_t h_particles_snapshot; // [hps]
// 
// 64-bit integer handle for a particles_snapshot (long in C#). All
// handles given will be strictly positive. Can be converted to a
// pointer to heap-allocated particles_snapshot. Since these are
// really typeless integers, use the hps Hungarian notation prefix to
// distinguish them.
// 
// h_particles_snapshot is used to avoid handing out raw pointers to
// C# programmers (I don't trust them). It also gives us the chance to
// run some additional checks.
// 
//     Error handling strategy:
// 
// The error handling for this library is really primitive. If any
// C#-visible function returns a negative number, then it has failed,
// with the return value being the negative line number where the
// failure was detected.

typedef h_particles_snapshot return_code;

#include <stddef.h>

// Structure for holding the position, color, and radius of an
// on-screen particle. These are produced by interpolating the state
// between simulation steps.
// 
// DO NOT add, remove, reorder, or otherwise mess with this
// struct. This is/will be used to communicate between C++, C#, Unity,
// and OpenGL and assumptions about the memory layout and size of
// this class are all over the place in this code.
struct visual_particle
{
    // Formerly struct particle_render_info
    float x = 0;
    float y = 0;
    float z = 0;
    float red = 0;
    float green = 0;
    float blue = 0;
    float radius = 0.0f;
};

extern "C" {

// Heap-allocate a new particle simulation snapshot, and return the
// handle to it.
h_particles_snapshot hpsBedrock_new_particles_snapshot();

// Deallocate the particle simulation with the given handle.
return_code Bedrock_free_hps(h_particles_snapshot);

// Make a deep copy of the given particles snapshot, and return a
// handle to the copy.
h_particles_snapshot hpsBedrock_copy_particles_snapshot(h_particles_snapshot);

// Given a snapshot of a particle simulation at step N, compute its
// snapshot for step N+1 using the given target cell counts, and write
// the new state to given particles snapshot in-place. As a convenience,
// if successful, this function returns the passed-in handle.
// 
// Long rationale: The idea is that the given particles snapshot was
// just created by hpsBedrock_copy_particles_snapshot. Conceptually, I
// wanted one function that took a read-only particles_snapshot and
// created an entirely new particles_snapshot with its state
// advanced. However, this was not possible, as only the thread that
// schedules a simulation step can safely copy the particles_snapshot,
// while only the actual worker thread can advance the snapshot's
// state.
return_code Bedrock_simulation_step_fn_inplace(
    h_particles_snapshot hps,
    float target_tumor_cells,
    float target_ctl_cells,
    float target_t_helper_cells,
    float target_treg_cells);

// Perform visual particle interpolation on the given
// particles_snapshot.  t=0 indicates that we want the previous step's
// data, and t=1 indicates current step's data. (t is clamped to [0,1]
// in case it is outside this range). The visual_particles are written
// to the passed array pointer, and no more than max_bytes are
// written.
// 
// Returns the number of visual_particles written, unless the passed
// array was too small, in which case a negative error code is
// returned. All other negative numbers indicate an error as usual.
return_code Bedrock_interpolate_visual_particles(
    h_particles_snapshot hps,
    float t,
    visual_particle* vp_ptr,
    size_t max_bytes);

// Calculate the max_bytes parameter needed to successfully store the
// result of Bedrock_interpolate_visual_particles, assuming the
// particles snapshot named by hps doesn't change in the meantime.
return_code Bedrock_interpolate_visual_particles_bytes_needed(
    h_particles_snapshot hps);

// Get the number of valid (allocated) particles_snapshots. Mostly
// used to check for memory leaks.
return_code Bedrock_get_hps_count();

// NOT threadsafe.
// 
// Force deleting all particles_snapshots and their handles.
// Obviously VERY dangerous.  This is needed because we can't do
// proper cleanup when we exit play mode in the editor, which can
// quickly leak away GBs of memory.
return_code Bedrock_hack_delete_all_particles_snapshots();

// *** Setters for the "constants" ***
// 
// These are here to avoid having to recompile to tweak a magic number.
// 
// NOTE: These technically are not thread safe, but, the consequences
// of a data race on them is not so dire. I expect these setters will
// be called just once as part of setup, if at all. (I just expose
// them here in case someone wants to muck with the magic numbers w/o
// recompiling the C++ library).
void Bedrock_set_particle_cells_per_cubic_unit(float arg);
void Bedrock_set_max_cells_per_voxel(float arg);
void Bedrock_set_soft_max_cells_per_voxel(float arg);
void Bedrock_set_max_cells_per_particle(float arg);
void Bedrock_set_min_cells_per_particle(float arg);
void Bedrock_set_undersize_removal_ratio(int32_t arg);
void Bedrock_set_ctl_death_probability(float arg);
void Bedrock_set_ctl_particle_cell_ratio(float arg);
void Bedrock_set_ctl_particle_distance_ratio(float arg);
void Bedrock_set_t_helper_particle_cell_ratio(float arg);
void Bedrock_set_t_helper_particle_distance_ratio(float arg);
void Bedrock_set_drawn_radius_scale_factor(float arg);
void Bedrock_set_particle_init_volume_scale_magic(float arg);
void Bedrock_set_cost_distance_from_origin_magic(float arg);
void Bedrock_set_desperation_distance_from_origin_magic(float arg);
void Bedrock_set_stochastic_vector_step_x(float arg);
void Bedrock_set_stochastic_vector_step_y(float arg);
void Bedrock_set_stochastic_vector_step_z(float arg);
void Bedrock_set_stochastic_tries(int32_t arg);
void Bedrock_set_stochastic_step_size_step(float arg);
void Bedrock_set_treg_hack_remove_treg_probability(double arg);
void Bedrock_set_tumor_red(float);
void Bedrock_set_tumor_green(float);
void Bedrock_set_tumor_blue(float);
void Bedrock_set_ctl_red(float);
void Bedrock_set_ctl_green(float);
void Bedrock_set_ctl_blue(float);
void Bedrock_set_t_helper_red(float);
void Bedrock_set_t_helper_green(float);
void Bedrock_set_t_helper_blue(float);
void Bedrock_set_treg_red(float);
void Bedrock_set_treg_green(float);
void Bedrock_set_treg_blue(float);

} // end extern "C"

#endif
