// Bedrock 4.0 particle simulation prototype C++11 implementation.
// Originally developed on GNU Linux; port to Windows, iOS, Android et
// freaking al. later.
// 
// By David Zhao Akeley (dza724`at`gmail.com) 2019-09-23
// 
// See bedrock_particles.hpp for discussion. PLEASE PLEASE skim & and
// keep current the documentation there if you make changes to this file.
// 
// PROPAGANDA: When I wrote this in C++ instead of C#, I knew that
// this would make the Bedrock program harder to port, risking
// reducing platform support. As a Linux user, my conscience requires
// me to request/beg that you continue providing this plugin on Linux
// builds, or at least disable the particle simulation on Linux
// instead of letting this plugin prevent support for that platform
// entirely (I believe there is a preprocessor variable for this in
// ParticleSimulation.cs -- this may also be useful for devices that
// are not powerful enough to run this code).
#include <algorithm>
#include <limits>
#include <math.h>
#include <memory>
#include <mutex>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
// needed for memcpy
#include <string.h>
#include <unordered_map>
#include <utility>
using std::move;
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

#include "bedrock_particles.h"

// *** Random math things that need to be at the top of the file. ***

// To reduce dependencies, I decided to switch from glm to
// implementing my own super-simple vectors. I only implement the
// operators I actually use.
struct Bedrock_vec3
{
    float x;
    float y;
    float z;
    
    Bedrock_vec3(float x_=0, float y_=0, float z_=0)
    {
        x = x_;
        y = y_;
        z = z_;
    }
    bool operator==(Bedrock_vec3 o) const
    {
        return x == o.x && y == o.y && z == o.z;
    }
};
using vec3 = Bedrock_vec3;

struct Bedrock_ivec3
{
    int32_t x;
    int32_t y;
    int32_t z;
    
    Bedrock_ivec3(int32_t x_=0, int32_t y_=0, int32_t z_=0)
    {
        x = x_;
        y = y_;
        z = z_;
    }
    operator Bedrock_vec3() const
    {
        return Bedrock_vec3(x, y, z);
    }
    bool operator==(Bedrock_ivec3 o) const
    {
        return x == o.x && y == o.y && z == o.z;
    }
};
using ivec3 = Bedrock_ivec3;

inline vec3 operator*(vec3 v, float s)
{
    return vec3(v.x * s, v.y * s, v.z * s);
}
inline vec3 operator*(float s, vec3 v) { return v * s; }
inline vec3& operator*=(vec3& v, float s) { v = v * s; return v; }

inline vec3 operator+(vec3 a, vec3 b)
{
    return vec3 (a.x + b.x, a.y + b.y, a.z + b.z);
}
inline vec3& operator+=(vec3& a, vec3 b) { a = a + b; return a; }

inline vec3 operator-(vec3 a, vec3 b)
{
    return vec3 (a.x - b.x, a.y - b.y, a.z - b.z);
}
inline vec3& operator-=(vec3& a, vec3 b) { a = a - b; return a; }

inline float clamp(float x, float low, float high)
{
    return std::max(low, std::min(high, x));
}

inline float cube_root(float n)
{
    return powf(n, 1.f/3.f);
}

// Return false iff any args are inf/NaN. NOTE: These checks will NOT WORK
// with unsafe compiler flags like -Ofast or -funsafe-math-optimizations !!!
//
// Disabled for gem5
inline bool is_real(vec3)
{
    // return v.x - v.x == 0 && 
    //        v.y - v.y == 0 && 
    //        v.z - v.z == 0;
    return true;
}

inline bool is_real(float)
{
    // return f - f == 0;
    return true;
}



// *** Exception used for my primitive line-number error reporting scheme. ***
struct particles_borked
{
    return_code error_code;
};

// Ersatz assert. If expr is false, throw an exception with the error
// code being the negative line number.
#define CHECK(expr) do { if (!(expr)) { \
    particles_borked CHECK_; \
    CHECK_.error_code = -__LINE__; \
    throw CHECK_; \
}} while (0);



// *** Particle simulation "constants" ***

// NOTE: If you are here to tweak the constants and recompile, you
// don't have to do it! Instead, you can use the Bedrock_set_[magic
// number name] functions. (Unless they don't work, in which case you
// should call/email and yell at me).

// Used to determine the radius of particles (modelled as axis-aligned
// cubes even though they're rendered as balls). Their 'radius' (half
// of cube edge length) is calculated such that the density of cells
// (cells / cubic unit) is this value.
static float particle_cells_per_cubic_unit = 6e6f;

// Target max cells in a voxel. If there are too many cells in a voxel
// (calculated by summing up cells of particles within the voxel),
// then the simulation attempts to move the particles out of the voxel.
static float max_cells_per_voxel = 3e6f;

// Similar to max_cells_per_voxel, but used when considering moving a
// particle back in towards the origin. Using a slightly lower limit
// avoids having us compressing particles together too aggressively,
// which caused jittery behaviour.
static float soft_max_cells_per_voxel = 2.5e6f;

// Maximum number of cells represented by a single particle. Particles
// are split up to respect this limit.
static float max_cells_per_particle = 2e6f;

// Minimum number of cells per particle. This is one-half of the
// maximum cells so that newly split particles are valid (if higher
// than 1/2, the 2 new cells split from an oversized particle would be
// below the size limit and immediately destroyed).
static float min_cells_per_particle = max_cells_per_particle * 0.5f;

// (Roughly) 1-in-undersize_removal_ratio undersize cells will be
// removed per simulation step.
static int32_t undersize_removal_ratio = 3;

// Probability per simulation step that a CTL particle attached to a
// tumor particle will die, destroying the connected tumor particle
static float ctl_death_probability = 0.15f;

// Ratio of (ctl cells represented by a ctl particle) / (tumor cells
// represented by particle the ctl is visually attached to)
static float ctl_particle_cell_ratio = 0.375f;
static float ctl_particle_radius_ratio = cube_root(ctl_particle_cell_ratio);

// Ratio of (tumor particle center -> atttached ctl particle center distance)
// / (tumor particle visual radius -- see drawn_radius_scale_factor).
static float ctl_particle_distance_ratio = 1.5f;

// Same for T Helpers.
static float t_helper_particle_cell_ratio = 0.45f;
static float t_helper_particle_radius_ratio = cube_root(t_helper_particle_cell_ratio);
static float t_helper_particle_distance_ratio = 2.4f;

// Drawn particle radius is scaled by this factor from the actual
// particle radius (as stored in the simulation)
static float drawn_radius_scale_factor = 0.6f;

// Weird magic numbers  
static float particle_init_volume_scale_magic = 1.2f;
static float cost_distance_from_origin_magic = 2.5e5f;
static float desperation_distance_from_origin_magic = 3.0f;
static vec3 stochastic_vector_step = vec3(7, 11, 13) * 0.0081f;

// Number of random displacement vectors considered per particle while
// trying to move a particle to minimize cost (overcrowding & distance
// from origin).
static int32_t stochastic_tries = 4;

// Related (but not equal) to the magnitude of the displacement
// vectors considered while moving a particle to minimize cost. At
// time of writing, this step size is modified based on a particle's
// desperation.
static float stochastic_step_size_step = 0.2f;

// TREG HACK
// 
// Late in development we noticed some really aesthetically
// unfortunate emergent behaviour: TReg particles tend to settle
// downward and clump together VERY conspicuously in the core of the
// particle simulation. We did not have time to investigate the real
// causes of this "brazil nut effect" behaviour. The only solution we
// could implement for now was to give each TReg particle a small
// chance of being destroyed each simulation step, triggering creation
// of new TReg particles that would then be evenly distributed.
static double treg_hack_remove_treg_probability = 0.025;

// Colors used to create visual particles of tumor/ctl/t_helper/treg type.
static float tumor_red = 0.0f;
static float tumor_green = 0.4f;
static float tumor_blue = 0.9f;
static float ctl_red = 1.0f;
static float ctl_green = 0.0f;
static float ctl_blue = 0.4f;
static float t_helper_red = 208/255.f;
static float t_helper_green = 240/255.f;
static float t_helper_blue = 192/255.f;
static float treg_red = 0.0f;
static float treg_green = 168/255.f;
static float treg_blue = 0.0f;



// *** Setters for the "constants" ***
// 
// These are here to avoid having to recompile to tweak a magic number.
// 
// NOTE: These technically are not thread safe, but, the consequences
// of a data race on them is not so dire. I expect these setters will
// be called just once as part of setup, if at all. (I just expose
// them here in case someone wants to muck with the magic numbers w/o
// recompiling the C++ library).

#define BEDROCK_DEFINE_SETTER(MAGIC_NUMBER_NAME) \
extern "C" void Bedrock_set_##MAGIC_NUMBER_NAME(decltype(MAGIC_NUMBER_NAME) arg) \
{ MAGIC_NUMBER_NAME = arg; }

BEDROCK_DEFINE_SETTER(particle_cells_per_cubic_unit);
BEDROCK_DEFINE_SETTER(max_cells_per_voxel);
BEDROCK_DEFINE_SETTER(soft_max_cells_per_voxel);
BEDROCK_DEFINE_SETTER(max_cells_per_particle);
BEDROCK_DEFINE_SETTER(min_cells_per_particle);
BEDROCK_DEFINE_SETTER(undersize_removal_ratio);
BEDROCK_DEFINE_SETTER(ctl_death_probability);
extern "C" void Bedrock_set_ctl_particle_cell_ratio(float arg)
{
    ctl_particle_cell_ratio = arg;
    ctl_particle_radius_ratio = cube_root(arg);
}
BEDROCK_DEFINE_SETTER(ctl_particle_distance_ratio);
extern "C" void Bedrock_set_t_helper_particle_radius_ratio(float arg)
{
    t_helper_particle_cell_ratio = arg;
    t_helper_particle_radius_ratio = cube_root(arg);
}
BEDROCK_DEFINE_SETTER(drawn_radius_scale_factor);
BEDROCK_DEFINE_SETTER(particle_init_volume_scale_magic);
BEDROCK_DEFINE_SETTER(cost_distance_from_origin_magic);
extern "C" void Bedrock_set_stochastic_vector_step_x(float x)
{
    stochastic_vector_step.x = x;
}
extern "C" void Bedrock_set_stochastic_vector_step_y(float y)
{
    stochastic_vector_step.y = y;
}
extern "C" void Bedrock_set_stochastic_vector_step_z(float z)
{
    stochastic_vector_step.z = z;
}
BEDROCK_DEFINE_SETTER(stochastic_tries);
BEDROCK_DEFINE_SETTER(stochastic_step_size_step);
BEDROCK_DEFINE_SETTER(treg_hack_remove_treg_probability);
BEDROCK_DEFINE_SETTER(tumor_red);
BEDROCK_DEFINE_SETTER(tumor_green);
BEDROCK_DEFINE_SETTER(tumor_blue);
BEDROCK_DEFINE_SETTER(ctl_red);
BEDROCK_DEFINE_SETTER(ctl_green);
BEDROCK_DEFINE_SETTER(ctl_blue);
BEDROCK_DEFINE_SETTER(t_helper_red);
BEDROCK_DEFINE_SETTER(t_helper_green);
BEDROCK_DEFINE_SETTER(t_helper_blue);
BEDROCK_DEFINE_SETTER(treg_red);
BEDROCK_DEFINE_SETTER(treg_green);
BEDROCK_DEFINE_SETTER(treg_blue);



// *** Structure definitions ***

// Helper struct for interpolate_display_position (for efficiency,
// avoid recalcalating bezier coefficients on every call).
struct bezier_coefficients
{
    float t, s0, s1, s2, s3;
    
    bezier_coefficients(float t_)
    {
        CHECK(0 <= t_ && t_ <= 1);
        t = t_;
        float t1 = 1-t;
        s0 = t1 * t1 * t1;
        s1 = (3 * t1) * (t1 * t);
        s2 = (3 * t1) * (t * t);
        s3 = t * t * t;
    }
};

// Main structure for representing a simulation particle, which itself
// represents a bundle of tumor cells, or a bundle of T reg cells.
struct sim_particle
{
    // Formerly struct particle.
    
    // Used as temporary storage for voxel fill functions (records
    // distribution of this particle's cells in the voxel grid).
    // See VOXEL FILL MODEL.
    struct voxel_fill;
    using voxel_fill_list = std::vector<voxel_fill>;
    voxel_fill_list tmp_voxel_fill_list;
    
    // Number of cells represented by this particle. Need not be integer.
    float cell_count = 0;
    
    // Current xyz position in this simulation step.
    vec3 position = vec3(0, 0, 0);
    
    // Xyz position in the previous simulation step. Used to smoothly
    // transition particle positions between simulation steps for
    // display purposes. Update with push_new_position.
    vec3 old_position = vec3(0, 0, 0);
    
    // Xyz position 2 simulation steps ago.
    vec3 even_older_position = vec3(0, 0, 0);
    
    // Number of cells in previous simulation step for this particle.
    // Update with push_cell_count.
    float old_cell_count = 0;
    
    // Derived from radii varibles with update_radius.
    float radius = 0;
    float old_radius = 0;
    
    // Iff has_ctl is true, then an extra ctl particle is drawn
    // attached to this particle, with the ctl's position offset from
    // the particle position by some scaled
    // ctl_particle_normalized_vector (which should have magnitude 1
    // as the name implies).
    vec3 ctl_particle_normalized_offset = vec3(0, 0, 0);
    bool has_ctl = false;
    
    // Iff true, a t helper is drawn attached to the particle, in the
    // same direction as ctl_particle_normalized_offset.
    bool has_t_helper = false;
    
    // Increases for every simulation step where the particle is
    // unhappy (cost > 0).
    uint32_t desperation = 0;
    
    
    // Use the passed arg as the new position, and populate the older
    // position fields as appropriate.
    void push_new_position(vec3 arg)
    {
        even_older_position = old_position;
        old_position = position;
        position = arg;
    }
    
    void push_cell_count(float arg)
    {
        old_cell_count = cell_count;
        cell_count = arg;
    }
    
    // Calculate the radius and old_radius based on the cell count.
    // See particle_cells_per_cubic_unit.
    void update_radius()
    {
        auto volume = cell_count / particle_cells_per_cubic_unit;
        radius = cube_root(volume) * 0.5f;
        
        volume = old_cell_count / particle_cells_per_cubic_unit;
        old_radius = cube_root(volume) * 0.5f;
    }
    
    // Given a bezier_coefficents struct constructed from an
    // interpolation parameter t (0 through 1), write the actual
    // displayed position of this particle to *particle_position.  t=0
    // indicates displayed position is the old position; t=1 indicates
    // new position. Radius is also interpolated, and written through
    // *particle_radius.
    void interpolate_display_position(
        bezier_coefficients bz,
        visual_particle* ri) const
    {
        vec3 p0 = old_position;
        vec3 p1 = old_position + (old_position - even_older_position);
        vec3 p2 = old_position;
        vec3 p3 = position;
        
        ri->radius = (bz.t * radius + (1 - bz.t) * old_radius)
                   * drawn_radius_scale_factor;
        
        vec3 pos = (bz.s0 * p0 + bz.s1 * p1) + (bz.s2 * p2 + bz.s3 * p3);
        ri->x = pos.x;
        ri->y = pos.y;
        ri->z = pos.z;
    }
    
    struct voxel_fill
    {
        // Integer 3-vector that identifies the voxel to position cells in. 
        ivec3 voxel_id_ivec;
        
        // Number of cells of this particle that is within the voxel named
        // by voxel_id_ivec.
        float cell_count_added;
        
        // Need to be paranoid and make sure the move constructors
        // etc. are noexcept to avoid unneeded copies on vector
        // reallocation. (If this doesn't make sense; don't worry. You
        // may want to read about std::move_if_noexcept.)
        voxel_fill(const voxel_fill&) noexcept = default;
        voxel_fill(voxel_fill&&) noexcept = default;
        voxel_fill& operator= (const voxel_fill&) noexcept = default;
        voxel_fill& operator= (voxel_fill&&) noexcept = default;
        ~voxel_fill() noexcept = default;
    };
};

// Arrays of particles that constitute a particle system. Also
// includes the voxel stochastic vector and random number generator,
// which are re-used in between states but otherwise don't correspond
// to any user-visible state.
struct particles_snapshot {
    // Formerly struct particle_system.
    std::vector<sim_particle> tumor_particles;
    std::vector<sim_particle> treg_particles;
    vec3 voxel_stochastic_vector = vec3(0,0,0);
    uint8_t voxel_stochastic_counter = 0;
    
    struct xorshift128p_state
    {
        uint64_t a = 0x19980724cafebabeu;
        uint64_t b = 0x20010106feef1f0fu;
    } rng;
};



// *** h_particles_snapshot map functions ***

// Map used for converting particles_snapshot handles to
// particles_snapshot instances. Use unique_ptr so that the actual
// snapshot objects won't move in memory unexpectedly during a rehash.
// (Actually, since unordered_map is node based this may not be needed,
// but I don't know the C++ standard well enough to count on it).
using hps_map_type = std::unordered_map<h_particles_snapshot,
                                        std::unique_ptr<particles_snapshot>>;
static hps_map_type& get_hps_map()
{
    static hps_map_type* hps_map = new hps_map_type;
    return *hps_map;
    // Use singleton pattern to avoid initialization order problems and
    // de-initailization crashes (the latter actually happened as the
    // hps_map destructor ran while some other thread was still operating
    // on it. Solved by just not running the destructor now.)
}

// Mutex used when accessing said map.
static std::mutex* hps_map_mutex = new std::mutex();

// Allocate a new particles_snapshot and return a handle to it. Write
// the new pointer through the optional argument.
inline h_particles_snapshot hps_alloc(particles_snapshot** out=nullptr)
{
    std::lock_guard<std::mutex> lock(*hps_map_mutex);
    auto& hps_map = get_hps_map();
    
    // We promised strictly positive hps in bedrock_particles.hpp
    static h_particles_snapshot hps_next = 1;
    static_assert(sizeof hps_next >= 8, "Need 64-bit h_particles_snapshot to avoid overflow.");
    
    auto* ps_raw_ptr = new particles_snapshot();
    std::unique_ptr<particles_snapshot> ps_ptr(ps_raw_ptr);
    
    auto hps = hps_next++;
    hps_map.emplace(hps, move(ps_ptr));
    if (out) *out = ps_raw_ptr;
    return hps;
}

// Get a pointer to a particles_snapshot given its handle. Returns a
// non-null pointer or throws. The pointer for a given hps is
// guaranteed not to change as long as hps remains valid (this is why
// std::unique_ptr is used).
inline particles_snapshot* ptr_from_hps(h_particles_snapshot hps)
{
    std::lock_guard<std::mutex> lock(*hps_map_mutex);
    auto& hps_map = get_hps_map();
    
    auto find_it = hps_map.find(hps);
    
    // Fail if handle was not valid.
    CHECK(find_it != hps_map.end());
    
    // Convert unique_ptr to raw pointer. If hps_map is working
    // correctly, ptr should not be null given that the
    // h_particle_snapshot was found in the map.
    auto* result = (*find_it).second.get();
    CHECK(result);
    
    return result;
}

// Check that the given particles_snapshot handle is valid. Throws
// particles_borked if it's not.
inline void check_hps(h_particles_snapshot hps)
{
    std::lock_guard<std::mutex> lock(*hps_map_mutex);
    auto& hps_map = get_hps_map();
    auto find_it = hps_map.find(hps);
    CHECK(find_it != hps_map.end());
}

// Deallocate the particles_snapshot with the given handle, and mark
// the handle as invalid.
inline void free_hps(h_particles_snapshot hps)
{
    std::lock_guard<std::mutex> lock(*hps_map_mutex);
    auto& hps_map = get_hps_map();
    
    // When a particles_snapshot is scheduled for deletion, it's not
    // actually deleted right away. It will be kept around for two
    // more free_hps calls using these pointers. This provides a bit
    // of robustness in case someone deletes a h_particles_snapshot
    // that is still being used by someone else. However, the hps
    // itself is made invalid immediately upon calling free_hps, so in
    // case this race condition does occur, we can check for it using
    // check_hps rather than crashing.
    static auto& ptr0 = *new std::unique_ptr<particles_snapshot>;
    static auto& ptr1 = *new std::unique_ptr<particles_snapshot>;
    
    auto find_it = hps_map.find(hps);
    
    // Error to free an already-invalid h_particles_snapshot.
    CHECK(find_it != hps_map.end());
    
    ptr0 = std::move(ptr1);
    ptr1 = std::move((*find_it).second);
    hps_map.erase(find_it);
}

inline size_t get_hps_count()
{
    std::lock_guard<std::mutex> lock(*hps_map_mutex);
    auto& hps_map = get_hps_map();
    return hps_map.size();
}

// NOT threadsafe.
// 
// Force deleting all particles_snapshots and their handles.
// Obviously VERY dangerous.  This is needed because we can't do
// proper cleanup when we exit play mode in the editor, which can
// quickly leak away GBs of memory.
inline void hack_delete_all_particles_snapshots()
{
    // Replace the mutex in case a thread was killed while it owned
    // the mutex.  (Yes, this leaks a mutex; C++ standard specifically
    // mentions it is undefined behaviour if a thread dies while
    // owning a mutex and this is the best I can do).
    hps_map_mutex = new std::mutex();
    auto& hps_map = get_hps_map();
    hps_map.clear();
}



// *** Random number generation ***

using rng_state = particles_snapshot::xorshift128p_state;

/* https://en.wikipedia.org/wiki/Xorshift */
/* The state must be seeded so that it is not all zero */
static uint64_t xorshift128p(rng_state* state)
{
    uint64_t t = state->a;
    uint64_t const s = state->b;
    state->a = s;
    t ^= t << 23; // a
    t ^= t >> 17; // b
    t ^= s ^ (s >> 26); // c
    state->b = t;
    return t + s;
}

// Return a random double in [low, high].
inline double
xorshift_randrange(rng_state* state, double low, double high)
{
    CHECK(is_real(low) && is_real(high));
    constexpr double normalizer = (1.0 / std::numeric_limits<uint64_t>::max());
    uint64_t n = xorshift128p(state);
    auto result = low + (high - low) * (n * normalizer);
    result = clamp(result, low, high);
    CHECK(is_real(result));
    return result;
}

// Generate a random 3-vector on the surface of a sphere of given radius.
static vec3 random_spherical_vector(
    rng_state* state,
    double radius)
{
    auto old_state = *state;
    
    // https://math.stackexchange.com/questions/1585975/how-to-generate-random-points-on-a-sphere
    // Use the fact that if you cut a sphere of a given radius with
    // two parallel planes, the area of the strip of spherical surface
    // between the planes depends only on the distance between the
    // planes, not on where they cut the sphere. -- Brian M. Scott
    auto theta = xorshift_randrange(state, 0, 2 * M_PI);
    auto z = xorshift_randrange(state, -radius, +radius);
    
    // Radius of the circular "slice" made by intersecting the sphere
    // with the plane z = z.
    auto slice_radius = sqrt(radius * radius - z * z);
    
    auto x = cos(theta) * slice_radius;
    auto y = sin(theta) * slice_radius;
    
    auto result = vec3(float(x), float(y), float(z));
    if (!is_real(result)) {
        fprintf(stderr, "old_state = %llu %llu\n", (unsigned long long) old_state.a, (unsigned long long) old_state.b);
        fprintf(stderr, "result = (%f %f %f)\n", result.x, result.y, result.z);
        fprintf(stderr, "theta = %f\n", theta);
        fprintf(stderr, "slice_radius = %f\n", slice_radius);
        fprintf(stderr, "x = %f\n", x);
        fprintf(stderr, "y = %f\n", y);
        fprintf(stderr, "z = %f (%.20e)\n", z, z);
        fprintf(stderr, "radius = %f (%.20e)\n", radius, radius);
        CHECK(0);
    }
    return result;
}

// Return a random prime number coprime to x, or return 1.
inline size_t random_coprime_prime_or_1(
    rng_state* state,
    size_t x)
{
    constexpr size_t primes[] = {
        2, 3, 5, 7,
        11, 13, 17, 19,
        23, 29, 31, 37,
        41, 43, 47, 53
    };
    constexpr auto prime_count = sizeof primes / sizeof primes[0];
    // Higher order xorshift "random" bits are better.
    size_t i = (xorshift128p(state) >> 50) % prime_count;
    size_t prime = primes[i];
    
    if (x % prime == 0) {
        return x % 59 == 0 ? 1 : 59;
    }
    return prime;
}



// *** Particle simulation utility functions ***

// VOXEL FILL MODEL
// 
// I will model the world as a grid of 1x1x1 unit voxels. Each step,
// we calcalate the total number of cells in the voxel by summing up
// the cell counts of the particles within the voxel (scaled by the
// amount of overlap for particles near the boundary). Particles try
// to move in/out of underfull/overfull voxels.
// 
// They are aligned on boundaries that are integers PLUS a
// voxel_stochastic_vector. This stochastic vector changes every step
// to avoid the aliasing problems that may occur if the voxels were
// always on integer boundaries.
// 
// We need to store these voxels in a dictionary and identify them
// uniquely. They are identified by integer 3-vectors (xyz), which are
// packed into a single integer (key_from_ivec3). The position of the
// voxel is between the corners (xyz) and (xyz + (1,1,1)), translated
// by the stochastic vector.
struct voxel
{
    // Used to sum up the cell count of this voxel.
    float cell_count = 0.0f;
};

using voxel_key_t = uint64_t;

static voxel_key_t voxel_key_from_ivec3(ivec3 arg)
{
    uint16_t x16 = arg.x;
    uint16_t y16 = arg.y;
    uint16_t z16 = arg.z;
    
    // Overflow check.
    CHECK(int16_t(x16) == arg.x);
    CHECK(int16_t(y16) == arg.y);
    CHECK(int16_t(z16) == arg.z);
    
    return voxel_key_t(x16) |
           voxel_key_t(y16) << 16 |
           voxel_key_t(z16) << 32;
}


// Dictionary type used for mapping voxel keys (see key_from_ivec3) to
// voxel structures (contains sum of voxel cell count).
using voxel_map = std::unordered_map<voxel_key_t, voxel>;

// Return a reference to the voxel in the voxel dictionary that
// corresponds to the voxel with lower corner voxel_id_ivec +
// voxel_stochastic_vector. See VOXEL FILL MODEL.
static voxel& get_voxel(
    voxel_map* map_ptr,
    ivec3 voxel_id_ivec)
{
    auto key = voxel_key_from_ivec3(voxel_id_ivec);
    auto& map = *map_ptr;
    
    auto end_it = end(map);
    
    auto it = map.find(key);
    
    // Key found -- return existing voxel.
    if (it != end_it) {
        return (*it).second;
    }
    // Key not found -- emplace a new voxel.
    else {
        voxel new_voxel;
        it = map.emplace_hint(it, key, move(new_voxel));
        return (*it).second;
    }
}

// Given the voxel dictionary and an integer 3-vector that identifies
// a voxel (see VOXEL FILL MODEL), return the number of cells in the
// named voxel.
inline float get_voxel_cell_count(
    const voxel_map& vmap,
    ivec3 voxel_id_ivec)
{
    auto key = voxel_key_from_ivec3(voxel_id_ivec);
    auto end_it = end(vmap);
    
    auto it = vmap.find(key);
    
    if (it != end_it) {
        return (*it).second.cell_count;
    }
    else {
        return 0;
    }
}

// We have to calculate the number of cells within each voxel. To do
// this, for each particle we calcalate the number of its cells that
// are in each voxel, and store this in a voxel_fill_list. The list may
// have multiple elements for particles near the boundary of a voxel --
// the particle's cells may be positioned in multiple voxels.
// 
// Moved into sim_particle struct due to incomplete type error.
using voxel_fill = sim_particle::voxel_fill;
using voxel_fill_list = sim_particle::voxel_fill_list;

inline void print_voxel_fill_list(voxel_fill_list const& vflst)
{
    printf("voxel_fill_list:\n");
    for (auto& vf : vflst) {
        auto xi = int(vf.voxel_id_ivec.x);
        auto yi = int(vf.voxel_id_ivec.y);
        auto zi = int(vf.voxel_id_ivec.z);
        auto cells = float(vf.cell_count_added);
        
        printf("(%i, %i, %i) %f\n", xi, yi, zi, cells);
    }
    printf("\n");
}

// Given a voxel_fill_list and an integer 3-vector that identifies
// a voxel (see VOXEL FILL MODEL), return the number of cells that
// would be added to the named voxel by the given voxel_fill_list.
inline float get_voxel_cell_count_added(
    const voxel_fill_list& fill_list,
    ivec3 voxel_id_ivec)
{
    auto voxel_id_ivec_equal = [voxel_id_ivec] (const voxel_fill& vf) -> bool
    {
        return vf.voxel_id_ivec == voxel_id_ivec;
    };
    
    auto begin_it = begin(fill_list);
    auto end_it = end(fill_list);
    
    auto find_it = find_if(begin_it, end_it, voxel_id_ivec_equal);
    bool found = (find_it != end_it);
    
    // Paranoid that I didn't use the standard container method correctly.
    CHECK(!found || (*find_it).voxel_id_ivec == voxel_id_ivec);
    
    return found ? (*find_it).cell_count_added : 0;
}

// Add the cells of a particle to the running sums of cells in the
// voxel dictionary, given the voxel_fill_list calculated for the
// particle. Also need to pass in the voxel_stochastic_vector as usual.
template <bool subtract=false>
inline void add_voxel_fill(
    voxel_map* map_ptr,
    const voxel_fill_list& fill_list)
{
    for (const voxel_fill& vf : fill_list) {
        voxel& v = get_voxel(map_ptr, vf.voxel_id_ivec);
        
        if (subtract) v.cell_count -= vf.cell_count_added;
        else v.cell_count += vf.cell_count_added;
    }
}

// Reverse operation of add_voxel_fill.
inline void subtract_voxel_fill(
    voxel_map* map_ptr,
    const voxel_fill_list& fill_list)
{
    add_voxel_fill<true> (map_ptr, fill_list);
}

// Calculate the intersection between a cubic particle of cells and
// the grid of voxels. Return a list of all voxels that contain cells
// represented by the particle, and the number of cells each.
// 
// The cube has the given position and "radius" (half edge-length) and
// has the same cell density as particles (particle_cells_per_cubic_unit).
static voxel_fill_list cube_make_voxel_fill_list(
    vec3 position,
    float cube_radius,
    vec3 voxel_stochastic_vector)
{
    voxel_fill_list result;

    // In case of NaN/inf.
    if (!is_real(position)) {
        fputs("cube_make_voxel_fill_list: non-real position?\n", stderr);
        return result;
    }
    if (!is_real(cube_radius)) {
        fputs("cube_make_voxel_fill_list: non-real radius?\n", stderr);
        return result;
    }
    
    // Cryptic variable names glossary:
    // r: "radius" (particle cube half-edge-width)
    // x, y, z: axis
    // L: low
    // H: high
    // i: integer
    // v: vector
    // t: translated (by -voxel_stochastic_vector)
    auto r = cube_radius;
    
    // Calculate the position of the lower and upper corners (low xyz
    // / high xyz) of the particle, which is modelled as an
    // axis-aligned cube of uniform cell density.
    vec3 rrr(r,r,r);
    vec3 vL = position - rrr;
    vec3 vH = position + rrr;
    // auto xL = vL.x;
    // auto yL = vL.y;
    // auto zL = vL.z;
    // auto xH = vH.x;
    // auto yH = vH.y;
    // auto zH = vH.z;
    
    // The particle intersects voxels with integer 3-vector
    // identifiers ranging from (xLi, yLi, zLi) to (xHi, yHi, zHi).
    // Remember that the grid of 1x1x1 voxels has boundaries that are
    // offset from the integer grid by voxel_stochastic_vector.
    using I = int32_t;
    auto vLt = vL - voxel_stochastic_vector;
    auto vHt = vH - voxel_stochastic_vector;
    auto xLi = I(floor(vLt.x));
    auto yLi = I(floor(vLt.y));
    auto zLi = I(floor(vLt.z));
    auto xHi = I(floor(vHt.x));
    auto yHi = I(floor(vHt.y));
    auto zHi = I(floor(vHt.z));
    
    // Hard to explain: these are the amounts "missing" between the
    // edge of the particle (cube) and the integer+stochastic
    // boundaries between voxels. See AkeleyLand/diagram-dont-delete.jpg
    auto xL_underfill = vLt.x - xLi;
    auto yL_underfill = vLt.y - yLi;
    auto zL_underfill = vLt.z - zLi;
    auto xH_underfill = 1 - (vHt.x - xHi);
    auto yH_underfill = 1 - (vHt.y - yHi);
    auto zH_underfill = 1 - (vHt.z - zHi);
    
    // Using the corners (xLi, yLi, zLi) and (xHi, yHi, zHi) as
    // boundary conditions, name (using an ivec3) every voxel that
    // intersects this particle and calculate the cell count this
    // particle contributes to each voxel. x/y/z_fill together
    // determine the proportion of a voxel that is filled with this particle
    // (1 for interior voxels; less for those at the boundary).
    // 
    // Again, look at AkeleyLand/diagram-dont-delete.jpg
    auto x_fill = 1.0f - xL_underfill;
    for (auto xi = xLi; xi <= xHi; ++xi) {
        x_fill -= (xi == xHi ? xH_underfill : 0);
        
        auto y_fill = 1.0f - yL_underfill;
        for (auto yi = yLi; yi <= yHi; ++yi) {
            y_fill -= (yi == yHi ? yH_underfill : 0);
            
            auto z_fill = 1.0f - zL_underfill;
            for (auto zi = zLi; zi <= zHi; ++zi) {
                z_fill -= (zi == zHi ? zH_underfill : 0);
                
                auto cell_count = x_fill * y_fill * z_fill
                                * particle_cells_per_cubic_unit;
                auto voxel_id_ivec = ivec3(xi, yi, zi);
                result.push_back(voxel_fill{voxel_id_ivec, cell_count});
                
                z_fill = 1.0f;
            }
            y_fill = 1.0f;
        }
        x_fill = 1.0f;
    }
    
    return result;
}

// Before calling, ensure that the radii (particle.update_radius) of
// the particle are correct.
// 
// Calculate the intersection between the passed particle and the grid
// of voxels. Return a list of all voxels that contain cells
// represented by the particle, and the number of cells each.
static voxel_fill_list make_voxel_fill_list(
    const sim_particle& p,
    vec3 voxel_stochastic_vector)
{
    return cube_make_voxel_fill_list(
        p.position,
        p.radius,
        voxel_stochastic_vector);
}


// *** Particle simulation step function implementation ***

// Perturb the voxel_stochastic_vector used for this particle system.
// 
// This is difficult to explain. Andrew figured out that actually
// randomizing the stochastic voxel vector every time was causing
// a lot of jitter, so instead we march it forward by a fixed
// small amount every simulation step.
static void update_voxel_stochastic_vector(particles_snapshot* ps)
{
    // The stochastic vector goes back-and-forth from 0 to 128 times
    // the stochastic step vector.
    auto counter = ++ps->voxel_stochastic_counter;
    float scalar = fabsf(float(int8_t(counter)));
    ps->voxel_stochastic_vector = stochastic_vector_step * scalar;
    // auto v = ps->voxel_stochastic_vector;
    // fprintf(stderr, "%f %f %f %f\n", scalar, v.x, v.y, v.z);
}

static float radius_cubed_from_cell_count(float cell_count)
{
    // Calculate radius of the cloud assuming (soft) maximum voxel
    // density.  Not 100% accurate because it doesn't account for the
    // fact that the particles themselves have radius.
    // 
    // V = 4/3 pi r^3
    // r^3 = 3/4 * V/pi
    // r = cube_root(3/4 * V/pi)
    // 
    // density = cell_count / V
    // V = cell_count / density
    // 
    // r = cube_root(3/4 * (cell_count/density) / pi)
    // 
    // BUT, it turns out that the cube of the radius is what I want
    // (explained later).
    auto density = soft_max_cells_per_voxel;
    auto V = cell_count / density;
    // Magic fudge factor.
    const auto magic = particle_init_volume_scale_magic;
    auto cloud_radius_cubed = float(0.75 * V / M_PI) * magic;
    return cloud_radius_cubed;
}

// Replaces the pointed-to vector of particles with an entirely new
// array of particles with (approximately) the expected cell
// count. The particle cloud produced has roughly the given radius
// (passed as the cube of radius). Mostly used to (re-)initialize the
// particle simulation when there are no particles.
// 
// Note: The motive for the explicit cloud_radius_cubed parameter is
// to support the Treg particle functionality bolted-on to the
// original particle simulation prototype. This was implemented as
// a complete particle array independent from the main tumor particles array.
// Originally, this function calculated the cloud radius itself from the
// given target cell count and the density constant. However, this was not
// suitable for the TReg particles cloud, since they usually were much fewer
// than the tumor particles, and would be all bunched-up in the center
// if their cloud radius was calculated from their particle count.
static void make_initial_particles(
    std::vector<sim_particle>* out,
    float target_cell_count,
    rng_state* rng,
    float cloud_radius_cubed)
{
    std::vector<sim_particle>& particles = *out;
    particles.clear();
    
    if ((long)target_cell_count <= 0) {
        particles.shrink_to_fit();
        return;
    }
    
    // I want to have an even distribution of particle sizes, from
    // minimum to maximum size. Calculate the number of particles
    // needed in order to maintain the expected average size.
    auto avg_cells_per_particle =
        double(min_cells_per_particle + max_cells_per_particle) * 0.5f;
    
    auto particle_count_f = target_cell_count / avg_cells_per_particle;
    auto particle_count = size_t(ceil(particle_count_f));
    
    // Special case: only one particle at the origin that contains all
    // the particles needed.
    if (particle_count <= 1) {
        particles.emplace_back();
        auto& particle = particles.front();
        
        // particle position *should* be zero due to default member
        // initializers.
        CHECK(particle.position == vec3(0,0,0));
        particle.cell_count = target_cell_count;
        particle.old_cell_count = target_cell_count;
        particle.update_radius();
        
        return;
    }
    
    // Avoid unnecessary vector reallocations.
    particles.reserve(particle_count);
    
    // Add particles to the simulation from smallest (min size) to
    // largest (max size), with sizes uniformly spaced.
    float cell_count_step = (max_cells_per_particle - min_cells_per_particle)
                          / (particle_count - 1);

    for (size_t p_idx = 0; p_idx < particle_count; ++p_idx) {
        auto cell_count = min_cells_per_particle + cell_count_step * p_idx;
        
        particles.emplace_back();
        auto& new_particle = particles.back();
        
        // Think of the spherical cloud of particles as a bunch of
        // spherical shells. We want the probability that a certain
        // shell is chosen for the next particle to be proportional to
        // the square of the shell's radius (because its surface area
        // is ~to square of its radius, and we want the density of the
        // cloud to be even).
        // 
        // I am not sure but I think sampling the cube of the radius from
        // a linear distribution, then cube rooting to get the actual radius
        // of the shell chosen for the next particle, is the correct way
        // to achieve this effect.
        auto r3 = xorshift_randrange(rng, 0, cloud_radius_cubed);
        auto radius = cube_root(r3);
        
        auto position = random_spherical_vector(rng, radius);
        
        new_particle.position = position;
        new_particle.old_position = position;
        new_particle.even_older_position = position;
        
        new_particle.cell_count = cell_count;
        new_particle.old_cell_count = cell_count;
        
        new_particle.update_radius();
    }
}

// Remove particles that are too small. The lost cells will be
// compensated for later. For the sake of visuals, this is done over
// the course of 2 simulation steps:
// 
// 1. Cell count set to zero, so that the particle visually implodes.
// 
// 2. Then, said particles, which should have 0 cell count, are
// actually removed from memory.
// 
// This also cycles the cell count history (used for visual
// interpolation; see sim_particle::old_radius).
template<bool treg_hack=false>
inline void remove_undersize_particles_and_push_cell_counts(
    std::vector<sim_particle>* particles_ptr,
    rng_state* rng=nullptr)
{
    size_t dst_idx = 0;
    auto& particles = *particles_ptr;
    size_t size = particles.size();
    size_t undersized_count = 0;
    bool treg_hack_triggered = false;
    
    for (size_t src_idx = 0; src_idx < size; ++src_idx) {
        auto& particle = particles.at(src_idx);
        
        // Remove from memory those particles with zero cell count by
        // not assigning them back to the particle array. Also take
        // this chance to remove particles that for some reason have
        // non-real position.
        if ((long)particle.cell_count == 0) {
            continue;
        }
        if (!is_real(particle.position)) {
            fputs("non-real particle position?\n", stderr);
            continue;
        }
        if (!is_real(particle.radius)) {
            fputs("non-real particle radius?\n", stderr);
            continue;
        }
        
        // See TREG HACK.
        if (treg_hack) {
            CHECK(rng != nullptr);
            auto random = xorshift_randrange(rng, 0, 1);
            treg_hack_triggered = (random < treg_hack_remove_treg_probability);
        }
        
        // Otherwise, remove 1-in-undersize_removal_ratio undersize
        // particles.  Note: for our purposes, removal is setting cell
        // count to 0, not actual removal from memory, since that's
        // what indicates removal to the user.
        bool undersized = (particle.cell_count < min_cells_per_particle);
        bool removed = false;
        
        if (undersized) {
            if (undersized_count % size_t(undersize_removal_ratio) == 0) {
                removed = true;
            }
            ++undersized_count;
        }
         
        if (removed || treg_hack_triggered) {
            particle.push_cell_count(0);
        } else {
            particle.push_cell_count(particle.cell_count);
        }
        
        particles[dst_idx++] = move(particle);
    }
    
    particles.resize(dst_idx);
}

// Scale the given particles array to match the target number of tumor
// cells by scaling the number of cells per particle (with some
// randomness).  Does NOT use push_cell_count; we are expecting the
// above function to have been run and done that for us already.
// 
// If, however, their are no particles (or they all have zero cell
// count), then instead this function (re-)creates the passed array of
// particles with the given target cell count (recreated as a particle
// cloud with given cubed radius), since scaling clearly will not have
// any effect in this case.
// 
// Returns the total number of cells AFTER the scaling is done.
static float scale_cell_counts_or_initialize(
    std::vector<sim_particle>* particles_ptr,
    float target_cell_count,
    rng_state* rng,
    float radius_cubed_for_initialization)
{
    auto& particles = *particles_ptr;
    float actual_total_cells = 0.0f;
    
    // Compute the total number of tumor cells represented
    // cosmetically by the particle simulation. This lets us see how
    // far off we are from the target cell count specified by targets.
    float total_cell_count = 0.0f;
    for (auto& particle : particles) {
        total_cell_count += particle.cell_count;
    }
    
    if ((long)total_cell_count == 0) {
        make_initial_particles(
            &particles,
            target_cell_count,
            rng,
            radius_cubed_for_initialization);
        
        for (auto& particle : particles) {
            total_cell_count += particle.cell_count;
        }
        return total_cell_count;
    }
    
    float volume_scale_factor = target_cell_count / total_cell_count;
    float volume_scale_factor_delta = volume_scale_factor - 1;
    
    
    for (auto& particle : particles) {
        auto random = xorshift_randrange(rng, 0.125, 1.875);
        if ((long)(volume_scale_factor_delta * 10) <= -3) {
        //if (volume_scale_factor_delta < -0.4f) {
            // Disable randomness in extreme cases to avoid negative
            // cell count.
            random = 1;
        }
        auto f = 1 + volume_scale_factor_delta * random;
        
        auto new_cell_count = particle.cell_count * f;
        particle.cell_count = new_cell_count;
        actual_total_cells += new_cell_count;
    }
    
    return actual_total_cells;
}

// Duplicate particles that are oversize and cut their cell count in
// half so that they appear to have split.
static void split_oversize_particles(
    std::vector<sim_particle>* particles_ptr,
    rng_state* rng)
{
    auto& particles = *particles_ptr;
    
    // Can't use iterators due to iterator invalidation.
    size_t original_size = particles.size();
    
    for (size_t idx = 0; idx < original_size; ++idx) {
        auto* ptr_particle = &particles.at(idx);
        
        if ((long)ptr_particle->cell_count > (long)max_cells_per_particle) {
            auto r = ptr_particle->radius;
            
            // Evenly distribute the cells between the "two" particles
            // (duplication will occur soon).
            ptr_particle->cell_count *= 0.5f;
            ptr_particle->update_radius();
            
            // Now duplicate.
            particles.push_back(*ptr_particle);
            auto& new_particle = particles.back();
            ptr_particle = &particles.at(idx); // In case of reallocation.
            
            // Give random offset to the two particles so they don't overlap.
            auto random_vector = random_spherical_vector(rng, r);
            ptr_particle->position += random_vector;
            new_particle.position -= random_vector;
            
            // Only the old particle keeps any attached th/ctl particle.
            new_particle.has_ctl = false;
            new_particle.has_t_helper = false;
        }
    }
}

// For every oversived particle in the passed treg particles array,
// split into two smaller particles, with the newer of the two placed
// at a position matching a randomly-selected tumor particle.
static void split_teleport_oversize_treg_particles(
    std::vector<sim_particle>* treg_particles_ptr,
    const std::vector<sim_particle>& tumor_particles,
    rng_state* rng)
{
    auto& treg_particles = *treg_particles_ptr;
    size_t size = treg_particles.size();
    for (size_t i = 0; i < size; ++i) {
        auto* particle_ptr = &treg_particles[i];
        if ((long)particle_ptr->cell_count > (long)max_cells_per_particle) {
            auto halved_cell_count = particle_ptr->cell_count * 0.5f;
            particle_ptr->cell_count = halved_cell_count;
            particle_ptr->update_radius();
            
            // Adding new particle invalidates particle_ptr.
            particle_ptr = nullptr;
            treg_particles.emplace_back();
            auto& new_particle = treg_particles.back();
            
            // Select random position for the new particle.
            size_t random_idx = size_t(xorshift_randrange(
                rng, 0, tumor_particles.size() + 1));
            const sim_particle* random_particle =
                random_idx < tumor_particles.size()
                ? &tumor_particles.at(random_idx) : nullptr;
            vec3 position = random_particle
                          ? random_particle->position : vec3(0,0,0);
            float teleport_offset_distance = random_particle
                                           ? random_particle->radius : 1.0f;
            position += random_spherical_vector(rng, teleport_offset_distance);
            
            // Set the new particle's positions + history.
            new_particle.position = position;
            new_particle.old_position = position;
            new_particle.even_older_position = position;
            
            // Make new TRegs fade into existence from 0 radius.
            new_particle.old_cell_count = 0;
            new_particle.cell_count = halved_cell_count;
            new_particle.update_radius();
        }
    }
}

// See VOXEL FILL MODEL. Look at the cell counts of the array of
// particles passed in, distribute them among the voxels, and add
// those counts to the counts maintained in the passed voxel_map.
// 
// This also updates particles' tmp_voxel_fill_list with information
// about their cell count intersections with the voxel grid.
static void add_cells_to_vmap_update_voxel_fill_lists(
    voxel_map* vmap_ptr,
    std::vector<sim_particle>* particles_ptr,
    vec3 voxel_stochastic_vector)
{
    voxel_map& vmap = *vmap_ptr;
    
    for (auto& particle : *particles_ptr) {
        CHECK(particle.cell_count >= 0);
        particle.update_radius();
        auto vflst = make_voxel_fill_list(particle, voxel_stochastic_vector);
        add_voxel_fill(&vmap, vflst);
        particle.tmp_voxel_fill_list = move(vflst);
    }
}

// Sanity check implementations: see that the total number of cells
// based on the voxel map matches the actual total number of cells in
// the particle simulation.
static float voxel_map_sanity_check_delta_proportion(
    const voxel_map& vmap,
    float expected_cell_count)
{
    auto delta = expected_cell_count;
    for (auto& kv_pair : vmap) {
        auto voxel_cell_count = kv_pair.second.cell_count;
        delta -= voxel_cell_count;
    }
    auto delta_proportion =
        expected_cell_count == 0 ? 0 : delta / expected_cell_count;
    
    return delta_proportion;
}

static void voxel_map_sanity_check(
    const voxel_map& vmap,
    float expected_cell_count)
{
    auto delta_proportion = voxel_map_sanity_check_delta_proportion(
        vmap,
        expected_cell_count);
    
    CHECK(-0.01f < delta_proportion && delta_proportion < 0.01f);
}

// Split into 2 functions in source code so I can see WHICH assertion
// was tripped.
static void voxel_map_sanity_check_generous(
    const voxel_map& vmap,
    float expected_cell_count)
{
    auto delta_proportion = voxel_map_sanity_check_delta_proportion(
        vmap,
        expected_cell_count);

    CHECK(-0.03f < delta_proportion && delta_proportion < 0.03f);
}

// Given the voxel map that's been filled with the cells of all
// particles in the simulation (see VOXEL FILL MODEL), return the cost
// of moving the given particle by the given displacement
// vector. Positive = Bad.
// 
// Assumes that the tmp_voxel_fill_list of the particle is accurate.
// This function writes to out_fill_list what the new voxel_fill_list
// for this particle should be if the considered displacement is
// actually done.
// 
// Be sure to check that cost_of_doing_nothing matches this function
// in case this cost function changes (I considered doing nothing to
// be a sufficiently special case for a new function).
static float cost_function(
    const sim_particle& p,
    vec3 displacement,
    const voxel_map& vmap,
    vec3 voxel_stochastic_vector,
    voxel_fill_list* out_fill_list)
{
    const auto& old_fill_list = p.tmp_voxel_fill_list;
    
    // Write as if the particle's future position is current position
    // + displacement even though we are only considering the move.
    vec3 new_position = p.position + displacement;
    
    *out_fill_list = 
        cube_make_voxel_fill_list(
            new_position,
            p.radius,
            voxel_stochastic_vector);
    
    // Calculate the number of excess cells in voxels that would be caused
    // by the particle's presence at new_position (worlds fail me).
    float overfill_cells = 0;
    float soft_overfill_cells = 0;
    for (const voxel_fill& vf : *out_fill_list) {
        // New cell count of this voxel is the amount for it recorded
        // in the voxel map, minus the original contribution of cells
        // by this particle (before the move out we're considering)
        // plus the new contribution of this particle at its new
        // position.
        ivec3 voxel_id_ivec = vf.voxel_id_ivec;
        float voxel_cell_count = get_voxel_cell_count(vmap, voxel_id_ivec);
        float cells_added = vf.cell_count_added;
        float cells_lost =
            get_voxel_cell_count_added(old_fill_list, voxel_id_ivec);
        voxel_cell_count += (cells_added - cells_lost);
        
        // Minimum is in case there's no overfill; maximum so that the
        // particle is not blamed for more than its contribution to
        // overfill.
        float voxel_overfill_cells = clamp(
            voxel_cell_count - max_cells_per_voxel, 0.0f, cells_added);

        float soft_voxel_overfill_cells = clamp(
            voxel_cell_count - soft_max_cells_per_voxel, 0.0f, cells_added);
        
        overfill_cells += voxel_overfill_cells;
        soft_overfill_cells += soft_voxel_overfill_cells;
    }
    
    auto magnitude = [] (vec3 v) -> float
    {
        return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
    };
    auto old_distance_from_origin = magnitude(p.position);
    auto new_distance_from_origin = magnitude(new_position);
    auto delta = new_distance_from_origin - old_distance_from_origin;
    
    // If moving the particle to the new position would still cause
    // overfill, bias in favor of displacements that point away from
    // the origin.
    if ((long)overfill_cells > 0) {
        auto desperation_magic = (0.5f + p.desperation)
                               * desperation_distance_from_origin_magic;
        auto distance_cost = -delta
                           * cost_distance_from_origin_magic
                           * desperation_magic;
        auto cost = std::max(0.0f, overfill_cells + distance_cost);
        return cost;
    }
    
    // Otherwise, bias in favor of displacements that point to the
    // origin.  However, here we use a somewhat lower limit of voxel
    // cell count, in order to avoid pathological jittering (particle
    // moves inwards to a new voxel, sees it's overfull and regrets
    // the decision, moves out again, then in again, etc.)
    if ((long)soft_overfill_cells <= 0) {
        return +delta * cost_distance_from_origin_magic;
    }
    
    return 0.0f;
}

// Baseline cost of not moving the referenced particle at all. Requires
// that the voxel_map has been filled with the cell counts of all the
// particles in the simulation and that p.tmp_voxel_fill_list is accurate.
static float cost_of_doing_nothing(
    const sim_particle& p,
    const voxel_map& vmap)
{
    const auto& fill_list = p.tmp_voxel_fill_list;
    
    // Calculate the number of excess cells in voxels caused by the
    // particle being in its current position.
    float overfill_cells = 0;
    for (const voxel_fill& vf : fill_list) {
        float cells_added = vf.cell_count_added;
        float voxel_cell_count = get_voxel_cell_count(vmap, vf.voxel_id_ivec);
        float voxel_overfill_cells = clamp(
            voxel_cell_count - max_cells_per_voxel, 0.0f, cells_added);
        
        overfill_cells += voxel_overfill_cells;
    }
    
    return overfill_cells;
}

// For every particle in the passed array, try to move it in a
// direction that minimizes its cost function. This requires read and
// write access to the voxel map: read, to see when particles are in
// areas that are too dense; write, to modify the map when we move
// particles (and therefore the cell densities they represent).
static void push_new_particle_positions_to_minimize_cost(
    voxel_map* vmap_ptr,
    std::vector<sim_particle>* particles_ptr,
    vec3 voxel_stochastic_vector,
    rng_state* rng)
{
    voxel_map& vmap = *vmap_ptr;
    for (auto& particle : *particles_ptr) {
        // Generate random displacement vectors and choose to
        // displace the particle by the one with lowest cost.
        // We also first consider doing no displacement at all.
        vec3 displacement(0,0,0);
        float cost = cost_of_doing_nothing(particle, vmap);
        voxel_fill_list try_fill_list;
        bool particle_moved = false;
        
        CHECK(stochastic_tries >= 0);
        
        for (size_t i = 0; i < size_t(stochastic_tries); ++i) {
            auto step_size =
                (particle.desperation + 0.5f) * stochastic_step_size_step;
            vec3 try_displacement =
                random_spherical_vector(rng, step_size);
            
            float try_cost = cost_function(
                particle,
                try_displacement,
                vmap,
                voxel_stochastic_vector,
                &try_fill_list);
            
            // Use < instead of <= because, all else being equal, we
            // would prefer to not move the particle to minimize jitter.
            if ((long)try_cost < (long)cost) {
                cost = try_cost;
                displacement = try_displacement;
                particle_moved = true;
            }
        }
        
        
        
        // If doing something other than staying still has a lower
        // cost, then we choose to move the particle. This requires
        // updating the voxel map with the new positions of the cells
        // contained in the moved particle.
        if (particle_moved) {
            const voxel_fill_list& old_fill_list = particle.tmp_voxel_fill_list;
            subtract_voxel_fill(&vmap, old_fill_list);
            add_voxel_fill(&vmap, try_fill_list);
            particle.tmp_voxel_fill_list = move(try_fill_list);
            ++particle.desperation;
        }
        else {
            // Failed...
            if ((long)(cost*100) > 0) {
                ++particle.desperation;
            }
            else {
                particle.desperation = 0;
            }
        }
        
        // Need to use push_new_position to preserve position history
        // for visual interpolation.
        particle.push_new_position(particle.position + displacement);
    }
}

// Tumor particles can have a smaller red particle attached to
// represent CTL. Remove some tumor particles with attached ctl to
// give visual indication of CTL = tumor cell death. Then randomly add
// CTL particles until we reach the target CTL cell count (if we can).
// This also adds T helper cells, although I am not too motivated to
// reach the target T Helper count.
static void update_ctl_t_helper_particles(
    std::vector<sim_particle>* tumor_particles_ptr,
    float target_ctl_cells,
    float target_t_helper_cells,
    rng_state* rng)
{
    auto& particles = *tumor_particles_ptr;
    float total_ctl_cells = 0.0f;
    float total_t_helper_cells = 0.0f;
    
    // Randomly remove some particles with CTL attached by setting
    // their cell counts to zero. They will be actually removed from
    // memory later in the undersize particle check step.
    // Q: How expensive are the branch mispredicts caused here?
    for (sim_particle& p : particles) {
        if (p.has_ctl) {
            auto random = xorshift_randrange(rng, 0, 1);
            if (random < ctl_death_probability) {
                p.cell_count = 0;
                p.radius = 0;
            }
            total_ctl_cells += p.cell_count * ctl_particle_cell_ratio;
        }
        if (p.has_t_helper) {
            total_t_helper_cells += p.cell_count * t_helper_particle_cell_ratio;
        }
    }
    
    size_t size = particles.size();
    size_t prime = random_coprime_prime_or_1(rng, size);
    size_t offset = size_t(xorshift_randrange(rng, 0, size-1));
    
    // Now add ctls randomly until we meet the target ctl cell count.
    for (size_t i = 0; i < size; ++i) {
        // Iterate through the particle list using a random prime
        // step.  This may be expensive (cache performance?) but
        // reduces memory-order artifacts.
        sim_particle& p = particles[(i * prime + offset) % size];
        
        bool needs_ctl = total_ctl_cells < target_ctl_cells;
        bool needs_t_helper = total_t_helper_cells < target_t_helper_cells;
        
        if (!p.has_ctl && needs_ctl) {
            auto particle_ctl_cells = p.cell_count * ctl_particle_cell_ratio;
            total_ctl_cells += particle_ctl_cells;
            
            p.has_ctl = true;
            auto offset = random_spherical_vector(rng, 1);
            p.ctl_particle_normalized_offset = offset;
        }
        
        // Can only add t helpers to cells with ctl.
        if (p.has_ctl && !p.has_t_helper && needs_t_helper) {
            p.has_t_helper = true;
            auto particle_t_helper_cells = p.cell_count
                                         * t_helper_particle_cell_ratio;
            total_t_helper_cells += particle_t_helper_cells;
        }
        
        if (!needs_ctl && !needs_t_helper) {
            break;
        }
    }
}

// Perform one step of the particle simulation with the given target
// cell counts, modifying the particles_snapshot in-place. See
// bedrock_particles.hpp
static void simulation_step_fn_inplace(
    particles_snapshot* ps,
    float target_tumor_cells,
    float target_ctl_cells,
    float target_t_helper_cells,
    float target_treg_cells)
{
    auto& tumor_particles = ps->tumor_particles;
    auto& treg_particles = ps->treg_particles;
    update_voxel_stochastic_vector(ps);
    rng_state* rng = &ps->rng;
    
    // Defensive programming: flush negative target cell counts to 0.
    if ((long)target_tumor_cells < 0) target_tumor_cells = 0;
    if ((long)target_ctl_cells < 0) target_ctl_cells = 0;
    if ((long)target_t_helper_cells < 0) target_t_helper_cells = 0;
    if ((long)target_treg_cells < 0) target_treg_cells = 0;
    
    // Remove too-small particles.
    remove_undersize_particles_and_push_cell_counts(&tumor_particles);
    remove_undersize_particles_and_push_cell_counts<true>(
        &treg_particles, // See TREG HACK
        rng);
    
    // Scale cell counts of particles to match expected cell
    // count. Deal with consequences of this later. Re-initializes the
    // particle cloud if needed (if 0 particles).
    float radius_cubed_for_initialization = radius_cubed_from_cell_count(
        target_tumor_cells + target_treg_cells);
    
    float actual_total_cells = scale_cell_counts_or_initialize(
        &tumor_particles,
        target_tumor_cells,
        rng,
        radius_cubed_for_initialization);
    
    actual_total_cells += scale_cell_counts_or_initialize(
        &treg_particles,
        target_treg_cells,
        rng,
        radius_cubed_for_initialization);
    
    // Time to go through and sum up the number of cells in each voxel,
    // which allows us to find areas of too-high density, and inner
    // regions with too-low density (which we can move particles into).
    voxel_map vmap;
    auto vsv = ps->voxel_stochastic_vector;
    add_cells_to_vmap_update_voxel_fill_lists(&vmap, &tumor_particles, vsv);
    add_cells_to_vmap_update_voxel_fill_lists(&vmap, &treg_particles, vsv);
    voxel_map_sanity_check(vmap, actual_total_cells);
    
    // Now try to reposition particles in a way that makes them happy.
    // This also calls particle::push_new_position to cycle position
    // history for visual interpolation between simulation steps.
    // Randomly decide between updating tumor particles first, and
    // updating tregs first.
    bool tumors_first = (xorshift128p(rng) >> 60) & 1;
    auto& first_particles = tumors_first ? tumor_particles : treg_particles;
    auto& second_particles = tumors_first ? treg_particles : tumor_particles;
    auto* f = &push_new_particle_positions_to_minimize_cost;
    f(&vmap, &first_particles, vsv, rng);
    f(&vmap, &second_particles, vsv, rng);
    
    // Use more generous delta for this check due to all the
    // adding/subtracting done when moving the particles around.
    voxel_map_sanity_check_generous(vmap, actual_total_cells);
    
    // Add / remove attached CTL particles (attached to tumor cells only).
    update_ctl_t_helper_particles(
        &tumor_particles,
        target_ctl_cells,
        target_t_helper_cells,
        rng);
    
    // Split up too-big tumor particles.
    split_oversize_particles(&tumor_particles, rng);
    
    // However, for oversive TReg particles, just split them into 2
    // smaller tregs but spawn the new particle in some random
    // place. This was needed because of noticable clumping with a
    // splitting-based approach.
    split_teleport_oversize_treg_particles(
        &treg_particles,
        tumor_particles,
        rng);
    
    // Clear temporary voxel_fill_list for performance (avoids
    // unneccessary copies).
    for (auto& particle : tumor_particles) {
        particle.tmp_voxel_fill_list.clear();
    }
    for (auto& particle : treg_particles) {
        particle.tmp_voxel_fill_list.clear();
    }
}

return_code Bedrock_simulation_step_fn_inplace(
    h_particles_snapshot hps,
    float target_tumor_cells,
    float target_ctl_cells,
    float target_t_helper_cells,
    float target_treg_cells)
{
    try {
        // Program defensively: since we are mutating the
        // particles_snapshot, we should have sole access to it.
        // So, we can move it to a local copy (which we know won't
        // dissappear from under us), modify its contents in place,
        // then write back after we verify the original pointer is
        // still valid. So, if someone deallocates hps while we are
        // running, we get a CHECK fail, instead of a hard crash.
        particles_snapshot* original_ps_ptr = ptr_from_hps(hps);
        particles_snapshot local_copy_ps = move(*original_ps_ptr);
        simulation_step_fn_inplace(
            &local_copy_ps,
            target_tumor_cells,
            target_ctl_cells,
            target_t_helper_cells,
            target_treg_cells);
        check_hps(hps);
        *original_ps_ptr = move(local_copy_ps);
        // Potential race between these last 2 steps, only if the potential
        // bug we're checking for happens *exactly* now.
        
        // Also note that the write-back does not occur if an
        // exception occurs. This is GOOD -- if there's something
        // subtly wrong that may crash the system, then we should get
        // occasional user reports of disappearing particle systems
        // that would warn us about it.
    }
    catch (particles_borked& pb) {
        return pb.error_code < 0 ? pb.error_code : -1234567;
    }
    catch (...) {
        return -1234567;
    }
    return hps;
}

// Convert a particle system into a list of spheres to draw
// (positions, color, radius). bz is used to visually interpolate
// between the particle system's current particle positions and sizes,
// and the positions and sizes of the previous simulation step.
// 
// Assumes that update_radius was already called (so that visual
// radius corresponds to cell count)
static std::vector<visual_particle>
make_visual_particle_list(const particles_snapshot& ps, bezier_coefficients bz)
{
    std::vector<visual_particle> vp_array;
    for (auto& particle : ps.tumor_particles) {
        auto count = particle.cell_count;
        float a = (count - min_cells_per_particle)
                / (max_cells_per_particle - min_cells_per_particle);
        a = 1 - clamp(a, 0, 1);
        
        visual_particle tumor_vp;
        
        tumor_vp.red = a * tumor_red;
        tumor_vp.green = a * tumor_green;
        tumor_vp.blue = a * tumor_blue;
        
        particle.interpolate_display_position(bz, &tumor_vp);
        
        vp_array.push_back(tumor_vp);
        
        if (particle.has_ctl) {
            visual_particle ctl_vp;
            float ctl_distance = tumor_vp.radius * ctl_particle_distance_ratio;
            vec3 ctl_offset_vector = ctl_distance
                                   * particle.ctl_particle_normalized_offset;
            
            ctl_vp.x = tumor_vp.x + ctl_offset_vector.x;
            ctl_vp.y = tumor_vp.y + ctl_offset_vector.y;
            ctl_vp.z = tumor_vp.z + ctl_offset_vector.z;
            ctl_vp.radius = tumor_vp.radius * ctl_particle_radius_ratio;
            ctl_vp.red = ctl_red;
            ctl_vp.green = ctl_green;
            ctl_vp.blue = ctl_blue;
            
            vp_array.push_back(ctl_vp);
        }
        
        if (particle.has_t_helper) {
            visual_particle t_helper_vp;
            float t_helper_distance = tumor_vp.radius * t_helper_particle_distance_ratio;
            vec3 t_helper_offset_vector = t_helper_distance
                                        * particle.ctl_particle_normalized_offset;
            
            t_helper_vp.x = tumor_vp.x + t_helper_offset_vector.x;
            t_helper_vp.y = tumor_vp.y + t_helper_offset_vector.y;
            t_helper_vp.z = tumor_vp.z + t_helper_offset_vector.z;
            t_helper_vp.radius = tumor_vp.radius * t_helper_particle_radius_ratio;
            t_helper_vp.red = t_helper_red;
            t_helper_vp.green = t_helper_green;
            t_helper_vp.blue = t_helper_blue;
            
            vp_array.push_back(t_helper_vp);
        }
    }
    
    for (auto& particle : ps.treg_particles) {
        visual_particle treg_vp;
        treg_vp.red = treg_red;
        treg_vp.green = treg_green;
        treg_vp.blue = treg_blue;
        
        particle.interpolate_display_position(bz, &treg_vp);
        
        vp_array.push_back(treg_vp);
    }
    return vp_array;
}

// See bedrock_particles.hpp
// 
// If you modify this function, make sure
// Bedrock_interpolate_visual_particles_bytes_needed
// gets updated as well.
return_code Bedrock_interpolate_visual_particles(
    h_particles_snapshot hps,
    float t,
    visual_particle* vp_ptr,
    size_t max_bytes)
{
    try {
        // We promised in bedrock_particles.hpp that t would be
        // clamped to [0, 1] by us if needed.
        if ((long)t < 0) t = 0;
        if ((long)t >= 1) t = 1;
        bezier_coefficients bz(t);
        
        // In principle using a vector is pointlessly inefficient when we
        // could just write directly to the output buffer, but I don't want
        // to risk screwing this up while porting from C++ to C# just to
        // save a few heap allocations.
        // 
        // (Actually this is really sad: eventually this data is
        // needed in a C# List which is basically the same thing as a
        // std::vector [only somewhat worse] yet we'll have to copy
        // this data two, maybe three times! More if we didn't have
        // enough room and have to retry!)
        // 
        // If this is actually a noticable performance issue and you
        // have an idea to fix it, go ahead and try it out (but please
        // recompile for linux if you do so.)
        auto vp_vector = make_visual_particle_list(*ptr_from_hps(hps), bz);
        size_t vp_count = vp_vector.size();
        auto bytes = vp_count * sizeof(vp_vector[0]);
        CHECK(bytes <= max_bytes);
        memcpy(vp_ptr, vp_vector.data(), bytes);
        
        return_code to_return = vp_count;
        CHECK(size_t(to_return) == vp_count); // overflow check.
        check_hps(hps);
        return to_return;
    }
    catch (particles_borked& pb) {
        return pb.error_code < 0 ? pb.error_code : -1234567;
    }
    catch (...) {
        return -1234567;
    }
}

return_code Bedrock_interpolate_visual_particles_bytes_needed(
    h_particles_snapshot hps)
{
    try {
        const particles_snapshot& ps = *ptr_from_hps(hps);
        size_t vp_count = ps.treg_particles.size();
        
        for (const sim_particle& sp : ps.tumor_particles) {
            vp_count++;
            vp_count += (sp.has_ctl ? 1 : 0);
            vp_count += (sp.has_t_helper ? 1 : 0);
        }
        
        check_hps(hps);
        return vp_count * sizeof(visual_particle);
    }
    catch (particles_borked& pb) {
        return pb.error_code < 0 ? pb.error_code : -1234567;
    }
    catch (...) {
        return -1234567;
    }
}

h_particles_snapshot hpsBedrock_new_particles_snapshot()
{
    try {
        return hps_alloc();
    }
    catch (particles_borked& pb) {
        return pb.error_code < 0 ? pb.error_code : -1234567;
    }
    catch (...) {
        return -1234567;
    }
}

return_code Bedrock_free_hps(h_particles_snapshot hps)
{
    try {
        free_hps(hps);
        return 0;
    }
    catch (particles_borked& pb) {
        return pb.error_code < 0 ? pb.error_code : -1234567;
    }
    catch (...) {
        return -1234567;
    }
}

return_code hpsBedrock_copy_particles_snapshot(h_particles_snapshot hps_src)
{
    try {
        particles_snapshot* src_ptr = ptr_from_hps(hps_src);
        particles_snapshot* dest_ptr = nullptr;
        auto hps_dest = hps_alloc(&dest_ptr);
        *dest_ptr = *src_ptr;
        return hps_dest;
    }
    catch (particles_borked& pb) {
        return pb.error_code < 0 ? pb.error_code : -1234567;
    }
    catch (...) {
        return -1234567;
    }
}

return_code Bedrock_get_hps_count()
{
    try {
        return return_code(get_hps_count());
    }
    catch (particles_borked& pb) {
        return pb.error_code < 0 ? pb.error_code : -1234567;
    }
    catch (...) {
        return -1234567;
    }
}

return_code Bedrock_hack_delete_all_particles_snapshots()
{
    try {
        hack_delete_all_particles_snapshots();
        return 0;
    }
    catch (particles_borked& pb) {
        return pb.error_code < 0 ? pb.error_code : -1234567;
    }
    catch (...) {
        return -1234567;
    }
}
