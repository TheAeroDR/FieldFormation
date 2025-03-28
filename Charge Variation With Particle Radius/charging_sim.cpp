//includes
#include <iostream>
#include <numbers>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <ctime>
#include <chrono>
#include <iomanip>
#include <filesystem>
#include <csignal>
#include <functional>


// Define class containing all model parameters, and initialise
class Params {
public:
    double box_length = 2000.0; //um
    double packing_density = 0.10;
    int init_particles = 250;
    int n_particles = 0;
    int n_eqil_steps = 10000000;
    int n_max_arrangement = 100000;
    double speed_l = 1.0; //um/s
    double speed_u = 2.0; //um/s
    double density = 1.40e-15; //1.40g/cm^3 from mgs1 data sheet converted to kg/um^3
    double charge_density = 1; //electrons/um^2
    double temperature = 300.00; //K
};

// Define class containing all particle parameters and some basic functions
class Particle {
public:
    double r; // um
    double mass; //kg
    int initial_high_e;
    int high_e;
    int low_e;
    double v[3]; //um/s
    double pos[3]; //um

    // Function to set particle mass
    void spherical_mass(double density){
        // Spherical volume x density
        mass = (4.0/3.0) * M_PI * std::pow(r, 3.0) * density;
    }
    
    // Function to set high energy electrons to start value 
    void init_high_e(double q_density){
        // Spherical surface area x charge density
        high_e = static_cast<int>(round(4.0 * M_PI * std::pow(r, 2.0) * q_density));
        initial_high_e = high_e;
    }

    // Function to randomly position particles in the box - ensureing they don't start within a radius of any edge (inside the wall)
    void init_position(double box){
        // Create random item
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis( 1.0001 * r, box - 1.0001 * r);

        // Generate random numbers from the object
        pos[0] = dis(gen);
        pos[1] = dis(gen);
        pos[2] = dis(gen);
    }

    // Function to initialise velocity based off the a random distribution with between u_l and u_u in each direction
    void old_init_velocity(double u_l, double u_u){
        // Create random items
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> sign(0, 1);
        std::uniform_real_distribution<> dis(u_l, u_u);

        // Generate random numbers
        v[0] = ((sign(gen) == 0) ? -1 : 1) * dis(gen);
        v[1] = ((sign(gen) == 0) ? -1 : 1) * dis(gen);
        v[2] = ((sign(gen) == 0) ? -1 : 1) * dis(gen);
    }

    // Function to initialise velocity based off a random distribution between u_l and u_u in magnitude in a random direction
    void init_velocity(double u_l, double u_u){
        // Create random item
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(-1.0, 1.0);
        std::uniform_real_distribution<> v_dis(u_l, u_u);

        // Generate random numbers
        double vx = dis(gen);
        double vy = dis(gen);
        double vz = dis(gen);

        // Determine magnitude, mean and difference
        double mag = std::sqrt(vx * vx + vy * vy + vz * vz);

        // Assign numbers based off the mean + half the difference, times the component in each direction
        v[0] = v_dis(gen) * (vx / mag);
        v[1] = v_dis(gen) * (vy / mag);
        v[2] = v_dis(gen) * (vz / mag);
    }

    // Function to initialise velocity based off a Maxwell-Boltzmann kinetic energy distribution magnitude in a random direction
    void alt_init_velocity(double temperature) {
        // Create random device
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(-1.0, 1.0);

        // Generate random numbers
        double vx = dis(gen);
        double vy = dis(gen);
        double vz = dis(gen);

        // Calculate magnitude
        double mag = std::sqrt(vx * vx + vy * vy + vz * vz);

        // Work out boltzmann rms velocity for initialisation temperature
        double boltzmann_constant = 1.38064852e-11; //  um^2 kg s^-2 K^-1
        double velocity_rms = std::sqrt((3.0 * boltzmann_constant * temperature) / (mass * 1e-18));

        // Scale the rms onto the random directional unit vector
        v[0] = (vx / mag) * velocity_rms;
        v[1] = (vy / mag) * velocity_rms;
        v[2] = (vz / mag) * velocity_rms;
    }      

    // Function to update position 
    void time_step(double time){
        // New position = old + distance moved in time step
        pos[0] = pos[0] + v[0] * time;
        pos[1] = pos[1] + v[1] * time;
        pos[2] = pos[2] + v[2] * time; 
    }

    // Function to get kinetic energy of particle
    double getKineticEnergy() const {
        double velocityMagnitudeSquared = 0.0;
        for (int i = 0; i < 3; ++i) {
            velocityMagnitudeSquared += v[i] * v[i];
        }
        return 0.5 * mass * velocityMagnitudeSquared;
    }
};

// Global variables to store params and parts
Params params;
std::vector<Particle> parts;

// Define a custom tqdm class for c++
class Tqdm {
public:
    Tqdm(int total) : total(total), n(0), step(0) {
        start_time = std::chrono::system_clock::now();
    }
 
    ~Tqdm() {
        auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - start_time);
        std::cout << "Total time: " << elapsed_seconds.count() << "s" << std::endl;
    }
 
    void update(int delta, int timestep) {
        n = delta;
        step = timestep;
        display();
    }
 
private:
    int total;
    int n;
    int step;
    std::chrono::time_point<std::chrono::system_clock> start_time;
 
    void display() {
        double progress = static_cast<double>(n) / total * 100;
 
        auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - start_time);
        double elapsed_time = elapsed_seconds.count();
        double avg_time_per_step = elapsed_time / n;
        double remaining_time = avg_time_per_step * (total - n);
        double iterations_per_second = step / elapsed_time;
 
        int hours = static_cast<int>(elapsed_time) / 3600;
        int minutes = static_cast<int>(elapsed_time) % 3600 / 60;
        int seconds = static_cast<int>(elapsed_time) % 60;
 
        int rem_hours = static_cast<int>(remaining_time) / 3600;
        int rem_minutes = static_cast<int>(remaining_time) % 3600 / 60;
        int rem_seconds = static_cast<int>(remaining_time) % 60;
 
        std::stringstream ss;
        ss << "[";
        int bar_width = 50;
        int completed = static_cast<int>(progress / 100 * bar_width);
        for (int i = 0; i < bar_width; ++i) {
            if (i < completed) ss << "=";
            else if (i == completed) ss << ">";
            else ss << " ";
        }
 
        ss << "] " << std::setw(6) << std::setfill(' ') << std::round(progress * 100) / 100 << "% ";
        int num_width = std::to_string(total).length();
        ss << std::setw(num_width) << std::setfill(' ') << n << "/" << total << " ";
        ss << "[" << std::setfill('0') << std::setw(2) << hours << ":" << std::setw(2) << minutes << ":" << std::setw(2) << seconds;
        ss << "<" << std::setw(2) << rem_hours << ":" << std::setw(2) << rem_minutes << ":" << std::setw(2) << rem_seconds;
        ss << ", " << std::round(iterations_per_second * 100) / 100 << " iter/s]\r";
 
        std::cout << ss.str();
        std::cout.flush();
    }
};

// Define struct to hold next collision info
struct MinimumInfo{
    int row;
    int col;
    double min_value;
};

// Function to initialise number of particles based on the suggested number, and the quantity rounding adding/sibtracting from this
void initializeNumParts(Params& params) {
    // Define particle size brackets from Raach et al, 10.1089/ast.2016.1544
    std::vector<double> quantity = {0.128315105, 0.183572807, 0.124311959, 0.087425835, 0.062477661, 0.051611981, 0.135821002, 0.208449496, 0.018014154};
    std::vector<double> particle_sizes = {1.0, 3.0, 6.0, 12.0, 23.5, 47.0, 94.0, 187.5, 375.0};
    std::vector<double> ranges = {1.0, 1.0, 2.0, 4.0, 7.5, 16.0, 31.0, 62.5, 125.0};

    //std::vector<double> quantity = {0.5, 0.5};
    //std::vector<double> particle_sizes = {1.0, 100.0};
    //std::vector<double> ranges = {1.0, 10.0};

    // Total is summation of rounded number of particles in each percentage quantity
    int tot = 0;
    for (int i = 0; i < quantity.size(); ++i) {
        tot += static_cast<int>(round(params.init_particles * quantity[i]));
    }

    params.n_particles = tot;
}

void initializeParticleRadii(std::vector<Particle>& parts, Params params) {
    // Define particle size brackets from Raach et al, 10.1089/ast.2016.1544
    std::vector<double> quantity = {0.128315105, 0.183572807, 0.124311959, 0.087425835, 0.062477661, 0.051611981, 0.135821002, 0.208449496, 0.018014154};
    std::vector<double> particle_sizes = {1.0, 3.0, 6.0, 12.0, 23.5, 47.0, 94.0, 187.5, 375.0};
    std::vector<double> ranges = {1.0, 1.0, 2.0, 4.0, 7.5, 16.0, 31.0, 62.5, 125.0};

    //std::vector<double> quantity = {0.5, 0.5};
    //std::vector<double> particle_sizes = {1.0, 100.0};
    //std::vector<double> ranges = {1.0, 10.0};

    // Zip diameter, quantity and range into array
    std::vector<std::vector<double>> psq;
    for (int i = 0; i < quantity.size(); ++i) {
        std::vector<double> row = {particle_sizes[i], round(params.init_particles * quantity[i]), ranges[i]};
        psq.push_back(row);
    }

    // loop through each row in psq and make a vector of radii
    std::vector<double> radii;
    for (const auto& psqr : psq) {
        double size = psqr[0];
        double q = static_cast<int>(psqr[1]);
        double rad_range = psqr[2];

        // Create random device
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(-1, 1);
        
        // calculate radius including random for q elements in each bracket
        for (int j = 0; j < q; ++j) {
            // need 0.5 because Raack data is diameter
            double radius = 0.5 * (size + (dis(gen) * rad_range));
            radii.push_back(radius); 
        }
    }

    // Assign calculated radius to each particle
    for (int i = 0; i < params.n_particles; ++i) {
        parts[i].r = radii[i];
    }
}

// Function to work out how full the box is
double calc_den(std::vector<Particle>& parts, double box, double n){
    // Box volume is cube
    double b_vol = std::pow(box,3);

    // add the spherical volume of each particle
    double p_vol = 0.0;
    for (int i = 0; i < n; ++i){
        p_vol += ((4/3) * M_PI * std::pow(parts[i].r,3));
    }

    // Return ratio of volumes
    return p_vol/b_vol;
}

// Function to work out how size of box for desired packing density
double size_box(std::vector<Particle>& parts, double packing_density, double n){
    // Box volume is cube
    // add the spherical volume of each particle
    double p_vol = 0.0;
    for (int i = 0; i < n; ++i){
        p_vol += ((4/3) * M_PI * std::pow(parts[i].r,3));
    }

    double box = 0.0;
    box = std::pow(p_vol/packing_density,1.0/3.0);
    // Return ratio of volumes
    return box;
}

// Function to check if particles overlap
bool overlap(const Particle& p1, const Particle& p2) {
    // Calculate the euclidean distance between points
    double dx = p1.pos[0] - p2.pos[0];
    double dy = p1.pos[1] - p2.pos[1];
    double dz = p1.pos[2] - p2.pos[2];
    double distance = std::sqrt(dx*dx + dy*dy + dz*dz);

    // Return TRUE if distance smaller than the sum of the radii
    return distance < p1.r + p2.r;
}

// Function to check for overlaps between all particles
bool checkForOverlaps(std::vector<Particle>& parts, int n_particles, Params params) {
    // Only look at upper triangular in a pass
    for (int i = 0; i < n_particles - 1; ++i) {
        for (int j = i + 1; j < n_particles; ++j) {
            // Check if overlap and if so shift particle, breaking at first overlap
            if (overlap(parts[i], parts[j])) {
                if (parts[i].r > parts[j].r){
                    parts[j].init_position(params.box_length);
                }
                else {
                    parts[i].init_position(params.box_length);
                }
                return true; // Overlap found
            }
        }
    }
    return false; // No overlap found
}

// Function to handle elastic collision of particles
void collision(Particle& p1, Particle& p2) {
    // Extract particle velocities and masses
    double u1[3];
    double u2[3];
    double v1[3];
    double v2[3];
    double pos1[3];
    double pos2[3];

    double m1 = p1.mass;
    double m2 = p2.mass;

    for (int i = 0; i < 3; ++i) {
        u1[i] = p1.v[i];
        pos1[i] = p1.pos[i];

        u2[i] = p2.v[i];
        pos2[i] = p2.pos[i];
    }

    // Calculate normal vector

    double n[3];
    double n_mag = std::sqrt(((pos2[0] - pos1[0]) * (pos2[0] - pos1[0])) + ((pos2[1] - pos1[1]) * (pos2[1] - pos1[1])) + ((pos2[2] - pos1[2]) * (pos2[2] - pos1[2])));

    for (int i = 0; i < 3; ++i) {
        n[i] = (pos2[i]-pos1[i]) / n_mag;
    }
    
    // Calculate effective mass
    double m_eff = 1 / ((1/m1) + (1/m2));

    // Calculate impact speed
    double v_imp = 0;

    for (int i = 0; i < 3; ++i) {
        v_imp += n[i] * (u1[i] - u2[i]);
    }

    // Calculate impulse magnitude
    double J;
    double e = 1; //coefficient restitution
    J = (1 + e) * m_eff * v_imp;

    // calcualte deltas
    double delta_v1[3];
    double delta_v2[3];

    for (int i = 0; i < 3; ++i) {
        delta_v1[i] = - (J/m1) * n[i];
        delta_v2[i] = + (J/m2) * n[i];

        v1[i] = u1[i] + delta_v1[i];
        v2[i] = u2[i] + delta_v2[i];

        p1.v[i] = v1[i];
        p2.v[i] = v2[i];
    }

    // // Check KE and momentum
    // double KE_before = 0.0;
    // double KE_after = 0.0;

    // double momentum_before[3] = {0.0};
    // double momentum_after[3] = {0.0};

    // for (int i = 0; i < 3; ++i) {
    //     KE_before += 0.5 * m1 * u1[i] * u1[i] + 0.5 * m2 * u2[i] * u2[i];

    //     momentum_before[i] = m1 * u1[i] + m2 * u2[i];

    //     momentum_after[i] = m1 * v1[i] + m2 * v2[i]; 

    //     KE_after += 0.5 * m1 * v1[i] * v1[i] + 0.5 * m2 * v2[i] * v2[i];
    // }

    // std::cout << "Initial Kinetic Energy: " << KE_before << std::endl;
    // std::cout << "Final Kinetic Energy: " << KE_after << std::endl;

    // // Check momentum conservation
    // std::cout << "Initial Momentum: (" << momentum_before[0] << ", " << momentum_before[1] << ", " << momentum_before[2] << ")" << std::endl;
    // std::cout << "Final Momentum: (" << momentum_after[0] << ", " << momentum_after[1] << ", " << momentum_after[2] << ")" << std::endl;
}

void collision_tom_derivative(Particle& p1, Particle& p2) {
    double u1[3];
    double u2[3];
    double v1[3];
    double v2[3];
    double pos1[3];
    double pos2[3];

    double m1 = p1.mass;
    double m2 = p2.mass;

    for (int i = 0; i < 3; ++i) {
        u1[i] = p1.v[i];
        pos1[i] = p1.pos[i];

        u2[i] = p2.v[i];
        pos2[i] = p2.pos[i];
    }

    double r_rel[3];
    double v_rel[3];
    double r_rel_mag = 0.0;
    double v_dot_r = 0.0;

    for (int i = 0; i < 3; ++i){
        r_rel[i] = pos2[i] - pos1[i];
        v_rel[i] = u2[i] - u1[i];
        
        r_rel_mag += r_rel[i] * r_rel[i];
        v_dot_r += r_rel[i] * v_rel[i];
    }
    
    r_rel_mag = std::sqrt(r_rel_mag);

    double v_rel_al = v_dot_r / r_rel_mag;

    double v_rel_al_aft = (v_rel_al * (m1 - m2)) + 2 * ((m2 * v_rel_al) / (m1 +m2));

    double d_v_rel_al = v_rel_al_aft - v_rel_al;

    for (int i = 0; i < 3; ++i){
        v1[i] = u1[i] + ((d_v_rel_al * r_rel[i]) / r_rel_mag);
        v2[i] = u2[i] - ((d_v_rel_al * r_rel[i]) / r_rel_mag);

        p1.v[i] = v1[i];
        p2.v[i] = v2[i];
    }
    
}

// Function to perform charge exchange in collision
void charge_exchange(Particle& p1, Particle& p2){
    // If p1 has high energy electrons, decrement p1_h and increment p2_l
    if (p1.high_e > 0){
        p1.high_e -= 1;
        p2.low_e += 1;
    }

    // If p2 has high energy electrons, decrement p2_h and increment p1_l
    if (p2.high_e > 0){
        p2.high_e -= 1;
        p1.low_e += 1;
    }
}

// Function to calculate number high energy electrons
int sum_high_e(std::vector<Particle>& parts, int num_parts){
    // Loop through all particles and sum
    int a = 0;
    for (int i = 0; i < num_parts; ++i){
        a += parts[i].high_e;
    }
    return a;
}

// Function to determine the time a given particle will take to hit a given wall
double timeUntilWallCollision(const Particle& p, double box, int wall_id) {
    // Initialise to inf -> i.e. particle doesn't collide
    double t = std::numeric_limits<double>::infinity();

    int ax_id = wall_id % 3;
    // bottom: x -> 0, y -> 1, z -> 2
    // top: x -> 3, y -> 4, z -> 5
    
    // For 'bottom' walls, the particle has to be moving towards them, i.e. negative velocty
    if (wall_id < 3){
        if (p.v[wall_id] < 0.0) {
            t = (p.r - p.pos[wall_id]) / p.v[wall_id];
        }
    }
    // For top walls, to be moving towards them the velocity must be positive
    else if (wall_id < 6 && wall_id >=3){
        
        if (p.v[ax_id] > 0.0) {
            t = (box - p.r - p.pos[ax_id]) / p.v[ax_id];        
        }
    }
    return t;
}

// Function to determine the time a given particle will take to hit a given wall
double timeUntilWallCollision_tom_derative(const Particle& p, double box, int wall_id) {   
    // bottom: x -> 0, y -> 1, z -> 2
    // top: x -> 3, y -> 4, z -> 5
    int ax_id = wall_id % 3;

    int upper_check = wall_id / 3; // 1 or 0

    double dist = (box * upper_check) - p.pos[ax_id] - (p.r * ((2.0 * upper_check) - 1.0));

    //std::cout << upper_check << std::endl;
    
    double t = dist / p.v[ax_id];

    if (t < 0.0){
        t = std::numeric_limits<double>::infinity();
    }

    return t;
}

// Function to determine time until a particle hits another particle
double timeUntilParticleCollision(const Particle& p1, const Particle& p2) {
    // From 10.1140/epje/s10189-022-00180-8

    // Difference in velocity
    double dv[3];
    for (int i = 0; i < 3; ++i) {
        dv[i] = p2.v[i] - p1.v[i];
    }

    // Difference in positions
    double dr[3];
    for (int i = 0; i < 3; ++i) {
        dr[i] = p2.pos[i] - p1.pos[i];
    }

    // Dot product of the two differences
    double dvdr = 0.0;
    for (int i = 0; i < 3; ++i) {
        dvdr += dv[i] * dr[i];
    }

    // Dot product of the velocity difference with itself
    double dv2 = 0.0;
    for (int i = 0; i < 3; ++i) {
        dv2 += dv[i] * dv[i];
    }

    // Dot product of the position difference with itself
    double dr2 = 0.0;
    for (int i = 0; i < 3; ++i) {
        dr2 += dr[i] * dr[i];
    }

    // Check for zero velocity difference
    if (dv2 == 0.0) {
        // Zero velocity difference indicates either parallel movement or one particle stationary
        // They either never collide or are already in collision
        return std::numeric_limits<double>::infinity(); // No collision
    }

    // Calculate discriminant of the time quadratic
    double discriminant = (dvdr * dvdr) - (dv2 * (dr2 - std::pow((p1.r + p2.r), 2)));

    // Collision doesn't occur if the discriminant is negative, or if the two particles are moving apart (dv.dr positive)
    if (discriminant >= 0.0 && dvdr < 0.0) {
        // Avoid division by zero and handle small dr2
        double timeToCollision = (-dvdr - std::sqrt(discriminant)) / dv2;
        if (timeToCollision > 0.0){
            return timeToCollision;
        }
    }

    return std::numeric_limits<double>::infinity(); // No collision
}

// Function to return the initial matrix of collision times
void initializeCollisionTimes(std::vector<std::vector<double>>& collisionTimes, std::vector<Particle>& parts, int n_particles, double box_length) {
    // Calculate collision times between particle pairs, upper triangular only
    for (int i = 0; i < n_particles - 1; ++i) {
        for (int j = i + 1; j < n_particles; ++j) {
            double time_to_collision = timeUntilParticleCollision(parts[i], parts[j]);
            if (time_to_collision > 0.0){
                collisionTimes[i][j] = time_to_collision;
            }
        }
    }

    // Calculate collision times with walls for each particle
    for (int i = 0; i < n_particles; ++i) {
        for (int j = n_particles; j < n_particles + 6; ++j) {
            double time_to_collision = timeUntilWallCollision(parts[i], box_length, j - n_particles);
            if (time_to_collision > 0.0){
                collisionTimes[i][j] = time_to_collision;
            }
        }
    }
}

// Function to sum the kinetic energy of particles
double calculateTotalKineticEnergy(std::vector<Particle>& parts, int numParticles) {
    // Loop through all particles summing
    double totalKE = 0.0;
    for (int i = 0; i < numParticles; ++i) {
        totalKE += parts[i].getKineticEnergy();
    }
    return totalKE;
}

// Function to return the row and column id of the minimum time in the collision time matrix
MinimumInfo findMinimum(const std::vector<std::vector<double>>& matrix){
    // Initialise values
    double min_value = std::numeric_limits<double>::infinity();
    int minRow = -1;
    int minCol = -1;
    // Loop through matrix to find the values when the value of the matrix is lowest
    for (int i = 0; i < matrix.size(); ++i){
        for (int j = 0; j < matrix[i].size(); ++j){
            if (matrix[i][j] < min_value){
                // Set values in the strucutre to the values of the current minimum
                min_value = matrix[i][j];
                minRow = i;
                minCol = j;
            }
        }
    }
    return {minRow, minCol, min_value};
}

// Function to subtract propegated time from the collision time matrix
void subtractPropagationTime(std::vector<std::vector<double>>& collisionTimes, double time) {
    for (int i = 0; i < collisionTimes.size(); ++i) {
        for (int j = 0; j < collisionTimes[i].size(); ++j) {
            collisionTimes[i][j] -= time;
        }
    }
}

// Function to recalculate collision time row for a single colliding particle
void recalculateCollisionTime(std::vector<std::vector<double>>& collisionTimes, std::vector<Particle>& parts, int n_particles, double box_length, int index) {

    // Recalculate collision times between particle pairs keeping upper triangular for simplicity
    for (int i = 0; i < index; ++i) {
        double time_to_collision = timeUntilParticleCollision(parts[i], parts[index]);
        if (time_to_collision > 0.0){
            collisionTimes[i][index] = time_to_collision;
        }
    }

    for (int j = index + 1; j < n_particles; ++j) {
        double time_to_collision = timeUntilParticleCollision(parts[index], parts[j]);
        if (time_to_collision > 0.0){
            collisionTimes[index][j] = time_to_collision;
        }
    }

    // Calculate collision times with walls for each particle
    for (int j = n_particles; j < n_particles + 6; ++j) {
        double time_to_collision = timeUntilWallCollision(parts[index], box_length, j - n_particles);
        if (time_to_collision > 0.0){
            collisionTimes[index][j] = time_to_collision;
        }
    }
}

// Function to handle collision between particle and wall (reflect boundary condition)
void handleWallCollision(Particle& p, double box_length, int wall_id) {
    int axis_id = wall_id % 3; // Determine the axis corresponding to the wall_id
    if (wall_id < 3) {
        // Particle hits the bottom wall
        if (p.pos[axis_id] - p.r <= 0.0 && p.v[axis_id] < 0.0) {
            p.pos[axis_id] = p.r; // Move the particle to the edge
            p.v[axis_id] = std::abs(p.v[axis_id]); // Reflect the velocity
        }
    }
    else{
        // Particle hits the upper wall
        if (p.pos[axis_id] + p.r >= box_length && p.v[axis_id] > 0.0) {
            p.pos[axis_id] = box_length - p.r; // Move the particle to the edge
            p.v[axis_id] = -std::abs(p.v[axis_id]); // Reflect the velocity
        }
    }
}

void handleWallCollision_tom_derivative(Particle& p, double box_length, int wall_id) {
    int axis_id = wall_id % 3; // Determine the axis corresponding to the wall_id
    int upper_check = wall_id / 3; // 1 or 0

    p.pos[axis_id] = (box_length * upper_check) - (p.r * ((2.0 * upper_check) - 1.0));

    p.v[axis_id] = -1.0 * p.v[axis_id];
}

// Function to create output file and write data to it
void writeOutput(const Params& params, const std::vector<Particle>& particles) {

    std::ofstream output_file("lack_model_output.txt");

    // FIO check
    if (!output_file.is_open()) {
        std::cerr << "Unable to open output file!" << std::endl;
        return;
    }
    
    // Parameter information
    output_file << "Parameters:\n";
    output_file << "box_length\t" << params.box_length << "\n";
    output_file << "n_particles\t" << params.n_particles << "\n";
    output_file << "n_eqil_steps\t" << params.n_eqil_steps << "\n";
    output_file << "speed_l\t" << params.speed_l << "\n";
    output_file << "speed_u\t" << params.speed_u << "\n";
    output_file << "density\t" << params.density << "\n";
    output_file << "charge_density\t" << params.charge_density << "\n";

    // Loop through all particles writing their final information
    output_file << "\nradii,masses,position(x),position(y),position(z),velocity(x),velocity(y),velocity(z),initial high e,final high e,final low e,charge\n";
    for (const Particle& p : particles) {
        output_file << p.r << ",";
        output_file << p.mass << ",";
        output_file << p.pos[0] << ",";
        output_file << p.pos[1] << ",";
        output_file << p.pos[2] << ",";
        output_file << p.v[0] << ",";
        output_file << p.v[1] << ",";
        output_file << p.v[2] << ",";
        output_file << p.initial_high_e << ",";
        output_file << p.high_e << ",";
        output_file << p.low_e << ",";
        output_file << p.initial_high_e - p.low_e  << "\n";
    }

    output_file.close();
}

// Function to display a matrix subset like a table
void displayMatrix(const std::vector<std::vector<double>>& matrix) {
    // Loop through matrix printing each element followed by a tab. At the end of the row, new line.
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j < matrix[1].size(); ++j) {
            std::cout << matrix[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

// Define the signal handler function
void signal_callback_handler(int signum){
    writeOutput(params, parts);
    exit(signum);
}

int main(){

    // Retrieve start time
    auto start = std::chrono::system_clock::now();

    // Precalculate number of particles based on distribution and rounding, then allocate object
    initializeNumParts(params);
    parts.resize(params.n_particles);

    // Set up signal handler for SIGINT (Ctrl+C)
    auto signal_handler = [](int signum) {
        // Inside the lambda, call the actual signal handler function
        signal_callback_handler(signum);
    };

    // Register the signal handler
    signal(SIGINT, signal_handler);
    
    // Initialise size of collision matrix, and prefil with inf
    std::vector<std::vector<double>> collisionTimes(params.n_particles, std::vector<double>(params.n_particles + 6, std::numeric_limits<double>::infinity()));

    // Initialise particle radii
    initializeParticleRadii(parts, params);

    params.box_length = size_box(parts, params.packing_density, params.n_particles);


    // Initialise other particle aspects
    for (int i = 0; i < params.n_particles; ++i) {
        //std::cout << parts[i].r << std::endl;
        parts[i].spherical_mass(params.density);
        parts[i].init_high_e(params.charge_density);
        parts[i].low_e = 0;
        parts[i].init_position(params.box_length);
        
        parts[i].init_velocity(params.speed_l, params.speed_u);
        //parts[i].old_init_velocity(params.speed_l, params.speed_u);
        //parts[i].alt_init_velocity(params.temperature);
    }

    // Check box fullness
    std::cout << "Box fullness: " << (calc_den(parts, params.box_length, params.n_particles)) << std::endl;
    
    // Check for overlaps
    int attempts = 0;
    while (checkForOverlaps(parts,params.n_particles,params) && attempts < params.n_max_arrangement){
        attempts++;
    }

    // If still overlaps, terminate simulation, and display message
    if (attempts >= params.n_max_arrangement) {
        std::cout << "Maximum number of attempts reached. Unable to find a non-overlapping configuration." << std::endl;    
    }

    // Otherwise continue with simulation
    else {
        std::cout << "No more overlaps. Initialization complete." << std::endl;
        // Check KE
        std::cout << "KE: " << calculateTotalKineticEnergy(parts, params.n_particles) << std::endl;

        // Initialise event calendar;
        initializeCollisionTimes(collisionTimes,parts,params.n_particles,params.box_length);
        //displayMatrix(collisionTimes);

        // Initialise minimum structure
        MinimumInfo collisionData;

        // Initialise step counter
        int step = 0;

        // Check in equilibriate loop
        std::cout << "Equilibriation Step Start" << std::endl;
        Tqdm ebar(params.n_eqil_steps);
        while (step<params.n_eqil_steps){
            // Find minimum time to event
            collisionData = findMinimum(collisionTimes);

            // Propagate to event
            for (int i =0; i < params.n_particles; ++i){
                parts[i].time_step(collisionData.min_value);
            }

            // Subtract time from event calendar
            subtractPropagationTime(collisionTimes, collisionData.min_value);

            // Handle collision
            if (collisionData.col >= params.n_particles) {
                // Wall collision
                //std::cout << collisionData.row << " " << " wall " << collisionData.col - params.n_particles << std::endl;
                handleWallCollision(parts[collisionData.row], params.box_length, collisionData.col - params.n_particles);

                // Recalculate collision times for particle involved
                recalculateCollisionTime(collisionTimes, parts, params.n_particles, params.box_length, collisionData.row);
            } 
            else{
                // Particle collision
                //std::cout << collisionData.row << " collision " << collisionData.col << std::endl;
                collision(parts[collisionData.row], parts[collisionData.col]);

                // Recalculate collision times for both particles involved
                recalculateCollisionTime(collisionTimes, parts, params.n_particles, params.box_length, collisionData.row);
                recalculateCollisionTime(collisionTimes, parts, params.n_particles, params.box_length, collisionData.col);

            }

            //collisionTimes[collisionData.row][collisionData.col] = std::numeric_limits<double>::infinity();

            // Increment step
            step +=1;

            // Update progress bar
            ebar.update(step, step);
        }
        std::cout << "Equilibriation Step Complete" << std::endl;

        // Initialise total high electron counter
        int total = 0;
        total = sum_high_e(parts, params.n_particles);
        std::cout << total << std::endl;

        // KE check
        std::cout << params.n_particles << std::endl;

        std::cout << "KE: " << calculateTotalKineticEnergy(parts, params.n_particles) << std::endl;

        // Set up progress bar
        int intit_total = total;
        Tqdm pbar(intit_total);

        // Loop while high energy electrons remain
        while (total != 0){        
            // Find Find minimum time to event
            collisionData = findMinimum(collisionTimes);

            // Propagate to event
            for (int i =0; i < params.n_particles; ++i){
                parts[i].time_step(collisionData.min_value);
            }

            // Subtract time from event calendar
            subtractPropagationTime(collisionTimes, collisionData.min_value);

            // Handle collision
            if (collisionData.col > params.n_particles-1) {
                // Wall collision
                // std::cout << collisionData.row << " " << " wall " << collisionData.col - params.n_particles << std::endl;
                handleWallCollision(parts[collisionData.row], params.box_length, collisionData.col - params.n_particles);

                // Recalculate collision times for particle involved
                recalculateCollisionTime(collisionTimes, parts, params.n_particles, params.box_length, collisionData.row);

            }
            else{
                // Particle collision
                // std::cout << collisionData.row << " collision " << collisionData.col << std::endl;
                collision(parts[collisionData.row], parts[collisionData.col]);

                // Exchange charges
                charge_exchange(parts[collisionData.row], parts[collisionData.col]);

                // Recalculate collision times for both particles involved
                recalculateCollisionTime(collisionTimes,parts,params.n_particles,params.box_length,collisionData.row);
                recalculateCollisionTime(collisionTimes,parts,params.n_particles,params.box_length,collisionData.col);

            }

            //collisionTimes[collisionData.row][collisionData.col] = std::numeric_limits<double>::infinity();

            step +=1;

            // Update total high energy electrons
            total = sum_high_e(parts, params.n_particles);

            // Update progress bar
            pbar.update(intit_total - total,step - params.n_eqil_steps);
        }

        // Write output
        writeOutput(params, parts);
    }

    // Final timing handling
    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "finished computation at " << std::ctime(&end_time) << "elapsed time: " << elapsed_seconds.count() << "s" << std::endl;
}