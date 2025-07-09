// CUDA Convergence Study for electromagnetic field computation
// Tests particle counts from 10^4 to 10^7 with 1000 runs each
// Only computes fields at single point (0,0,1)
// Compile: nvcc -O3 -arch=sm_89 -std=c++17 convergence_cuda.cu -o convergence_cuda.exe

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cmath>
#include <random>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <chrono>

// Physical constants
namespace PhysicalConstants {
    constexpr double pi = 3.14159265358979323846;
}

// Error checking macro
#define CUDA_CHECK(call) \
    do { \
        cudaError_t error = call; \
        if (error != cudaSuccess) { \
            std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__ << " - " << cudaGetErrorString(error) << std::endl; \
            exit(1); \
        } \
    } while(0)

// GPU particle data structure
struct ParticleDataGPU {
    double *x, *y, *z;
    double *q;
    double *wx, *wy, *wz;
    int num_particles;
    
    ParticleDataGPU(int n) : num_particles(n) {
        size_t size = n * sizeof(double);
        CUDA_CHECK(cudaMalloc(&x, size));
        CUDA_CHECK(cudaMalloc(&y, size));
        CUDA_CHECK(cudaMalloc(&z, size));
        CUDA_CHECK(cudaMalloc(&q, size));
        CUDA_CHECK(cudaMalloc(&wx, size));
        CUDA_CHECK(cudaMalloc(&wy, size));
        CUDA_CHECK(cudaMalloc(&wz, size));
    }
    
    ~ParticleDataGPU() {
        cudaFree(x);
        cudaFree(y);
        cudaFree(z);
        cudaFree(q);
        cudaFree(wx);
        cudaFree(wy);
        cudaFree(wz);
    }
    
    void copyFromCPU(const std::vector<double>& h_x, const std::vector<double>& h_y, const std::vector<double>& h_z,
                     const std::vector<double>& h_q, const std::vector<double>& h_wx, const std::vector<double>& h_wy, const std::vector<double>& h_wz) {
        size_t size = num_particles * sizeof(double);
        CUDA_CHECK(cudaMemcpy(x, h_x.data(), size, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(y, h_y.data(), size, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(z, h_z.data(), size, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(q, h_q.data(), size, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(wx, h_wx.data(), size, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(wy, h_wy.data(), size, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(wz, h_wz.data(), size, cudaMemcpyHostToDevice));
    }
};

// CUDA kernel for single point field computation
__global__ void compute_single_point_field_kernel(
    const double* __restrict__ particle_x,
    const double* __restrict__ particle_y,
    const double* __restrict__ particle_z,
    const double* __restrict__ particle_q,
    const double* __restrict__ particle_wx,
    const double* __restrict__ particle_wy,
    const double* __restrict__ particle_wz,
    double point_x,
    double point_y,
    double point_z,
    double* __restrict__ result,  // [Bx, By, Bz, Ex, Ey, Ez]
    int num_particles,
    double center_x,
    double center_y,
    double v_x,
    double v_y)
{
    // Constants
    const double mu_0 = 1.2566370621219e-6;
    const double eps_0 = 8.854187812813e-12;
    const double q_e = 1.60217663e-19;
    const double pi = 3.14159265358979323846;
    
    const double mu_0_qe_div_4pi = (mu_0 * q_e) / (4.0 * pi);
    const double qe_div_4pi_eps0 = q_e / (4.0 * pi * eps_0);
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    
    // Shared memory for reduction (use larger data type for better precision)
    __shared__ double sdata[6][256]; // Bx, By, Bz, Ex, Ey, Ez
    
    int tid = threadIdx.x;
    
    // Initialize shared memory
    for (int i = 0; i < 6; ++i) {
        sdata[i][tid] = 0.0;
    }
    
    // Each thread processes multiple particles
    for (int i = idx; i < num_particles; i += stride) {
        const double dx = point_x - particle_x[i];
        const double dy = point_y - particle_y[i];
        const double dz = point_z - particle_z[i];
        
        const double r_sq = dx*dx + dy*dy + dz*dz;
        const double r_mag = sqrt(r_sq);
        const double r3_inv = 1.0 / fmax(r_sq * r_mag, 1e-12);
        
        // Magnetic field calculation
        const double r_x = particle_x[i] - center_x;
        const double r_y = particle_y[i] - center_y;
        const double r_z = particle_z[i] - 0.0;
        
        const double wxr_x = particle_wy[i] * r_z - particle_wz[i] * r_y;
        const double wxr_y = particle_wz[i] * r_x - particle_wx[i] * r_z;
        const double wxr_z = particle_wx[i] * r_y - particle_wy[i] * r_x;

        const double v_total_x = wxr_x + v_x;
        const double v_total_y = wxr_y + v_y;
        const double v_total_z = wxr_z;
        
        const double Bx = v_total_y * dz - v_total_z * dy;
        const double By = v_total_z * dx - v_total_x * dz;
        const double Bz = v_total_x * dy - v_total_y * dx;
        
        const double B_factor = mu_0_qe_div_4pi * particle_q[i] * r3_inv;
        const double E_factor = qe_div_4pi_eps0 * particle_q[i] * r3_inv;
        
        sdata[0][tid] += B_factor * Bx;
        sdata[1][tid] += B_factor * By;
        sdata[2][tid] += B_factor * Bz;
        sdata[3][tid] += E_factor * dx;
        sdata[4][tid] += E_factor * dy;
        sdata[5][tid] += E_factor * dz;
    }
    
    __syncthreads();
    
    // Reduction in shared memory
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            for (int i = 0; i < 6; ++i) {
                sdata[i][tid] += sdata[i][tid + s];
            }
        }
        __syncthreads();
    }
    
    // Write result for this block using higher precision accumulation
    if (tid == 0) {
        for (int i = 0; i < 6; ++i) {
            // Use a more precise atomic add for very small numbers
            double val = sdata[i][0];
            if (val != 0.0) {  // Only add non-zero contributions
                atomicAdd(&result[i], val);
            }
        }
    }
}

// CUDA kernel for particle timestep update (for convergence study)
__global__ void update_particles_convergence_kernel(
    double* __restrict__ particle_x,
    double* __restrict__ particle_y,
    const double* __restrict__ particle_wz,
    double dt,
    double center_x,
    double center_y,
    double velocity_x,
    double velocity_y,
    int num_particles)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i >= num_particles) return;
    
    // Rotational update using angular velocity directly
    const double theta = particle_wz[i] * dt;  // wz is already negative for clockwise
    const double cos_theta = cos(theta);
    const double sin_theta = sin(theta);
    
    const double dx = particle_x[i] - center_x;
    const double dy = particle_y[i] - center_y;
    
    particle_x[i] = center_x + (dx * cos_theta - dy * sin_theta);
    particle_y[i] = center_y + (dx * sin_theta + dy * cos_theta);
    
    // Translational update
    particle_x[i] += velocity_x * dt;
    particle_y[i] += velocity_y * dt;
}

// PDF class for grain diameter distribution
class CustomPDF {
private:
    double h_t;
    double a;
    double b;
    mutable std::mt19937 rng;

    double heaviside(double x, double y) const {
        return 0.5 * std::tanh(x - y) + 0.5;
    }

    std::vector<double> compute_pdf(const std::vector<double>& d_array, double h) const {
        double mu1 = 4.7849;
        double sigma1 = 0.7022;
        double mu2 = 1.3386;
        double sigma2 = 0.9418;
        
        std::vector<double> result(d_array.size());
        
        for (size_t i = 0; i < d_array.size(); ++i) {
            double log_d = std::log(d_array[i]);
            
            double term1 = (1.0 / (sigma1 * std::sqrt(2 * PhysicalConstants::pi))) *
                          std::exp(-((log_d - mu1) * (log_d - mu1)) / (2 * sigma1 * sigma1));
            
            double term2 = (1.0 / (sigma2 * std::sqrt(2 * PhysicalConstants::pi))) *
                          std::exp(-((log_d - mu2) * (log_d - mu2)) / (2 * sigma2 * sigma2));
            
            double H = heaviside(h, h_t);
            result[i] = (term1 * (1 - H) + term2 * H) / d_array[i];
        }
        return result;
    }

public:
    CustomPDF(double h_t_, double a_, double b_, unsigned int seed = std::random_device{}())
        : h_t(h_t_), a(a_), b(b_), rng(seed) {
    }

    std::vector<double> rvs(int size, double h) const {
        int batch_size = std::max(10000, size * 10);
        double max_pdf = 100.0;
        double envelope = max_pdf * (1.0 / (b - a));

        std::uniform_real_distribution<double> uniform_a_b(a, b);
        std::uniform_real_distribution<double> uniform_0_1(0.0, 1.0);

        std::vector<double> samples;
        samples.reserve(size);

        while ((int)samples.size() < size) {
            std::vector<double> y_batch(batch_size);
            std::vector<double> u_batch(batch_size);
            for (int i = 0; i < batch_size; ++i) {
                y_batch[i] = uniform_a_b(rng);
                u_batch[i] = uniform_0_1(rng);
            }

            std::vector<double> pdf_vals = compute_pdf(y_batch, h);
            
            for (int i = 0; i < batch_size && (int)samples.size() < size; ++i) {
                if (u_batch[i] < (pdf_vals[i] / envelope)) {
                    samples.push_back(y_batch[i]);
                }
            }
        }

        return samples;
    }
};

// Particle structure
struct SimpleParticle {
    double r, q, R, Vt;
    double x, y, z, wx, wy, wz;
    
    double charge_profile() {
        double q;
        double ed = 2;
        if(ed == 1){
            q = 12.53 * std::pow(r, 2) - 3.73e4 * std::pow(r, -0.8);
        }
        else if(ed == 0.1){
            q = 1.253 * std::pow(r, 2) - 3952 * std::pow(r, -0.8);
        }
        else if(ed == 0.2){
            q = 2.507 * std::pow(r, 2) - 7101 * std::pow(r, -0.8);
        }
        else if(ed == 0.3){
            q = 3.759 * std::pow(r, 2) - 1.185e4 * std::pow(r, -0.8);
        }
        else if(ed == 0.4){
            q = 5.009 * std::pow(r, 2) - 1.76e4 * std::pow(r, -0.8);
        }
        else if(ed == 0.5){
            q = 6.272 * std::pow(r, 2) - 2.027e4 * std::pow(r, -0.8);
        }
        else if(ed == 0.45){
            q = 5.643 * std::pow(r, 2) - 1.936e4 * std::pow(r, -0.8);
        }
        else if(ed == 2){
            q = 25.07 * std::pow(r, 2) - 5.464e4 * std::pow(r, -0.8);
        }
        else if(ed == 2.5){
            q = 31.33 * std::pow(r, 2) - 6.76e4 * std::pow(r, -0.8);
        }
        else if(ed == 3){
            q = 37.6 * std::pow(r, 2) - 1.039e5 * std::pow(r, -0.8);
        }
        else if(ed == 10){
            q = 125.3 * std::pow(r, 2) - 1.585e5 * std::pow(r, -0.8);
        }
        else if(ed == 20){
            q = 250.7 * std::pow(r, 2) - 5.748e5 * std::pow(r, -0.8);
        }
        else{
            std::cerr << "Invalid ed value for charge profile." << std::endl;
            q = 0.0;
        }
        return q;
    }
    
    double dd_vel_profile(double V0, double D) {
        int switching = 3; // Hollow core
        
        if (switching == 0){
            double n = 2.0;
            double a = 2.0;
            double rbar = R/(D/4.0);
            return V0 * (rbar/pow(std::pow(rbar,a*n)+1.0,1/n));
        }
        else if (switching == 1){
            if (R<=(D/2.0)){
                return V0 * (R/(D/4.0));
            }
            else{
                return V0 * (D/4.0)/(R);
            }
        }
        else if (switching == 2){
            double rbar = R/(D/4.0);
            return V0 / rbar * (1.0 - std::exp(-std::pow(rbar, 2.0)));
        }
        else if (switching == 3){
            if (R<(D/4.0)){
                return 0;
            }
            else{
                return V0 * (D/4.0)/(R);
            }
        }
        else{
            std::cerr << "Invalid switching value for velocity profile." << std::endl;
            return 0.0;
        }
    }
};

// Helper functions
double rad2height(double D) {
    return std::pow((D / 2.0) * (1.0 / 0.913), 2.0);
}

double dd_volume(double D, double H) {
    double outer_radius = D / 2.0;
    double inner_radius = D / 4.0;
    return H * PhysicalConstants::pi * (std::pow(outer_radius, 2.0) - std::pow(inner_radius, 2.0));
}

double tan_vel_profile(double H) {
    double g0 = 9.81;
    double Ro = 0.375;
    double deltaT = 8.0;
    double T0 = 305;
    double b0 = (g0 * deltaT)/T0;
    double K = 0.5; // Hollow core vortex
    double coeff = b0/(K-(pow(Ro,2.0)/2.0));
    return pow(coeff * H,0.5);
}

std::vector<double> grain_diameter(double h, int nparts) {
    CustomPDF custom_dist(2, 0.0001, 1000.0);
    return custom_dist.rvs(nparts, h);
}

double dd_rad_profile(double D) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dis(0.0, 1.0);

    double rand_val = dis(gen);
    return (D / 2.0) * (1.0 - 0.5 * rand_val);
}

std::vector<double> thirdbiased(int N, double min, double max) {
    double power = std::log(1.0/3.0)/ std::log(0.5);
    std::vector<double> t(N);
    for (int i = 0; i < N; ++i) {
        t[i] = std::pow((double)i / (N - 1), power);
    }
    for (int i = 0; i < N; ++i) {
        t[i] = min + (max - min) * t[i];
    }
    return t;
}

std::vector<SimpleParticle> initializeParticles(double D, int N_Sim, int N_AtHeight, unsigned int seed) {
    const double dust_rho = 2.0e6;
    const double dust_loading = 0.296;

    double max_height = rad2height(D);
    double volume = dd_volume(D, max_height);
    double V0 = tan_vel_profile(max_height);
    double dust_mass = dust_loading * volume;
    double av_mass = dust_rho * (4.0 / 3.0) * PhysicalConstants::pi * std::pow(3e-6, 3.0);

    double N = std::ceil(dust_mass / av_mass);
    double clump_factor = N / static_cast<double>(N_Sim);

    int N_heights = N_Sim / N_AtHeight;

    std::vector<double> heights = thirdbiased(N_heights, 0.0, max_height);
    std::vector<SimpleParticle> particles(N_Sim);
    std::mt19937 gen(seed);
    std::uniform_real_distribution<> angle_dist(0.0, 2.0 * PhysicalConstants::pi);

    for (int i = 0; i < N_heights; ++i) {
        std::vector<double> diameters = grain_diameter(heights[i], N_AtHeight);

        for (int j = 0; j < N_AtHeight; ++j) {
            int idx = i * N_AtHeight + j;
            if (idx >= N_Sim) break;
            
            SimpleParticle& p = particles[idx];

            p.z = heights[i];
            p.r = diameters[j] / 2.0;
            p.R = dd_rad_profile(D);
            p.Vt = p.dd_vel_profile(V0, D);
            p.wz = -1.0 * (p.Vt / p.R);

            double theta_offset = angle_dist(gen);
            p.x = -102.0 + p.R * std::cos(theta_offset);
            p.y = 0.0 + p.R * std::sin(theta_offset);
            p.q = -1.0 * clump_factor * p.charge_profile();
            
            p.wx = 0.0;
            p.wy = 0.0;
        }
    }

    return particles;
}

class ConvergenceSolver {
private:
    std::unique_ptr<ParticleDataGPU> gpu_particles;
    double *d_result;
    int max_particles;
    
public:
    ConvergenceSolver(int max_n_particles) : max_particles(max_n_particles) {
        gpu_particles = std::make_unique<ParticleDataGPU>(max_n_particles);
        CUDA_CHECK(cudaMalloc(&d_result, 6 * sizeof(double))); // Bx, By, Bz, Ex, Ey, Ez
    }
    
    ~ConvergenceSolver() {
        cudaFree(d_result);
    }
    
    void setParticles(const std::vector<double>& h_x, const std::vector<double>& h_y, const std::vector<double>& h_z,
                      const std::vector<double>& h_q, const std::vector<double>& h_wx, const std::vector<double>& h_wy, 
                      const std::vector<double>& h_wz) {
        
        int n = h_x.size();
        size_t size = n * sizeof(double);
        CUDA_CHECK(cudaMemcpy(gpu_particles->x, h_x.data(), size, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(gpu_particles->y, h_y.data(), size, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(gpu_particles->z, h_z.data(), size, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(gpu_particles->q, h_q.data(), size, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(gpu_particles->wx, h_wx.data(), size, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(gpu_particles->wy, h_wy.data(), size, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(gpu_particles->wz, h_wz.data(), size, cudaMemcpyHostToDevice));
    }
    
    std::vector<double> computeSinglePointField(int num_particles, double point_x, double point_y, double point_z,
                                               double center_x, double center_y, double v_x, double v_y) {
        
        // Clear result
        CUDA_CHECK(cudaMemset(d_result, 0, 6 * sizeof(double)));
        
        int threads_per_block = 256;
        int blocks = std::min(128, (num_particles + threads_per_block - 1) / threads_per_block);  // Reduced blocks for better precision
        
        compute_single_point_field_kernel<<<blocks, threads_per_block>>>(
            gpu_particles->x, gpu_particles->y, gpu_particles->z,
            gpu_particles->q, gpu_particles->wx, gpu_particles->wy, gpu_particles->wz,
            point_x, point_y, point_z,
            d_result,
            num_particles, center_x, center_y, v_x, v_y);
        
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
        
        std::vector<double> h_result(6);
        CUDA_CHECK(cudaMemcpy(h_result.data(), d_result, 6 * sizeof(double), cudaMemcpyDeviceToHost));
        
        return h_result;
    }
    
    void updateParticles(int num_particles, double dt, double center_x, double center_y, double velocity_x, double velocity_y) {
        int threads_per_block = 256;
        int blocks = std::min(256, (num_particles + threads_per_block - 1) / threads_per_block);
        
        update_particles_convergence_kernel<<<blocks, threads_per_block>>>(
            gpu_particles->x, gpu_particles->y, gpu_particles->wz,
            dt, center_x, center_y, velocity_x, velocity_y, num_particles);
        
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
    }
};

struct ConvergenceResult {
    int particle_count;
    std::vector<int> timesteps;  // Which timesteps were tested
    std::vector<std::vector<double>> Bx_values, By_values, Bz_values;  // [timestep][run]
    std::vector<std::vector<double>> Ex_values, Ey_values, Ez_values;  // [timestep][run]
    std::vector<double> Bx_mean, By_mean, Bz_mean;  // Mean for each timestep
    std::vector<double> Ex_mean, Ey_mean, Ez_mean;  // Mean for each timestep
    std::vector<double> Bx_std, By_std, Bz_std;    // Std for each timestep
    std::vector<double> Ex_std, Ey_std, Ez_std;    // Std for each timestep
    
    void initializeForTimesteps(const std::vector<int>& test_timesteps, int num_runs) {
        timesteps = test_timesteps;
        int num_times = test_timesteps.size();
        
        Bx_values.resize(num_times);
        By_values.resize(num_times);
        Bz_values.resize(num_times);
        Ex_values.resize(num_times);
        Ey_values.resize(num_times);
        Ez_values.resize(num_times);
        
        for (int t = 0; t < num_times; ++t) {
            Bx_values[t].reserve(num_runs);
            By_values[t].reserve(num_runs);
            Bz_values[t].reserve(num_runs);
            Ex_values[t].reserve(num_runs);
            Ey_values[t].reserve(num_runs);
            Ez_values[t].reserve(num_runs);
        }
        
        Bx_mean.resize(num_times);
        By_mean.resize(num_times);
        Bz_mean.resize(num_times);
        Ex_mean.resize(num_times);
        Ey_mean.resize(num_times);
        Ez_mean.resize(num_times);
        Bx_std.resize(num_times);
        By_std.resize(num_times);
        Bz_std.resize(num_times);
        Ex_std.resize(num_times);
        Ey_std.resize(num_times);
        Ez_std.resize(num_times);
    }

};

int main() {
    try {
        // Check CUDA device
        int device_count;
        CUDA_CHECK(cudaGetDeviceCount(&device_count));
        if (device_count == 0) {
            std::cerr << "No CUDA devices found!" << std::endl;
            return 1;
        }
        CUDA_CHECK(cudaSetDevice(0));
        
        std::cout << "CUDA Convergence Study Starting..." << std::endl;
        
        // Parameters
        double D = 7.0;
        int N_AtHeight = 100;
        double center_x = -102.0, center_y = 0.0;
        double velocity_x = 5.38, velocity_y = 0.0;
        
        // Test point (0, 0, 1)
        double test_x = 0.0, test_y = 0.0, test_z = 1.0;
        
        // Particle counts to test: 10^3, 10^4, 10^5, 10^6, 10^7
        std::vector<int> particle_counts = {1000, 10000, 100000, 1000000, 10000000};
        int num_runs = 1000;
        
        // Time simulation parameters
        double dt = 5e-1;        // 0.5 second timesteps
        double T_end = 38.0;     // 38 seconds total
        int timesteps = static_cast<int>(T_end / dt);  // 76 timesteps
        
        // Test times (sample every 5 timesteps to reduce data volume)
        std::vector<int> test_timesteps;
        for (int t = 0; t < timesteps; t += 5) {
            test_timesteps.push_back(t);
        }
        test_timesteps.push_back(timesteps - 1); // Always include final timestep
        
        std::cout << "Testing " << test_timesteps.size() << " time points over " << T_end << " seconds" << std::endl;
        
        // Initialize solver with maximum particle count
        ConvergenceSolver solver(*std::max_element(particle_counts.begin(), particle_counts.end()));
        
        std::vector<ConvergenceResult> results;
        
        for (int N_Sim : particle_counts) {
            std::cout << "\nTesting " << N_Sim << " particles (" << num_runs << " runs)..." << std::endl;
            
            ConvergenceResult result;
            result.particle_count = N_Sim;
            result.initializeForTimesteps(test_timesteps, num_runs);
            
            auto start_time = std::chrono::high_resolution_clock::now();
            
            for (int run = 0; run < num_runs; ++run) {
                // Initialize particles with different random seed for each run
                auto particles = initializeParticles(D, N_Sim, N_AtHeight, run + 1234);
                
                // Convert to GPU format
                std::vector<double> h_x, h_y, h_z, h_q, h_wx, h_wy, h_wz;
                for (const auto& p : particles) {
                    h_x.push_back(p.x);
                    h_y.push_back(p.y);
                    h_z.push_back(p.z);
                    h_q.push_back(p.q);
                    h_wx.push_back(p.wx);
                    h_wy.push_back(p.wy);
                    h_wz.push_back(p.wz);
                }
                
                solver.setParticles(h_x, h_y, h_z, h_q, h_wx, h_wy, h_wz);
                
                // Simulate through time and record fields at test timesteps
                double current_center_x = center_x;
                double current_center_y = center_y;
                int current_timestep = 0;  // Track current simulation time for this run
                
                for (size_t t_idx = 0; t_idx < test_timesteps.size(); ++t_idx) {
                    int target_timestep = test_timesteps[t_idx];
                    
                    // Advance particles from current_timestep to target_timestep
                    int steps_to_advance = target_timestep - current_timestep;
                    for (int step = 0; step < steps_to_advance; ++step) {
                        solver.updateParticles(N_Sim, dt, current_center_x, current_center_y, velocity_x, velocity_y);
                        current_center_x += velocity_x * dt;
                        current_center_y += velocity_y * dt;
                        current_timestep++;
                    }
                    
                    // Compute field at test point for this timestep
                    auto field = solver.computeSinglePointField(N_Sim, test_x, test_y, test_z, 
                                                               current_center_x, current_center_y, velocity_x, velocity_y);
                    
                    result.Bx_values[t_idx].push_back(field[0]);
                    result.By_values[t_idx].push_back(field[1]);
                    result.Bz_values[t_idx].push_back(field[2]);
                    result.Ex_values[t_idx].push_back(field[3]);
                    result.Ey_values[t_idx].push_back(field[4]);
                    result.Ez_values[t_idx].push_back(field[5]);
                }
                
                if ((run + 1) % 100 == 0) {
                    std::cout << "  Completed " << run + 1 << "/" << num_runs << " runs" << std::endl;
                }
            }
            
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
            
            results.push_back(result);
            
            std::cout << "  Completed in " << duration.count() << " seconds" << std::endl;
        }
        
        // Final summary message
        std::cout << "\nConvergence study completed!" << std::endl;
             
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
