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

// Kahan summation device function
__device__ void kahan_add(double& sum, double& compensation, double value) {
    double y = value - compensation;
    double t = sum + y;
    compensation = (t - sum) - y;
    sum = t;
}


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
    double* __restrict__ block_results,  // Changed: per-block results
    int num_particles,
    double center_x,
    double center_y,
    double v_x,
    double v_y)
{
    // Physical constants with high precision
    const double mu_0 = 1.25663706212e-6;
    const double eps_0 = 8.8541878128e-12;
    const double q_e = 1.602176634e-19;
    const double pi = 3.14159265358979323846;
    
    const double mu_0_qe_div_4pi = (mu_0 * q_e) / (4.0 * pi);
    const double qe_div_4pi_eps0 = q_e / (4.0 * pi * eps_0);
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    
    // Kahan summation variables for each thread
    double sum_B[3] = {0.0, 0.0, 0.0};
    double comp_B[3] = {0.0, 0.0, 0.0};
    double sum_E[3] = {0.0, 0.0, 0.0};
    double comp_E[3] = {0.0, 0.0, 0.0};
    
    // Process particles assigned to this thread
    for (int i = idx; i < num_particles; i += stride) {
        const double dx = point_x - particle_x[i];
        const double dy = point_y - particle_y[i];
        const double dz = point_z - particle_z[i];
        
        const double r_sq = dx*dx + dy*dy + dz*dz;
        
        const double r_mag = sqrt(r_sq);
        const double r3_inv = 1.0 / (r_sq * r_mag);
        
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
        
        // Use Kahan summation for better precision
        kahan_add(sum_B[0], comp_B[0], B_factor * Bx);
        kahan_add(sum_B[1], comp_B[1], B_factor * By);
        kahan_add(sum_B[2], comp_B[2], B_factor * Bz);
        kahan_add(sum_E[0], comp_E[0], E_factor * dx);
        kahan_add(sum_E[1], comp_E[1], E_factor * dy);
        kahan_add(sum_E[2], comp_E[2], E_factor * dz);
    }
    
    // Block-level reduction with Kahan summation
    __shared__ double sdata[6][256];
    __shared__ double scomp[6][256];
    
    int tid = threadIdx.x;
    
    // Store thread results in shared memory
    if (tid < 256) {
        sdata[0][tid] = sum_B[0];
        sdata[1][tid] = sum_B[1];
        sdata[2][tid] = sum_B[2];
        sdata[3][tid] = sum_E[0];
        sdata[4][tid] = sum_E[1];
        sdata[5][tid] = sum_E[2];
        
        // Initialize compensation
        for (int i = 0; i < 6; ++i) {
            scomp[i][tid] = comp_B[i % 3];
            if (i >= 3) scomp[i][tid] = comp_E[i - 3];
        }
    }
    
    __syncthreads();
    
    // Kahan reduction in shared memory
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s && tid + s < 256) {
            for (int i = 0; i < 6; ++i) {
                kahan_add(sdata[i][tid], scomp[i][tid], sdata[i][tid + s]);
            }
        }
        __syncthreads();
    }
    
    // Write block results instead of atomic add
    if (tid == 0) {
        int block_idx = blockIdx.x;
        for (int i = 0; i < 6; ++i) {
            block_results[i * gridDim.x + block_idx] = sdata[i][0];
        }
    }
}
// CUDA kernel for particle timestep update (for convergence study)
// Replace lines 208-230 with this improved kernel:
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
        
        // Get current position relative to current center
        double dx = particle_x[i] - center_x;
        double dy = particle_y[i] - center_y;

        double theta = particle_wz[i] * dt;
        double cos_theta = cos(theta);
        double sin_theta = sin(theta);

        double rotated_dx = dx * cos_theta - dy * sin_theta;
        double rotated_dy = dx * sin_theta + dy * cos_theta;

        double new_center_x = center_x + velocity_x * dt;
        double new_center_y = center_y + velocity_y * dt;

        // Final position = new center + rotated relative position
        particle_x[i] = new_center_x + rotated_dx;
        particle_y[i] = new_center_y + rotated_dy;
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
    double *d_block_results;  // Add this for per-block results
    int max_particles;
    int max_blocks;
    
public:
    ConvergenceSolver(int max_n_particles) : max_particles(max_n_particles) {
        gpu_particles = std::make_unique<ParticleDataGPU>(max_n_particles);
        CUDA_CHECK(cudaMalloc(&d_result, 6 * sizeof(double)));
        
        // Allocate for per-block results
        max_blocks = (max_n_particles + 255) / 256;
        CUDA_CHECK(cudaMalloc(&d_block_results, 6 * max_blocks * sizeof(double)));
    }
    
    ~ConvergenceSolver() {
        cudaFree(d_result);
        cudaFree(d_block_results);
    }
    
    void setParticles(const std::vector<double>& h_x, const std::vector<double>& h_y, const std::vector<double>& h_z,
                  const std::vector<double>& h_q, const std::vector<double>& h_wx, const std::vector<double>& h_wy, 
                  const std::vector<double>& h_wz) {
    
    int n = h_x.size();
    size_t size = n * sizeof(double);
    
    // Copy existing data
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
        
        // Clear results
        CUDA_CHECK(cudaMemset(d_result, 0, 6 * sizeof(double)));
        CUDA_CHECK(cudaMemset(d_block_results, 0, 6 * max_blocks * sizeof(double)));
        
        int threads_per_block = 256;
        int blocks = (num_particles + threads_per_block - 1) / threads_per_block;
        
        // Launch kernel with block results
        compute_single_point_field_kernel<<<blocks, threads_per_block>>>(
            gpu_particles->x, gpu_particles->y, gpu_particles->z,
            gpu_particles->q, gpu_particles->wx, gpu_particles->wy, gpu_particles->wz,
            point_x, point_y, point_z,
            d_block_results,  // Pass block results buffer
            num_particles, center_x, center_y, v_x, v_y);
        
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
        
        // Copy block results to CPU
        std::vector<double> h_block_results(6 * blocks);
        CUDA_CHECK(cudaMemcpy(h_block_results.data(), d_block_results, 
                             6 * blocks * sizeof(double), cudaMemcpyDeviceToHost));
        
        // High-precision CPU summation using Kahan summation
        std::vector<double> final_result(6, 0.0);
        std::vector<double> compensation(6, 0.0);
        
        for (int b = 0; b < blocks; ++b) {
            for (int component = 0; component < 6; ++component) {
                double value = h_block_results[component * blocks + b];
                
                // Kahan summation on CPU
                double y = value - compensation[component];
                double t = final_result[component] + y;
                compensation[component] = (t - final_result[component]) - y;
                final_result[component] = t;
            }
        }
        
        return final_result;
    }
    
     void updateParticles(int num_particles, double dt, double center_x, double center_y, double velocity_x, double velocity_y) {
        int threads_per_block = 256;
        int blocks = (num_particles + threads_per_block - 1) / threads_per_block;
        
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

void writeTimeseriesFile(const ConvergenceResult& result) {
    std::string filename = "convergence_" + std::to_string(result.particle_count) + "_particles_timeseries.csv";
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing" << std::endl;
        return;
    }
    
    // Write header
    file << "timestep,time,run,Bx,By,Bz,Ex,Ey,Ez" << std::endl;
    
    // Write data - one row per timestep per run
    for (size_t t = 0; t < result.timesteps.size(); ++t) {
        double time_seconds = result.timesteps[t] * 0.5; // dt = 0.5
        
        for (int run = 0; run < result.Bx_values[t].size(); ++run) {
            file << std::scientific << std::setprecision(12);
            file << result.timesteps[t] << "," << time_seconds << "," << (run + 1);
            file << "," << result.Bx_values[t][run] << "," << result.By_values[t][run] << "," << result.Bz_values[t][run];
            file << "," << result.Ex_values[t][run] << "," << result.Ey_values[t][run] << "," << result.Ez_values[t][run];
            file << std::endl;
        }
    }
    
    file.close();
    std::cout << "  Timeseries saved to: " << filename << std::endl;
}

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
        
        std::cout << "Testing " << test_timesteps.size() << " time points over " << T_end << " seconds" << std::endl;
        
        // Initialize solver with maximum particle count
        ConvergenceSolver solver(*std::max_element(particle_counts.begin(), particle_counts.end()));
        
        std::vector<ConvergenceResult> results;
        
        for (int N_Sim : particle_counts) {
            std::cout << "\nTesting " << N_Sim << " particles (" << num_runs << " runs)..." << std::endl;
            
            size_t free_mem, total_mem;
            CUDA_CHECK(cudaMemGetInfo(&free_mem, &total_mem));
            std::cout << "GPU Memory - Total: " << total_mem / (1024*1024) << " MB, " 
                    << "Free: " << free_mem / (1024*1024) << " MB" << std::endl;

            // Add this before each particle count test:
            std::cout << "Required memory for " << N_Sim << " particles: " 
                    << (N_Sim * 7 * sizeof(double)) / (1024*1024) << " MB" << std::endl;

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
                
                // EFFICIENT: Single pass through all timesteps
                double current_center_x = center_x;
                double current_center_y = center_y;
                size_t next_target_idx = 0;  // Index of next target timestep to record
                
                for (int timestep = 0; timestep <= timesteps; ++timestep) {
                    // Check if this is a target timestep
                    if (next_target_idx < test_timesteps.size() && timestep == test_timesteps[next_target_idx]) {
                        // Compute field at this target timestep
                        auto field = solver.computeSinglePointField(N_Sim, test_x, test_y, test_z, 
                                                                current_center_x, current_center_y, velocity_x, velocity_y);
                        
                        result.Bx_values[next_target_idx].push_back(field[0]);
                        result.By_values[next_target_idx].push_back(field[1]);
                        result.Bz_values[next_target_idx].push_back(field[2]);
                        result.Ex_values[next_target_idx].push_back(field[3]);
                        result.Ey_values[next_target_idx].push_back(field[4]);
                        result.Ez_values[next_target_idx].push_back(field[5]);
                        
                        next_target_idx++;
                    }
                    
                    // Advance particles to next timestep
                    if (timestep < timesteps) {
                        solver.updateParticles(N_Sim, dt, current_center_x, current_center_y, velocity_x, velocity_y);
                        current_center_x += velocity_x * dt;
                        current_center_y += velocity_y * dt;
                    }
                }
                
                if ((run + 1) % 100 == 0) {
                    std::cout << "  Completed " << run + 1 << "/" << num_runs << " runs" << std::endl;
                }
            }
            
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
            
            writeTimeseriesFile(result);

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
