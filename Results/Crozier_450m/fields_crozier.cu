// CUDA version for electromagnetic field computation
// to compile: nvcc -O3 -arch=sm_89 -std=c++17 fields_cuda.cu -o fields_cuda.exe
// Note: RTX 4000 Ada uses compute capability 8.9 (sm_89)

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

// Updated CUDA kernel for field computation with Kahan summation
__global__ void compute_fields_kernel(
    const double* __restrict__ particle_x,
    const double* __restrict__ particle_y,
    const double* __restrict__ particle_z,
    const double* __restrict__ particle_q,
    const double* __restrict__ particle_wx,
    const double* __restrict__ particle_wy,
    const double* __restrict__ particle_wz,
    const double* __restrict__ grid_x,
    const double* __restrict__ grid_y,
    const double* __restrict__ grid_z,
    double* __restrict__ output_data,
    int num_particles,
    int num_grid_points,
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
    
    int grid_idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (grid_idx >= num_grid_points) return;
    
    const double px = grid_x[grid_idx];
    const double py = grid_y[grid_idx];
    const double pz = grid_z[grid_idx];
    
    // Kahan summation variables
    double sum_B[3] = {0.0, 0.0, 0.0};
    double comp_B[3] = {0.0, 0.0, 0.0};
    double sum_E[3] = {0.0, 0.0, 0.0};
    double comp_E[3] = {0.0, 0.0, 0.0};
    
    // Each thread processes all particles for one grid point
    for (int i = 0; i < num_particles; ++i) {
        const double dx = px - particle_x[i];
        const double dy = py - particle_y[i];
        const double dz = pz - particle_z[i];
        
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
    
    // Store results (x, y, z, Bx, By, Bz, Ex, Ey, Ez)
    int base_idx = grid_idx * 9;
    output_data[base_idx + 0] = px;
    output_data[base_idx + 1] = py;
    output_data[base_idx + 2] = pz;
    output_data[base_idx + 3] = sum_B[0];
    output_data[base_idx + 4] = sum_B[1];
    output_data[base_idx + 5] = sum_B[2];
    output_data[base_idx + 6] = sum_E[0];
    output_data[base_idx + 7] = sum_E[1];
    output_data[base_idx + 8] = sum_E[2];
}

// Updated CUDA kernel for particle timestep update
__global__ void update_particles_kernel(
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

// Host wrapper class
class CUDAFieldSolver {
private:
    std::unique_ptr<ParticleDataGPU> gpu_particles;
    double *d_grid_x, *d_grid_y, *d_grid_z;
    double *d_output;
    int num_particles;
    int max_grid_points;
    
public:
    CUDAFieldSolver(int n_particles, int max_grid) 
        : num_particles(n_particles), max_grid_points(max_grid) {
        
        // RTX 4000 Ada has 20GB VRAM - can handle full grid without batching
        int actual_gpu_grid_size = max_grid;
        
        // Allocate GPU memory
        gpu_particles = std::make_unique<ParticleDataGPU>(n_particles);
        
        CUDA_CHECK(cudaMalloc(&d_grid_x, actual_gpu_grid_size * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_grid_y, actual_gpu_grid_size * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_grid_z, actual_gpu_grid_size * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_output, actual_gpu_grid_size * 9 * sizeof(double)));
        
        std::cout << "CUDA initialized for RTX 4000 Ada (20GB VRAM): " << actual_gpu_grid_size << " grid points" << std::endl;
    }
    
    ~CUDAFieldSolver() {
        cudaFree(d_grid_x);
        cudaFree(d_grid_y);
        cudaFree(d_grid_z);
        cudaFree(d_output);
    }
    
    void setParticles(const std::vector<double>& h_x, const std::vector<double>& h_y, const std::vector<double>& h_z,
                      const std::vector<double>& h_q, const std::vector<double>& h_wx, const std::vector<double>& h_wy, 
                      const std::vector<double>& h_wz) {
        
        gpu_particles->copyFromCPU(h_x, h_y, h_z, h_q, h_wx, h_wy, h_wz);
    }
    
    void computeFields(const std::vector<double>& h_grid_x, const std::vector<double>& h_grid_y, 
                       const std::vector<double>& h_grid_z, std::vector<double>& h_output,
                       double center_x, double center_y, double v_x, double v_y) {
        
        int num_grid = h_grid_x.size();
        
        // RTX 4000 Ada with 20GB VRAM - no batching needed for 1M+ grid points
        int max_grid_per_batch = max_grid_points;
        
        if (num_grid > max_grid_per_batch) {
            std::cout << "Large grid (" << num_grid << " points) - using batched computation..." << std::endl;
            
            h_output.resize(num_grid * 9);
            int num_batches = (num_grid + max_grid_per_batch - 1) / max_grid_per_batch;
            
            for (int batch = 0; batch < num_batches; ++batch) {
                int start_idx = batch * max_grid_per_batch;
                int end_idx = std::min(start_idx + max_grid_per_batch, num_grid);
                int batch_size = end_idx - start_idx;
                
                std::cout << "  Batch " << batch + 1 << "/" << num_batches 
                          << " (" << batch_size << " points)" << std::endl;
                
                // Create batch vectors
                std::vector<double> batch_grid_x(h_grid_x.begin() + start_idx, h_grid_x.begin() + end_idx);
                std::vector<double> batch_grid_y(h_grid_y.begin() + start_idx, h_grid_y.begin() + end_idx);
                std::vector<double> batch_grid_z(h_grid_z.begin() + start_idx, h_grid_z.begin() + end_idx);
                
                // Copy batch grid points to GPU
                CUDA_CHECK(cudaMemcpy(d_grid_x, batch_grid_x.data(), batch_size * sizeof(double), cudaMemcpyHostToDevice));
                CUDA_CHECK(cudaMemcpy(d_grid_y, batch_grid_y.data(), batch_size * sizeof(double), cudaMemcpyHostToDevice));
                CUDA_CHECK(cudaMemcpy(d_grid_z, batch_grid_z.data(), batch_size * sizeof(double), cudaMemcpyHostToDevice));
                
                // Launch kernel for this batch
                int threads_per_block = 256;
                int blocks = (batch_size + threads_per_block - 1) / threads_per_block;
                
                compute_fields_kernel<<<blocks, threads_per_block>>>(
                    gpu_particles->x, gpu_particles->y, gpu_particles->z,
                    gpu_particles->q, gpu_particles->wx, gpu_particles->wy, gpu_particles->wz,
                    d_grid_x, d_grid_y, d_grid_z,
                    d_output,
                    num_particles, batch_size, center_x, center_y, v_x, v_y);
                
                CUDA_CHECK(cudaGetLastError());
                CUDA_CHECK(cudaDeviceSynchronize());
                
                // Copy batch results back
                CUDA_CHECK(cudaMemcpy(h_output.data() + start_idx * 9, d_output, batch_size * 9 * sizeof(double), cudaMemcpyDeviceToHost));
            }
        } else {
            // Original single-batch computation for smaller grids
            CUDA_CHECK(cudaMemcpy(d_grid_x, h_grid_x.data(), num_grid * sizeof(double), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_grid_y, h_grid_y.data(), num_grid * sizeof(double), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_grid_z, h_grid_z.data(), num_grid * sizeof(double), cudaMemcpyHostToDevice));
            
            int threads_per_block = 256;
            int blocks = (num_grid + threads_per_block - 1) / threads_per_block;
            
            compute_fields_kernel<<<blocks, threads_per_block>>>(
                gpu_particles->x, gpu_particles->y, gpu_particles->z,
                gpu_particles->q, gpu_particles->wx, gpu_particles->wy, gpu_particles->wz,
                d_grid_x, d_grid_y, d_grid_z,
                d_output,
                num_particles, num_grid, center_x, center_y, v_x, v_y);
            
            CUDA_CHECK(cudaGetLastError());
            CUDA_CHECK(cudaDeviceSynchronize());
            
            h_output.resize(num_grid * 9);
            CUDA_CHECK(cudaMemcpy(h_output.data(), d_output, num_grid * 9 * sizeof(double), cudaMemcpyDeviceToHost));
        }
    }
    
    void updateParticles(double dt, double center_x, double center_y, double velocity_x, double velocity_y) {
        int threads_per_block = 256;
        int blocks = (num_particles + threads_per_block - 1) / threads_per_block;
        
        update_particles_kernel<<<blocks, threads_per_block>>>(
            gpu_particles->x, gpu_particles->y, gpu_particles->wz,
            dt, center_x, center_y, velocity_x, velocity_y,
            num_particles);
        
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
    }
    
    void getParticlePositions(std::vector<double>& h_x, std::vector<double>& h_y) {
        h_x.resize(num_particles);
        h_y.resize(num_particles);
        CUDA_CHECK(cudaMemcpy(h_x.data(), gpu_particles->x, num_particles * sizeof(double), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(h_y.data(), gpu_particles->y, num_particles * sizeof(double), cudaMemcpyDeviceToHost));
    }
};

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
        else if(ed == 1.8){
            q = 22.55 * std::pow(r, 2) - 5.656e4 * std::pow(r, -0.8);
        }
        else if(ed == 2){
            q = 25.07 * std::pow(r, 2) - 5.464e4 * std::pow(r, -0.8);
        }
        else if(ed == 2.2){
            q = 27.56 * std::pow(r, 2) - 8.025e4 * std::pow(r, -0.8);
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
    double K = 0.0;
    int switching = 3; // 0 for Vasitas, 1 for Rankine, 2 for Burgers, 3 for Hollow core
    if (switching == 0){
        K = PhysicalConstants::pi / 2.0; // Vasitas vortex
    }
    else if (switching == 1){
        K = 1.0; // Rankine vortex
    }
    else if (switching == 2){
        K = 1.702; // Burgers vortex
    }
    else if (switching == 3){
        K = 0.5; // Hollow core vortex
    }
    else{
        std::cerr << "Invalid switching value for velocity profile." << std::endl;
    }
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

// Generate a truly random seed combining multiple sources
unsigned int generateRandomSeed() {
    // Combine multiple sources of randomness
    std::random_device rd;
    auto time_seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    auto steady_seed = std::chrono::steady_clock::now().time_since_epoch().count();
    
    // Mix the seeds using a simple hash-like operation
    unsigned int mixed_seed = static_cast<unsigned int>(
        rd() ^ 
        (time_seed & 0xFFFFFFFF) ^ 
        ((steady_seed >> 32) & 0xFFFFFFFF) ^
        (steady_seed & 0xFFFFFFFF)
    );
    
    return mixed_seed;
}

std::vector<SimpleParticle> initializeParticles(double D, int N_Sim, int N_AtHeight, double centre_x, double centre_y, unsigned int seed = std::random_device{}()) {
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

    std::cout << "CUDA Simulation: " << clump_factor << " clump factor, " << max_height << " max height, seed: " << seed << std::endl;

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
            p.x = centre_x + p.R * std::cos(theta_offset);
            p.y = centre_y + p.R * std::sin(theta_offset);
            p.q = -1.0 * clump_factor * p.charge_profile();
            
            p.wx = 0.0;
            p.wy = 0.0;
        }
    }

    return particles;
}

void writeCSV(const std::vector<double>& output, int num_points, const std::string& filename) {
    std::ofstream file(filename);
    file.precision(12);
    file << std::scientific;
    file << "x,y,z,Bx,By,Bz,Ex,Ey,Ez\n";
    
    for (int i = 0; i < num_points; ++i) {
        int base = i * 9;
        file << output[base + 0] << "," << output[base + 1] << "," << output[base + 2] << ","
             << output[base + 3] << "," << output[base + 4] << "," << output[base + 5] << ","
             << output[base + 6] << "," << output[base + 7] << "," << output[base + 8] << "\n";
    }
}

int main(int argc, char* argv[]) {
    try {
        // Parse command line arguments for seed
        unsigned int seed = std::random_device{}(); // default seed
        if (argc > 1) {
            seed = static_cast<unsigned int>(std::stoul(argv[1]));
        }
        
        // Check CUDA device
        int device_count;
        CUDA_CHECK(cudaGetDeviceCount(&device_count));
        if (device_count == 0) {
            std::cerr << "No CUDA devices found!" << std::endl;
            return 1;
        }
        CUDA_CHECK(cudaSetDevice(0));
        
        // Create output directory
        std::string output_dir = "Crozier_2e_CUDA";
        std::filesystem::create_directory(output_dir);
        
        // Problem parameters
        double D = 20.0;
        int N_Sim = 10000000;
        int N_AtHeight = 100;

        // Simulation parameters
        double dt = 1.0;
        double T_end = 1.0;
        int timesteps = static_cast<int>(T_end / dt);
        double center_x = 0.0, center_y = 0.0;
        double velocity_x = -4.00, velocity_y = 0.0;
        
        // Grid parameters
        double grid_spacing = 0.2;
        int grid_size = 82;
        
        std::cout << "CUDA Simulation: " << N_Sim << " particles, " << grid_size << "x" << 2 << " grid" << std::endl;
        
        // Use seed from command line argument
        auto particles = initializeParticles(D, N_Sim, N_AtHeight, center_x, center_y, seed);
        
        // Save particles to CSV
        std::ofstream particle_file(output_dir + "/particles.csv");
        particle_file << "x,y,z,r,q,R,Vt,wz\n";
        for (const auto& p : particles) {
            particle_file << p.x << "," << p.y << "," << p.z << ","
                          << p.r << "," << p.q << "," << p.R << ","
                          << p.Vt << "," << p.wz << "\n";
        }
        particle_file.close();
        std::cout << "Particles saved to: " << output_dir << "/particles.csv" << std::endl;

        // Create grid points
        std::vector<double> grid_x, grid_y, grid_z;
        
        for (int i = 0; i < grid_size; ++i) {
        //for (int i = 0; i < int(floor(grid_size/2))+1; ++i) {
            //double y = -500.0 + i * grid_spacing; // generate grid from -500 to 500
            double z = 0.0 + i * grid_spacing; // generate grid from 0 to 40
            double x = -450.0; 
                double y = 0.0; // fixed y position
                grid_x.push_back(x);
                grid_y.push_back(y);
                grid_z.push_back(z);
                for (int j = 0; j < 2; ++j) {
                    
                }
        }
        // Convert particles to GPU format
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
        
        // Initialize CUDA solver
        CUDAFieldSolver solver(N_Sim, grid_x.size());
        solver.setParticles(h_x, h_y, h_z, h_q, h_wx, h_wy, h_wz);
        
        
        std::cout << "Starting " << timesteps << " timesteps..." << std::endl;
        
        for (int t = 0; t < timesteps; ++t) {
            auto start_time = std::chrono::high_resolution_clock::now();
            
            std::vector<double> output;
            solver.computeFields(grid_x, grid_y, grid_z, output, center_x, center_y, velocity_x, velocity_y);
            
            auto end_time = std::chrono::high_resolution_clock::now();
            
            // Write output file
            std::ostringstream filename;
            filename << output_dir << "/field_t" << std::setw(5) << std::setfill('0') << t << ".csv";
            writeCSV(output, grid_x.size(), filename.str());
            
            // Update particles
            solver.updateParticles(dt, center_x, center_y, velocity_x, velocity_y);
            
            // Update dust devil center
            center_x += velocity_x * dt;
            center_y += velocity_y * dt;
            
            // Progress reporting
            if ((t + 1) % 10 == 0 || t == 0) {
                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
                std::cout << "Step " << t+1 << "/" << timesteps << ": " << duration.count() << "ms" << std::endl;
            }
        }
        
        std::cout << "CUDA simulation completed successfully!" << std::endl;
        std::cout << "Output files written to: " << output_dir << "/" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
