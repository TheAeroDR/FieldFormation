#ifndef SENSITIVITY_PARAMETERS_H
#define SENSITIVITY_PARAMETERS_H

#include <string>

// Structure to hold all sensitivity parameters
struct SensitivityParams {
    // Bias factors (1.0 = no bias, 1.1 = +10%, 0.9 = -10%)
    double rad_profile_bias = 1.0;      // dd_rad_profile bias (inner=0.9, outer=1.1)
    double theta_bias = 1.0;            // Initial theta bias
    double height_bias = 1.0;           // rad2height bias
    double charge_bias = 1.0;           // Particle charge bias
    double electron_density_bias = 1.0; // Electron density bias
    double eyewall_velocity_bias = 1.0; // V0 bias
    double dust_loading_bias = 1.0;     // dust_loading bias
    double dust_density_bias = 1.0;     // dust_rho bias
    double avg_particle_size_bias = 1.0;// Average particle size bias
    double translational_velocity_bias = 1.0; // Translational velocity bias
    double diameter_bias = 1.0;         // Diameter bias (for particle size)
    double particles_per_height_bias = 1.0; // Bias for number of particles per height
    double height_median_bias = 1/3; // Bias for height median
    
    // Default constructor
    SensitivityParams() = default;
    
    // Constructor for specific sensitivity test
    SensitivityParams(const std::string& test_name) {
        if (test_name == "rad_profile_inner_10") {
            rad_profile_bias = 0.9;
        }
        else if (test_name == "rad_profile_outer_10") {
            rad_profile_bias = 1.1;
        }
        else if (test_name == "theta_plus_10") {
            theta_bias = 1.1;
        }
        else if (test_name == "theta_minus_10") {
            theta_bias = 0.9;
        }
        else if (test_name == "height_plus_10") {
            height_bias = 1.1;
        }
        else if (test_name == "height_minus_10") {
            height_bias = 0.9;
        }
        else if (test_name == "charge_plus_10") {
            charge_bias = 1.1;
        }
        else if (test_name == "charge_minus_10") {
            charge_bias = 0.9;
        }
        else if (test_name == "electron_density_plus_10") {
            electron_density_bias = 1.1;
        }
        else if (test_name == "electron_density_minus_10") {
            electron_density_bias = 0.9;
        }
        else if (test_name == "eyewall_velocity_plus_10") {
            eyewall_velocity_bias = 1.1;
        }
        else if (test_name == "eyewall_velocity_minus_10") {
            eyewall_velocity_bias = 0.9;
        }
        else if (test_name == "dust_loading_plus_10") {
            dust_loading_bias = 1.1;
        }
        else if (test_name == "dust_loading_minus_10") {
            dust_loading_bias = 0.9;
        }
        else if (test_name == "dust_density_plus_10") {
            dust_density_bias = 1.1;
        }
        else if (test_name == "dust_density_minus_10") {
            dust_density_bias = 0.9;
        }
        else if (test_name == "avg_particle_size_plus_10") {
            avg_particle_size_bias = 1.1;
        }
        else if (test_name == "avg_particle_size_minus_10") {
            avg_particle_size_bias = 0.9;
        }
        else if (test_name == "translational_velocity_plus_10") {
            translational_velocity_bias = 1.1;
        }
        else if (test_name == "translational_velocity_minus_10") {
            translational_velocity_bias = 0.9;
        }
        else if (test_name == "diameter_plus_10") {
            diameter_bias = 1.1;
        }
        else if (test_name == "diameter_minus_10") {
            diameter_bias = 0.9;
        }
        else if (test_name == "particles_per_height_plus") {
            particles_per_height_bias = 1.25;
        }
        else if (test_name == "particles_per_height_minus") {
            particles_per_height_bias = 0.8;
        }
        else if (test_name == "height_median_plus") {
            height_median_bias = 0.5;
        }
        else if (test_name == "height_median_minus") {
            height_median_bias = 0.1;
        }
        else {
            // Default to no bias if test name is unknown
            rad_profile_bias = 1.0;
            theta_bias = 1.0;
            height_bias = 1.0;
            charge_bias = 1.0;
            electron_density_bias = 1.0;
            eyewall_velocity_bias = 1.0;
            dust_loading_bias = 1.0;
            dust_density_bias = 1.0;
            avg_particle_size_bias = 1.0;
            translational_velocity_bias = 1.0;
            diameter_bias = 1.0;
            particles_per_height_bias = 1.0;
        }
    }
    
    std::string getDescription() const {
        std::string desc = "Sensitivity Parameters:\n";
        if (rad_profile_bias != 1.0) desc += "  Radial profile bias: " + std::to_string(rad_profile_bias) + "\n";
        if (theta_bias != 1.0) desc += "  Theta bias: " + std::to_string(theta_bias) + "\n";
        if (height_bias != 1.0) desc += "  Height bias: " + std::to_string(height_bias) + "\n";
        if (charge_bias != 1.0) desc += "  Charge bias: " + std::to_string(charge_bias) + "\n";
        if (electron_density_bias != 1.0) desc += "  Electron density bias: " + std::to_string(electron_density_bias) + "\n";
        if (eyewall_velocity_bias != 1.0) desc += "  Eyewall velocity bias: " + std::to_string(eyewall_velocity_bias) + "\n";
        if (dust_loading_bias != 1.0) desc += "  Dust loading bias: " + std::to_string(dust_loading_bias) + "\n";
        if (dust_density_bias != 1.0) desc += "  Dust density bias: " + std::to_string(dust_density_bias) + "\n";
        if (avg_particle_size_bias != 1.0) desc += "  Avg particle size bias: " + std::to_string(avg_particle_size_bias) + "\n";
        if (translational_velocity_bias != 1.0) desc += "  Translational velocity bias: " + std::to_string(translational_velocity_bias) + "\n";
        if (diameter_bias != 1.0) desc += "  Diameter bias: " + std::to_string(diameter_bias) + "\n";
        if (particles_per_height_bias != 1.0) desc += "  Particles per height bias: " + std::to_string(particles_per_height_bias) + "\n";
        if (height_median_bias != 1.0) desc += "  Height median bias: " + std::to_string(height_median_bias) + "\n";
        return desc;
    }
};

#endif // SENSITIVITY_PARAMETERS_H
