#include "tes_cpp/pulse.hpp"

#include <exception>
#include <iomanip>
#include <iostream>
#include <sstream>

int main(int argc, char** argv) {
    if (argc == 3 && std::string(argv[1]) == "--dump-linearization") {
        try {
            const tes_cpp::LinearizationSummary summary = tes_cpp::inspect_linearization(argv[2]);
            std::cout << std::setprecision(17)
                      << "{\"current_A\":" << summary.current_A
                      << ",\"tau_el_s\":" << summary.tau_el_s
                      << ",\"loop_gain\":" << summary.loop_gain
                      << ",\"tau_i_s\":" << summary.tau_i_s
                      << ",\"g_abs_tes_W_per_K\":" << summary.g_abs_tes_W_per_K
                      << ",\"tes_boundary_rate_per_s\":" << summary.tes_boundary_rate_per_s
                      << ",\"tes_hanging_rate_per_s\":" << summary.tes_hanging_rate_per_s
                      << ",\"tes_intrinsic_thermal_diag_per_s\":" << summary.tes_intrinsic_thermal_diag_per_s
                      << ",\"n_abs\":" << summary.n_abs
                      << ",\"hanging\":" << (summary.hanging ? "true" : "false")
                      << ",\"tes1_time_block\":[";
            for (std::size_t i = 0; i < summary.tes1_time_block.size(); ++i) {
                if (i) std::cout << ',';
                std::cout << summary.tes1_time_block[i];
            }
            std::cout << "],\"tes2_time_block\":[";
            for (std::size_t i = 0; i < summary.tes2_time_block.size(); ++i) {
                if (i) std::cout << ',';
                std::cout << summary.tes2_time_block[i];
            }
            std::cout << "]}\n";
            return 0;
        } catch (const std::exception& error) {
            std::cerr << "posi2pulse: " << error.what() << '\n';
            return 1;
        }
    }
    if ((argc != 5 && argc != 7) || std::string(argv[3]) != "--positions" ||
        (argc == 7 && std::string(argv[5]) != "--threads")) {
        std::cerr << "Usage: posi2pulse <input.json> <pulses.json> --positions 1,2,3 [--threads N]\n";
        return 2;
    }
    std::vector<int> positions;
    std::stringstream items(argv[4]);
    std::string item;
    while (std::getline(items, item, ',')) {
        if (item.empty()) { std::cerr << "empty position\n"; return 2; }
        positions.push_back(std::stoi(item));
    }
    try {
        std::size_t threads = 0;
        if (argc == 7) {
            const long long value = std::stoll(argv[6]);
            if (value < 1) throw std::runtime_error("threads must be at least 1");
            threads = static_cast<std::size_t>(value);
        }
        tes_cpp::generate_pulses_json(argv[1], positions, argv[2], threads);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "posi2pulse: " << error.what() << '\n';
        return 1;
    }
}
