#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace tes_cpp {

struct Pulse {
    int position{};
    std::vector<double> time;
    std::vector<double> ch0;
    std::vector<double> ch1;
};

// Positions are one-based absorber-block indices: 1 <= position <= n_abs.
std::vector<Pulse> generate_pulses(
    const std::string& input_json_path,
    const std::vector<int>& positions);

void write_pulses_json(
    const std::vector<Pulse>& pulses,
    const std::string& input_json_path,
    const std::string& output_path);

// Generate pulses directly into one JSON file.  Waveforms are kept only while
// their worker is computing them; completed waveforms are staged in temporary
// files and merged in the same order as ``positions``.
// A value of zero uses std::thread::hardware_concurrency().
void generate_pulses_json(
    const std::string& input_json_path,
    const std::vector<int>& positions,
    const std::string& output_path,
    std::size_t thread_count = 0);

}  // namespace tes_cpp
