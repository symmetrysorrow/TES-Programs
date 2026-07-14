#include "tes_cpp/pulse.hpp"

#include <exception>
#include <iostream>
#include <sstream>

int main(int argc, char** argv) {
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
