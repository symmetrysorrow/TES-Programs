#include "tes_cpp/converter.hpp"

#include <exception>
#include <iostream>
#include <string>
#include <filesystem>

namespace {
void usage() {
    std::cerr << "Usage: dump2event <dumpall.dat> <event.json> --input-energy <MeV> [--save-all] [--full-energy-list <path>]\n";
}
}

int main(int argc, char** argv) {
    if (argc < 5) { usage(); return 2; }
    tes_cpp::dump2event::Options options;
    std::string full_energy_list;
    for (int i = 3; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--input-energy" && i + 1 < argc) options.input_energy = std::stod(argv[++i]);
        else if (arg == "--save-all") options.save_all = true;
        else if (arg == "--full-energy-list" && i + 1 < argc) full_energy_list = argv[++i];
        else { usage(); return 2; }
    }
    try {
        const auto result = tes_cpp::dump2event::read_dump(argv[1], options);
        if (std::filesystem::path(argv[2]).extension() == ".h5" || std::filesystem::path(argv[2]).extension() == ".hdf5")
            tes_cpp::dump2event::write_event_hdf5(result.batch, argv[2]);
        else tes_cpp::dump2event::write_event_json(result.batch, argv[2]);
        if (!full_energy_list.empty()) tes_cpp::dump2event::write_event_ids(result.full_energy_event_ids, full_energy_list);
        std::cout << "Wrote " << result.batch.size() << " events to " << argv[2] << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dump2event: " << error.what() << '\n';
        return 1;
    }
}
