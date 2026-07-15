#pragma once

#include <map>
#include <string>
#include <vector>

namespace tes_cpp::dump2event {

struct EventInfo {
    int ityp{};
    std::vector<double> x, y, z, energy;
    std::vector<double> x_deposit, y_deposit, z_deposit, energy_deposit;
};

using History = std::map<int, EventInfo>;
using Batch = std::map<int, History>;

struct Options {
    // Unit must match the energy recorded in dumpall.dat (MeV in the current app).
    double input_energy{};
    // false: retain only full-energy events. true: retain every event and report
    // the full-energy event IDs separately.
    bool save_all{false};
    double energy_tolerance{0.001};
};

struct Result {
    Batch batch;
    std::vector<int> full_energy_event_ids;
};

// Parses a PHITS dumpall.dat file. Throws std::runtime_error on invalid input or I/O.
Result read_dump(const std::string& dump_path, const Options& options);

// Writes the stable event.json schema used by the existing application.
void write_event_json(const Batch& batch, const std::string& output_path);
void write_event_hdf5(const Batch& batch, const std::string& output_path);
void write_event_ids(const std::vector<int>& event_ids, const std::string& output_path);

}  // namespace tes_cpp::dump2event
