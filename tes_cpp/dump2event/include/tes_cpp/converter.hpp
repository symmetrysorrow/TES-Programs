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

struct HistorySummary {
    int event_id{};
    int particle_records{};
    bool had_secondary_particle{};
    int secondary_particles{};
    int secondary_photons{};
    int secondary_electrons{};
    int energy_cutoff_photons{};
    int energy_cutoff_electrons{};
    int energy_cutoff_particles{};
    int leaked_particles{};
    int leaked_photons{};
    int leaked_primary_photons{};
    int leaked_secondary_photons{};
    double leaked_energy{};
    int reaction_records{};
    double total_deposit{};
    bool energy_consistent{};
    bool full_energy{};
    bool fully_contained{};
    bool stored_event{};
};

struct Options {
    // Unit must match the energy recorded in dumpall.dat (MeV in the current app).
    double input_energy{};
    // false: retain only full-energy events. true: retain every event and report
    // the full-energy event IDs separately.
    bool save_all{false};
    // Explicit full-energy-only mode.  With save_all=false this is already
    // the effective retention policy, but the flag makes the initial
    // full-energy analysis explicit at the CLI/API.
    bool full_energy_only{false};
    // Ask the parser to return the IDs of full-energy histories even when
    // only those histories are retained.
    bool collect_full_energy_ids{false};
    // Stop after this many source histories when running a diagnostic pass.
    // Zero means read the complete dump file.
    int max_histories{0};
    // Tolerance used for the reconstructed-deposit diagnostic. Full-energy
    // classification is based on the absence of leakage terminations.
    double energy_tolerance{0.005};
    // Secondary charged particles below this energy are treated as locally
    // deposited.  The current PHITS/EGS5 dump files contain electron tracks
    // down to 1 keV, so the default must match that transport cut-off.
    double electron_min_energy{0.001};
    double photon_min_energy{0.001};
};

struct Result {
    Batch batch;
    std::vector<int> full_energy_event_ids;
    std::vector<HistorySummary> history_summaries;
};

// Parses a PHITS dumpall.dat file. Throws std::runtime_error on invalid input or I/O.
Result read_dump(const std::string& dump_path, const Options& options);

// Writes the stable event.json schema used by the existing application.
void write_event_json(const Batch& batch, const std::string& output_path);
void write_event_hdf5(const Batch& batch, const std::string& output_path);
void write_event_ids(const std::vector<int>& event_ids, const std::string& output_path);
void write_history_summary(const std::vector<HistorySummary>& summaries, const std::string& output_path);

}  // namespace tes_cpp::dump2event
