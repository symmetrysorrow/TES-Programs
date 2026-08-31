#include "tes_cpp/converter.hpp"

#include <exception>
#include <iostream>
#include <string>
#include <filesystem>

namespace {
void usage() {
    std::cerr << "Usage: dump2event <dumpall.dat> <event.json> --input-energy <MeV>\n"
                 "Normal options:\n"
                 "  --save-all                     retain all histories\n"
                 "  --full-energy-only             retain only histories without NCOL=12 leakage\n"
                 "  --full-energy-list <path>      write full-energy event IDs\n"
                 "Diagnostic/development options:\n"
                 "  --history-summary <path>       write per-history diagnostics\n"
                 "  --max-histories <N>             stop after N source histories\n"
                 "  --summary-only                 do not write an event file\n";
}
}

int main(int argc, char** argv) {
    if (argc < 5) { usage(); return 2; }
    tes_cpp::dump2event::Options options;
    std::string full_energy_list;
    std::string history_summary;
    bool summary_only = false;
    for (int i = 3; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--input-energy" && i + 1 < argc) options.input_energy = std::stod(argv[++i]);
        else if (arg == "--save-all") options.save_all = true;
        else if (arg == "--full-energy-only") options.full_energy_only = true;
        else if (arg == "--full-energy-list" && i + 1 < argc) full_energy_list = argv[++i];
        else if (arg == "--history-summary" && i + 1 < argc) history_summary = argv[++i];
        else if (arg == "--max-histories" && i + 1 < argc) options.max_histories = std::stoi(argv[++i]);
        else if (arg == "--summary-only") summary_only = true;
        else { usage(); return 2; }
    }
    if (options.full_energy_only) options.save_all = false;
    options.collect_full_energy_ids = !full_energy_list.empty();
    if (summary_only && history_summary.empty()) { usage(); return 2; }
    try {
        const auto result = tes_cpp::dump2event::read_dump(argv[1], options);
        if (!summary_only) {
            if (std::filesystem::path(argv[2]).extension() == ".h5" || std::filesystem::path(argv[2]).extension() == ".hdf5")
                tes_cpp::dump2event::write_event_hdf5(result.batch, argv[2]);
            else tes_cpp::dump2event::write_event_json(result.batch, argv[2]);
        }
        if (!full_energy_list.empty()) tes_cpp::dump2event::write_event_ids(result.full_energy_event_ids, full_energy_list);
        if (!history_summary.empty()) tes_cpp::dump2event::write_history_summary(result.history_summaries, history_summary);
        if (summary_only)
            std::cout << "Wrote " << result.history_summaries.size() << " history summaries to " << history_summary << '\n';
        else
            std::cout << "Wrote " << result.batch.size() << " events to " << argv[2] << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dump2event: " << error.what() << '\n';
        return 1;
    }
}
