#include "tes_cpp/dumpall_record_parser.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace tes_cpp::dump2event;

int main(int argc, char** argv) try {
    if (argc != 4) {
        std::cerr << "Usage: dumpall_transport_export <ASCII dumpall> <events.csv> <tracks.csv>\n";
        return 2;
    }
    std::ifstream input(argv[1], std::ios::binary);
    std::ofstream events(argv[2], std::ios::binary);
    std::ofstream tracks(argv[3], std::ios::binary);
    if (!input || !events || !tracks) throw std::runtime_error("cannot open input/output");
    events << "event_id\n";
    tracks << "event_id,particle_id,ityp,ncol,source_record_id,x_cm,y_cm,z_cm,energy_MeV,time,region_id\n"
           << std::setprecision(17);
    int event_id = -1;
    const auto summary = parse_dumpall_records(input, [&](const DumpallRecord& record) {
        if (record.ncol == 4) {
            event_id = record.nocas;
            events << event_id << '\n';
        } else if (record.ncol >= 10 && record.ncol <= 15) {
            if (event_id < 0) throw std::runtime_error("transport record precedes NCOL=4 history header");
            tracks << event_id << ',' << record.particle_id << ',' << record.ityp << ',' << record.ncol << ','
                   << record.line_start << ',' << record.x << ',' << record.y << ',' << record.z << ','
                   << record.energy << ',' << record.time << ',' << record.region_before << '\n';
        }
        if (!events || !tracks) throw std::runtime_error("failed writing streaming transport export");
    }, argv[1]);
    if (summary.terminal_ncol != 2) throw std::runtime_error("Dumpall did not end with NCOL=2");
    events.flush();
    tracks.flush();
    if (!events || !tracks) throw std::runtime_error("failed flushing transport export");
    std::cout << "records=" << summary.records << " histories=" << summary.histories
              << " terminal_ncol=" << summary.terminal_ncol << '\n';
    return 0;
} catch (const std::exception& error) {
    std::cerr << "dumpall_transport_export: " << error.what() << '\n';
    return 1;
}
