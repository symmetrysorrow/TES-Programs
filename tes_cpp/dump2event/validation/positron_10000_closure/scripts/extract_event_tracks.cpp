#include "tes_cpp/dumpall_record_parser.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>

using namespace tes_cpp::dump2event;

int main(int argc, char** argv) try {
    if (argc != 4) {
        std::cerr << "usage: extract_event_tracks <dumpall> <event-id-list> <csv-out>\n";
        return 2;
    }
    std::ifstream ids(argv[2]);
    std::set<int> wanted;
    int id = 0;
    while (ids >> id) wanted.insert(id);
    if (wanted.empty() || !ids.eof()) throw std::runtime_error("invalid/empty event ID list");
    std::ifstream input(argv[1], std::ios::binary);
    std::ofstream output(argv[3], std::ios::binary);
    if (!input || !output) throw std::runtime_error("cannot open input/output");
    output << "NOCAS,source_record_id,NCOL,particle_id,material_id,ITYP,JCOLL,KCOLL,E_MeV,EC_MeV,WT,OLDWT,has_QS,QS_MeV_per_cm,region_before,region_after,x_cm,y_cm,z_cm,xc_cm,yc_cm,zc_cm,secondary_index,child_ITYP,child_state,child_charge,child_KF,child_total_energy_GeV,child_rest_mass_GeV,child_kinetic_energy_MeV,child_weight,child_x_cm,child_y_cm,child_z_cm\n"
           << std::setprecision(17);
    int nocas = -1;
    std::set<int> seen;
    const auto result = parse_dumpall_records(input, [&](const DumpallRecord& r) {
        if (r.ncol == 4) nocas = r.nocas;
        if (!wanted.count(nocas)) return;
        seen.insert(nocas);
        const auto write_parent = [&]() {
            output << nocas << ',' << r.line_start << ',' << r.ncol << ',' << r.particle_id << ','
                   << r.material_id << ',' << r.ityp << ',' << r.jcoll << ',' << r.kcoll << ','
                   << r.energy << ',' << r.ec << ',' << r.weight << ',' << r.oldwt << ',' << (r.has_qs ? 1 : 0) << ',' << r.qs << ','
                   << r.region_before << ',' << r.region_after << ',' << r.x << ',' << r.y << ','
                   << r.z << ',' << r.xc << ',' << r.yc << ',' << r.zc << ',';
        };
        if (r.secondaries.empty()) {
            write_parent();
            output << ",,,,,,,,,,,\n";
        } else {
            for (std::size_t i = 0; i < r.secondaries.size(); ++i) {
                const auto& s = r.secondaries[i];
                write_parent();
                output << i << ',' << s.jclusts[3] << ',' << s.jclusts[4] << ',' << s.jclusts[5]
                       << ',' << s.jclusts[7] << ',' << s.qclusts[4] << ',' << s.qclusts[5]
                       << ',' << s.qclusts[7] << ',' << s.qclusts[8] << ',' << s.qclusts[10]
                       << ',' << s.qclusts[11] << ',' << s.qclusts[12] << '\n';
            }
        }
    }, argv[1]);
    output.flush();
    if (!output) throw std::runtime_error("failed writing CSV");
    std::cout << "records=" << result.records << " histories=" << result.histories
              << " terminal_ncol=" << result.terminal_ncol << " selected_histories=" << seen.size()
              << " requested_histories=" << wanted.size() << '\n';
    if (seen != wanted || result.terminal_ncol != 2) return 1;
    return 0;
} catch (const std::exception& e) {
    std::cerr << "extract_event_tracks: " << e.what() << '\n';
    return 1;
}
