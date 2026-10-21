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
        std::cerr << "usage: extract_selected_dumpall <dumpall> <line-id-list> <csv-out>\n";
        return 2;
    }
    std::ifstream id_file(argv[2]);
    std::set<std::size_t> wanted;
    std::size_t id = 0;
    while (id_file >> id) wanted.insert(id);
    if (!id_file.eof()) throw std::runtime_error("invalid source record ID list");
    std::ifstream in(argv[1], std::ios::binary);
    std::ofstream out(argv[3], std::ios::binary);
    if (!in || !out) throw std::runtime_error("cannot open input/output");
    out << "event_id_NOCAS,source_record_id,ncol,particle_id,material_id,ityp,jcoll,kcoll,E_MeV,EC_MeV,WT,OLDWT,region_before,region_after,x_cm,y_cm,z_cm,xc_cm,yc_cm,zc_cm,child_index,child_ityp,child_state,child_charge,child_kf,total_energy_GeV,rest_mass_GeV,kinetic_energy_MeV,child_weight,child_x_cm,child_y_cm,child_z_cm\n"
        << std::setprecision(17);
    int nocas = -1;
    std::set<std::size_t> seen;
    const auto summary = parse_dumpall_records(in, [&](const DumpallRecord& r) {
        if (r.ncol == 4) nocas = r.nocas;
        if (!wanted.count(r.line_start)) return;
        seen.insert(r.line_start);
        const std::size_t child_count = r.secondaries.empty() ? 1 : r.secondaries.size();
        for (std::size_t c = 0; c < child_count; ++c) {
            out << nocas << ',' << r.line_start << ',' << r.ncol << ',' << r.particle_id << ','
                << r.material_id << ',' << r.ityp << ',' << r.jcoll << ',' << r.kcoll << ','
                << r.energy << ',' << r.ec << ',' << r.weight << ',' << r.oldwt << ','
                << r.region_before << ',' << r.region_after << ',' << r.x << ',' << r.y << ','
                << r.z << ',' << r.xc << ',' << r.yc << ',' << r.zc << ',';
            if (r.secondaries.empty()) {
                out << ",,,,,,,,,,,,\n";
            } else {
                const auto& s = r.secondaries[c];
                out << c << ',' << s.jclusts[3] << ',' << s.jclusts[4] << ',' << s.jclusts[5]
                    << ',' << s.jclusts[7] << ',' << s.qclusts[4] << ',' << s.qclusts[5]
                    << ',' << s.qclusts[7] << ',' << s.qclusts[8] << ',' << s.qclusts[10]
                    << ',' << s.qclusts[11] << ',' << s.qclusts[12] << '\n';
            }
        }
    }, argv[1]);
    out.flush();
    if (!out) throw std::runtime_error("failed writing CSV");
    std::cout << "records=" << summary.records << " histories=" << summary.histories
              << " terminal_ncol=" << summary.terminal_ncol << " selected_source_records="
              << seen.size() << " requested_source_records=" << wanted.size() << '\n';
    if (seen != wanted) return 1;
    return summary.terminal_ncol == 2 ? 0 : 1;
} catch (const std::exception& e) {
    std::cerr << "extract_selected_dumpall: " << e.what() << '\n';
    return 1;
}
