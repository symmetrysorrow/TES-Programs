#include "tes_cpp/dumpall_record_parser.hpp"
#include "tes_cpp/physical_deposit.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace tes_cpp::dump2event;

int main(int argc, char** argv) try {
    if (argc != 4) {
        std::cerr << "Usage: dump2event_physical_reconstruct <ASCII dumpall> "
                     "<egs5_stopping_power_coefficients.csv> <ledger.csv>\n"
                     "Scope: PHITS 3.341; Pb region 1000; EGS5 medium 2; 1 keV cutoff; unit weights.\n";
        return 2;
    }
    std::ifstream dump(argv[1], std::ios::binary);
    if (!dump) throw std::runtime_error("cannot open Dumpall file");
    std::ofstream out(argv[3]);
    if (!out) throw std::runtime_error("cannot open output ledger");
    const auto table = Egs5StoppingTable::read_csv(argv[2], 2, "PB");
    PhysicalDepositContext context{"3.341", 1000, "PB", 2, 0.001, true};
    PhysicalDepositLedger ledger(context, table);
    out << "addition_id,event_id,particle_id,source_line,region_id,ncol,jcoll,parent_ityp,parent_E_MeV,"
           "parent_EC_MeV,weight,old_weight,process,status,E_deposit_MeV,position_status,"
           "x,y,z,x_end,y_end,z_end,secondary_index,reason\n" << std::setprecision(17);
    std::size_t next_id = 0;
    int event_id = 0;
    std::size_t reconstructed = 0, unresolved_count = 0, unsupported_count = 0;
    const auto summary = parse_dumpall_records(dump, [&](const DumpallRecord& record) {
        if (record.ncol == 4) event_id = record.nocas;
        if (record.ncol < 10 || record.ncol > 15 || record.region_before != context.region_id) return;
        if (record.material_id != context.medium)
            throw std::runtime_error("region 1000 maps to material ID " +
                std::to_string(record.material_id) + ", expected EGS5 medium 2");

        const auto write_result = [&](DepositInput input, const char* process,
                                      double x, double y, double z,
                                      double xe, double ye, double ze, int secondary_index) {
            input.addition_id = std::to_string(++next_id);
            input.event_id = event_id;
            input.region_id = record.region_before;
            input.source_line = record.line_start;
            const auto result = ledger.add(input);
            out << result.addition_id << ',' << result.event_id << ',' << result.particle_id << ','
                << result.source_line << ',' << result.region_id << ',' << input.ncol << ','
                << input.jcoll << ',' << input.particle_type << ',' << input.energy_mev << ','
                << input.energy_c_mev << ',' << input.transport_weight << ',' << input.old_weight << ','
                << process << ',' << result.status << ',';
            if (result.energy_mev) out << *result.energy_mev;
            out << ',' << result.position_status << ',';
            if (result.position_status != "none")
                out << x << ',' << y << ',' << z << ',' << xe << ',' << ye << ',' << ze
                    << ',' << secondary_index;
            else out << ",,,,,," << secondary_index;
            out << ",\""
                << result.reason << "\"\n";
            if (!out) throw std::runtime_error("failed writing physical deposit ledger");
            if (result.status == "reconstructed") ++reconstructed;
            else ++unresolved_count;
        };

        const double chord = std::hypot(std::hypot(record.xc - record.x, record.yc - record.y),
                                        record.zc - record.z);
        bool primary_written = false;
        if ((record.ncol == 13 || record.ncol == 14) && record.ityp == 12 && record.jcoll == 14) {
            DepositInput input;
            input.particle_id = record.particle_id; input.particle_type = record.ityp;
            input.ncol = record.ncol; input.jcoll = record.jcoll;
            input.process = DepositProcess::egs5_electron_continuous;
            input.energy_mev = record.energy; input.energy_c_mev = record.ec;
            input.transport_weight = record.weight; input.old_weight = record.oldwt;
            input.segment_length_cm = chord;
            write_result(input, "egs5_electron_continuous", record.x, record.y, record.z,
                         record.xc, record.yc, record.zc, -1);
            primary_written = true;
        } else if (record.ncol == 10 && record.ityp == 12) {
            DepositInput input;
            input.particle_id = record.particle_id; input.particle_type = record.ityp;
            input.ncol = record.ncol; input.jcoll = record.jcoll;
            input.process = DepositProcess::egs5_electron_boundary;
            input.energy_mev = record.energy; input.energy_c_mev = record.ec;
            input.transport_weight = record.weight; input.old_weight = record.oldwt;
            input.segment_length_cm = chord;
            write_result(input, "egs5_electron_boundary", record.x, record.y, record.z,
                         record.xc, record.yc, record.zc, -1);
            primary_written = true;
        } else if (record.ncol == 11 && record.ityp == 12) {
            DepositInput input;
            input.particle_id = record.particle_id; input.particle_type = record.ityp;
            input.ncol = record.ncol; input.process = DepositProcess::egs5_electron_cutoff;
            input.energy_mev = record.energy; input.transport_weight = record.weight;
            input.old_weight = record.oldwt; input.segment_length_cm = chord;
            write_result(input, "egs5_electron_cutoff", record.x, record.y, record.z,
                         record.xc, record.yc, record.zc, -1);
            primary_written = true;
        }

        if (!primary_written && (record.ncol == 10 || record.ncol == 13 || record.ncol == 14) &&
            record.energy == record.ec) {
            DepositInput input;
            input.particle_id = record.particle_id; input.particle_type = record.ityp;
            input.ncol = record.ncol; input.jcoll = record.jcoll;
            input.process = DepositProcess::exact_zero;
            input.energy_mev = record.energy; input.energy_c_mev = record.ec;
            input.transport_weight = record.weight; input.old_weight = record.oldwt;
            write_result(input, "exact_zero", record.x, record.y, record.z,
                         record.xc, record.yc, record.zc, -1);
        }

        if (!primary_written && (record.ncol == 10 || record.ncol == 11 ||
                                 record.ncol == 13 || record.ncol == 14) &&
            (record.energy != record.ec || record.ncol == 11)) {
            DepositInput input;
            input.particle_id = record.particle_id; input.particle_type = record.ityp;
            input.ncol = record.ncol; input.jcoll = record.jcoll;
            input.transport_weight = record.weight; input.old_weight = record.oldwt;
            write_result(input, "unsupported_primary", record.x, record.y, record.z,
                         record.xc, record.yc, record.zc, -1);
            ++unsupported_count;
        }

        if ((record.ncol == 13 || record.ncol == 14) && record.nclsts > 0) {
            for (std::size_t i = 0; i < record.secondaries.size(); ++i) {
                const auto& child = record.secondaries[i];
                const double child_energy = child.qclusts[7];
                if (child.jclusts[4] >= 0 || child.jclusts[5] == 0 || child_energy <= 0.0) continue;
                if (child_energy <= kDeadSecondaryScoreGuardMeV) {
                    DepositInput input;
                    input.particle_id = record.particle_id; input.particle_type = record.ityp;
                    input.ncol = record.ncol; input.jcoll = record.jcoll;
                    input.transport_weight = record.weight; input.old_weight = record.oldwt;
                    input.secondary_particle_type = child.jclusts[3];
                    input.secondary_status = child.jclusts[4]; input.secondary_energy_mev = child_energy;
                    input.secondary_weight = child.qclusts[8];
                    input.secondary_xyz_cm = {child.qclusts[10], child.qclusts[11], child.qclusts[12]};
                    write_result(input, "unsupported_secondary_threshold", child.qclusts[10],
                                 child.qclusts[11], child.qclusts[12], child.qclusts[10],
                                 child.qclusts[11], child.qclusts[12], static_cast<int>(i));
                    ++unsupported_count;
                    continue;
                }
                if (child.jclusts[3] != 12) {
                    DepositInput input;
                    input.particle_id = record.particle_id; input.particle_type = record.ityp;
                    input.ncol = record.ncol; input.jcoll = record.jcoll;
                    input.transport_weight = record.weight; input.old_weight = record.oldwt;
                    input.secondary_particle_type = child.jclusts[3];
                    input.secondary_status = child.jclusts[4]; input.secondary_energy_mev = child_energy;
                    input.secondary_weight = child.qclusts[8];
                    input.secondary_xyz_cm = {child.qclusts[10], child.qclusts[11], child.qclusts[12]};
                    write_result(input, "unsupported_secondary", child.qclusts[10], child.qclusts[11],
                                 child.qclusts[12], child.qclusts[10], child.qclusts[11],
                                 child.qclusts[12], static_cast<int>(i));
                    ++unsupported_count;
                    continue;
                }
                DepositInput input;
                input.particle_id = record.particle_id; input.particle_type = record.ityp;
                input.ncol = record.ncol; input.jcoll = record.jcoll;
                input.process = DepositProcess::dead_electron_secondary;
                input.transport_weight = record.weight; input.old_weight = record.oldwt;
                input.secondary_particle_type = child.jclusts[3];
                input.secondary_status = child.jclusts[4]; input.secondary_energy_mev = child_energy;
                input.secondary_weight = child.qclusts[8];
                input.secondary_xyz_cm = {child.qclusts[10], child.qclusts[11], child.qclusts[12]};
                write_result(input, "dead_electron_secondary", child.qclusts[10], child.qclusts[11],
                             child.qclusts[12], child.qclusts[10], child.qclusts[11],
                             child.qclusts[12], static_cast<int>(i));
            }
        }

    });
    std::cout << "records=" << summary.records << " histories=" << summary.histories
              << " additions=" << next_id << " reconstructed=" << reconstructed
              << " unresolved=" << unresolved_count << " unsupported=" << unsupported_count
              << " status=EXPERIMENTAL_PARTIAL terminator=NCOL=" << summary.terminal_ncol << '\n';
    return unresolved_count == 0 ? 0 : 1;
} catch (const std::exception& error) {
    std::cerr << "dump2event_physical_reconstruct: " << error.what() << '\n';
    return 1;
}
