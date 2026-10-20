#pragma once

#include <array>
#include <cstddef>
#include <functional>
#include <istream>
#include <vector>

namespace tes_cpp::dump2event {

struct DumpallSecondary {
    int iclust{};
    std::array<int, 9> jclusts{};
    std::array<double, 13> qclusts{};
    std::array<int, 3> jcounts{};
};

struct DumpallRecord {
    int ncol{};
    std::size_t line_start{}, line_end{};
    int particle_id{}, material_id{}, ityp{};
    int ktyp{}, jtyp{}, mtyp{};
    double rtyp{};
    double oldwt{};
    double qs{};
    bool has_qs{};
    int region_before{}, region_after{}, lattice_level_before{}, lattice_level_after{};
    std::vector<std::array<int, 5>> lattice_before, lattice_after;
    double costha{};
    std::array<double, 3> angular_direction{};
    int nsurf{};
    int particle_name{};
    std::array<int, 3> particle_counters{};
    double weight{};
    std::array<double, 3> direction{};
    double energy{}, time{}, x{}, y{}, z{};
    double ec{}, tc{}, xc{}, yc{}, zc{};
    std::array<double, 3> spin{};
    int nzst{};
    int nclsts{};
    int mathz{}, mathn{}, jcoll{}, kcoll{};
    std::vector<DumpallSecondary> secondaries;
    int nocas{}, nobch{};
    double rcasc{}, rsouin{};
};

struct DumpallParseSummary {
    std::size_t records{};
    std::size_t histories{};
    std::array<std::size_t, 16> by_ncol{};
    int terminal_ncol{};
};

// Parses PHITS ASCII Dumpall records using the write order in PHITS 3.341
// src/ovly15.f:dumpdat. Callback is called only for data records (NCOL 4-15).
// Requires a normal NCOL=2 end marker; malformed/truncated records throw with
// their starting line and NCOL context.
DumpallParseSummary parse_dumpall_records(std::istream& input,
    const std::function<void(const DumpallRecord&)>& callback);

}  // namespace tes_cpp::dump2event
