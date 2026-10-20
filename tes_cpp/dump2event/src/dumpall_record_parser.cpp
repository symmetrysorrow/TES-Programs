#include "tes_cpp/dumpall_record_parser.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace tes_cpp::dump2event {
namespace {

std::vector<double> parse_row(const std::string& line, std::size_t line_no,
                              std::size_t start, int ncol, std::size_t expected) {
    std::vector<double> values;
    const char* p = line.c_str();
    while (*p) {
        while (*p == ' ' || *p == '\t' || *p == '\r') ++p;
        if (!*p) break;
        char* end = nullptr;
        const double value = std::strtod(p, &end);
        if (end == p || !std::isfinite(value))
            throw std::runtime_error("Dumpall line " + std::to_string(line_no) +
                " (record at line " + std::to_string(start) + ", NCOL=" +
                std::to_string(ncol) + "): invalid or non-finite numeric token");
        values.push_back(value);
        p = end;
    }
    if (values.size() != expected)
        throw std::runtime_error("Dumpall line " + std::to_string(line_no) +
            " (record at line " + std::to_string(start) + ", NCOL=" +
            std::to_string(ncol) + "): expected " + std::to_string(expected) +
            " columns, found " + std::to_string(values.size()));
    return values;
}

std::vector<double> next_row(std::istream& input, std::string& line, std::size_t& line_no,
                             std::size_t start, int ncol, std::size_t expected) {
    if (!std::getline(input, line))
        throw std::runtime_error("Dumpall EOF inside record at line " +
            std::to_string(start) + " (NCOL=" + std::to_string(ncol) + ")");
    ++line_no;
    return parse_row(line, line_no, start, ncol, expected);
}

int integer(double value, const std::string& field, std::size_t start, int ncol) {
    if (std::trunc(value) != value || value < std::numeric_limits<int>::min() ||
        value > std::numeric_limits<int>::max())
        throw std::runtime_error("Dumpall record at line " + std::to_string(start) +
            " (NCOL=" + std::to_string(ncol) + "): non-integer " + field);
    return static_cast<int>(value);
}

}  // namespace

DumpallParseSummary parse_dumpall_records(std::istream& input,
    const std::function<void(const DumpallRecord&)>& callback) {
    DumpallParseSummary summary;
    std::string line;
    std::size_t line_no = 0;
    bool saw_start = false;
    bool saw_end = false;
    while (std::getline(input, line)) {
        ++line_no;
        const auto marker = parse_row(line, line_no, line_no, 0, 1);
        const int ncol = integer(marker[0], "NCOL", line_no, 0);
        if (ncol == 1) {
            if (saw_start || summary.records != 0)
                throw std::runtime_error("unexpected NCOL=1 at line " + std::to_string(line_no));
            saw_start = true;
            continue;
        }
        if (!saw_start) throw std::runtime_error("Dumpall is missing initial NCOL=1 marker");
        if (ncol == 2) {
            summary.terminal_ncol = 2;
            saw_end = true;
            break;
        }
        if (ncol == 3) {
            summary.terminal_ncol = 3;
            continue;
        }
        if (ncol < 4 || ncol > 15)
            throw std::runtime_error("unsupported NCOL=" + std::to_string(ncol) +
                                     " at Dumpall line " + std::to_string(line_no));

        DumpallRecord record;
        record.ncol = ncol;
        record.line_start = line_no;
        if (ncol == 4) {
            const auto history = next_row(input, line, line_no, record.line_start, ncol, 4);
            record.nocas = integer(history[0], "NOCAS", record.line_start, ncol);
            record.nobch = integer(history[1], "NOBCH", record.line_start, ncol);
            record.rcasc = history[2];
            record.rsouin = history[3];
            ++summary.histories;
        }
        const auto particle = next_row(input, line, line_no, record.line_start, ncol, 8);
        record.particle_id = integer(particle[0], "particle ID", record.line_start, ncol);
        record.material_id = integer(particle[1], "material ID", record.line_start, ncol);
        record.ityp = integer(particle[2], "ITYP", record.line_start, ncol);
        record.ktyp = integer(particle[3], "KTYP", record.line_start, ncol);
        record.jtyp = integer(particle[4], "JTYP", record.line_start, ncol);
        record.mtyp = integer(particle[5], "MTYP", record.line_start, ncol);
        record.rtyp = particle[6];
        record.oldwt = particle[7];
        if (record.ityp == 12 || record.ityp == 13) {
            record.qs = next_row(input, line, line_no, record.line_start, ncol, 1)[0];
            record.has_qs = true;
        }
        const auto regions = next_row(input, line, line_no, record.line_start, ncol, 4);
        record.region_before = integer(regions[0], "region before", record.line_start, ncol);
        record.region_after = integer(regions[1], "region after", record.line_start, ncol);
        record.lattice_level_before = integer(regions[2], "lattice level before", record.line_start, ncol);
        record.lattice_level_after = integer(regions[3], "lattice level after", record.line_start, ncol);
        if (record.lattice_level_before < 0 || record.lattice_level_after < 0)
            throw std::runtime_error("negative lattice level at Dumpall line " +
                                     std::to_string(record.line_start));
        for (int i = 0; i < record.lattice_level_before; ++i) {
            const auto lattice = next_row(input, line, line_no, record.line_start, ncol, 5);
            std::array<int, 5> indices{};
            for (std::size_t k = 0; k < lattice.size(); ++k)
                indices[k] = integer(lattice[k], "lattice-before index", record.line_start, ncol);
            record.lattice_before.push_back(indices);
        }
        for (int i = 0; i < record.lattice_level_after; ++i) {
            const auto lattice = next_row(input, line, line_no, record.line_start, ncol, 5);
            std::array<int, 5> indices{};
            for (std::size_t k = 0; k < lattice.size(); ++k)
                indices[k] = integer(lattice[k], "lattice-after index", record.line_start, ncol);
            record.lattice_after.push_back(indices);
        }
        const auto geometry = next_row(input, line, line_no, record.line_start, ncol, 5);
        record.costha = geometry[0];
        record.angular_direction = {geometry[1], geometry[2], geometry[3]};
        record.nsurf = integer(geometry[4], "NSURF", record.line_start, ncol);
        const auto name_counts = next_row(input, line, line_no, record.line_start, ncol, 4);
        record.particle_name = integer(name_counts[0], "particle name", record.line_start, ncol);
        for (std::size_t k = 0; k < record.particle_counters.size(); ++k)
            record.particle_counters[k] = integer(name_counts[k + 1], "particle counter",
                                                   record.line_start, ncol);
        const auto transport = next_row(input, line, line_no, record.line_start, ncol, 4);
        record.weight = transport[0];
        record.direction = {transport[1], transport[2], transport[3]};
        const auto e = next_row(input, line, line_no, record.line_start, ncol, 5);
        record.energy = e[0]; record.time = e[1]; record.x = e[2]; record.y = e[3]; record.z = e[4];
        const auto ec = next_row(input, line, line_no, record.line_start, ncol, 5);
        record.ec = ec[0]; record.tc = ec[1]; record.xc = ec[2]; record.yc = ec[3]; record.zc = ec[4];
        const auto spin = next_row(input, line, line_no, record.line_start, ncol, 3);
        record.spin = {spin[0], spin[1], spin[2]};
        record.nzst = integer(next_row(input, line, line_no, record.line_start, ncol, 1)[0],
                              "NZST", record.line_start, ncol);
        if (ncol == 13 || ncol == 14) {
            record.nclsts = integer(next_row(input, line, line_no, record.line_start, ncol, 1)[0],
                                    "NCLSTS", record.line_start, ncol);
            if (record.nclsts < 0)
                throw std::runtime_error("negative NCLSTS at Dumpall line " +
                                         std::to_string(record.line_start));
            if (record.nclsts > 0) {
                const auto collision = next_row(input, line, line_no, record.line_start, ncol, 4);
                record.mathz = integer(collision[0], "MATHZ", record.line_start, ncol);
                record.mathn = integer(collision[1], "MATHN", record.line_start, ncol);
                record.jcoll = integer(collision[2], "JCOLL", record.line_start, ncol);
                record.kcoll = integer(collision[3], "KCOLL", record.line_start, ncol);
                record.secondaries.reserve(static_cast<std::size_t>(record.nclsts));
                for (int i = 0; i < record.nclsts; ++i) {
                    DumpallSecondary secondary;
                    secondary.iclust = integer(next_row(input, line, line_no, record.line_start,
                                                         ncol, 1)[0], "ICLUSTS", record.line_start, ncol);
                    const auto j = next_row(input, line, line_no, record.line_start, ncol, 9);
                    for (std::size_t k = 0; k < j.size(); ++k)
                        secondary.jclusts[k] = integer(j[k], "JCLUSTS", record.line_start, ncol);
                    const auto q = next_row(input, line, line_no, record.line_start, ncol, 13);
                    std::copy(q.begin(), q.end(), secondary.qclusts.begin());
                    const auto counts = next_row(input, line, line_no, record.line_start, ncol, 3);
                    for (std::size_t k = 0; k < counts.size(); ++k)
                        secondary.jcounts[k] = integer(counts[k], "JCOUNT", record.line_start, ncol);
                    record.secondaries.push_back(secondary);
                }
            }
        }
        record.line_end = line_no;
        ++summary.records;
        ++summary.by_ncol[static_cast<std::size_t>(ncol)];
        callback(record);
    }
    if (!saw_start) throw std::runtime_error("empty Dumpall or missing initial NCOL=1 marker");
    if (!saw_end) throw std::runtime_error("Dumpall reached EOF without required NCOL=2 end marker");
    return summary;
}

}  // namespace tes_cpp::dump2event
