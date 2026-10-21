#include "tes_cpp/physical_deposit.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <utility>

namespace tes_cpp::dump2event {

ScoreGateDecision classify_dead_secondary_score(double kinetic_energy_mev,
                                                  std::optional<double> eletb1_mev) noexcept {
    if (!std::isfinite(kinetic_energy_mev) || kinetic_energy_mev < 0.0)
        return ScoreGateDecision::invalid;
    if (!eletb1_mev || !std::isfinite(*eletb1_mev) || *eletb1_mev <= 0.0)
        return ScoreGateDecision::unresolved;
    return kinetic_energy_mev > *eletb1_mev ? ScoreGateDecision::scored :
                                              ScoreGateDecision::not_scored;
}
namespace {

std::string upper(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::toupper(ch)); });
    return value;
}

double number(const std::string& value, std::size_t line, const char* column) {
    std::size_t used = 0;
    double parsed = 0.0;
    try { parsed = std::stod(value, &used); }
    catch (const std::exception&) {
        throw std::runtime_error("stopping table line " + std::to_string(line) +
                                 ": invalid " + column);
    }
    if (used != value.size() || !std::isfinite(parsed))
        throw std::runtime_error("stopping table line " + std::to_string(line) +
                                 ": invalid/non-finite " + column);
    return parsed;
}

std::vector<std::string> csv_fields(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) fields.push_back(field);
    if (!line.empty() && line.back() == ',') fields.emplace_back();
    return fields;
}

bool finite_nonnegative(double value) { return std::isfinite(value) && value >= 0.0; }

DepositResult unresolved(const DepositInput& input, std::string why) {
    return {input.addition_id, input.event_id, input.particle_id, input.region_id,
            "unresolved", std::move(why), std::nullopt, "unresolved", input.source_line};
}

}  // namespace

Egs5StoppingTable Egs5StoppingTable::read_csv(const std::string& path, int expected_medium,
                                               const std::string& expected_material) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open EGS5 stopping table: " + path);
    Egs5StoppingTable table;
    table.medium_ = expected_medium;
    table.material_ = upper(expected_material);
    std::string line;
    std::size_t line_number = 0;
    if (!std::getline(input, line)) throw std::runtime_error("empty EGS5 stopping table");
    ++line_number;
    const auto header = csv_fields(line);
    const std::vector<std::string> expected = {"medium", "medname", "lelke", "eke0", "eke1",
        "ededx0", "ededx1", "pdedx0", "pdedx1", "phits_density", "egs_density"};
    if (header != expected) throw std::runtime_error("unexpected EGS5 stopping table columns");
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty()) continue;
        const auto f = csv_fields(line);
        if (f.size() != 11) throw std::runtime_error("stopping table line " +
            std::to_string(line_number) + ": expected 11 columns");
        const int medium = static_cast<int>(number(f[0], line_number, "medium"));
        const auto material = upper(f[1]);
        if (medium != expected_medium) continue;
        if (material != table.material_)
            throw std::runtime_error("stopping table medium maps to an unexpected material");
        const int index = static_cast<int>(number(f[2], line_number, "lelke"));
        Egs5StoppingRow row{medium, material, index,
            number(f[3], line_number, "eke0"), number(f[4], line_number, "eke1"),
            number(f[5], line_number, "ededx0"), number(f[6], line_number, "ededx1"),
            number(f[7], line_number, "pdedx0"), number(f[8], line_number, "pdedx1"),
            number(f[9], line_number, "phits_density"), number(f[10], line_number, "egs_density")};
        if (index <= 0 || row.phits_density <= 0.0 || row.egs_density <= 0.0)
            throw std::runtime_error("stopping table contains invalid index or density");
        if (!table.rows_.emplace(index, row).second)
            throw std::runtime_error("duplicate EGS5 stopping-table energy index " + std::to_string(index));
    }
    if (table.rows_.empty() || table.rows_.find(1) == table.rows_.end())
        throw std::runtime_error("stopping table is missing required index 1");
    table.phits_density_ = table.rows_.begin()->second.phits_density;
    table.egs_density_ = table.rows_.begin()->second.egs_density;
    for (const auto& [index, row] : table.rows_) {
        (void)index;
        if (row.phits_density != table.phits_density_ || row.egs_density != table.egs_density_)
            throw std::runtime_error("inconsistent density metadata in stopping table");
    }
    return table;
}

double Egs5StoppingTable::stopping_power(double kinetic_energy_mev, int medium) const {
    if (medium != medium_) throw std::runtime_error("requested EGS5 medium does not match table");
    if (!std::isfinite(kinetic_energy_mev) || kinetic_energy_mev < 0.0)
        throw std::runtime_error("invalid kinetic energy for stopping-power lookup");
    const double log_energy = std::log(std::max(kinetic_energy_mev, 0.001));
    const auto& base = rows_.at(1);
    const double raw_index = base.eke1 * log_energy + base.eke0;
    if (!std::isfinite(raw_index) || raw_index < 1.0 ||
        raw_index > static_cast<double>(std::numeric_limits<int>::max()))
        throw std::runtime_error("EGS5 stopping-power index is outside table bounds");
    const int index = static_cast<int>(raw_index); // PHITS/EGS5 integer truncation
    const auto found = rows_.find(index);
    if (found == rows_.end()) throw std::runtime_error("EGS5 stopping-power table lacks index " +
                                                       std::to_string(index));
    const double value = found->second.ededx1 * log_energy + found->second.ededx0;
    if (!std::isfinite(value) || value < 0.0)
        throw std::runtime_error("EGS5 stopping-power table produced invalid stopping power");
    return value;
}

double Egs5StoppingTable::positron_stopping_power(double kinetic_energy_mev, int medium) const {
    if (medium != medium_) throw std::runtime_error("requested EGS5 medium does not match table");
    if (!std::isfinite(kinetic_energy_mev) || kinetic_energy_mev < 0.0)
        throw std::runtime_error("invalid kinetic energy for positron stopping-power lookup");
    const double log_energy = std::log(std::max(kinetic_energy_mev, 0.001));
    const auto& base = rows_.at(1);
    const double raw_index = base.eke1 * log_energy + base.eke0;
    if (!std::isfinite(raw_index) || raw_index < 1.0 ||
        raw_index > static_cast<double>(std::numeric_limits<int>::max()))
        throw std::runtime_error("EGS5 positron stopping-power index is outside table bounds");
    const int index = static_cast<int>(raw_index); // PHITS/EGS5 integer truncation
    const auto found = rows_.find(index);
    if (found == rows_.end()) throw std::runtime_error("EGS5 stopping-power table lacks positron index " +
                                                       std::to_string(index));
    const double value = found->second.pdedx1 * log_energy + found->second.pdedx0;
    if (!std::isfinite(value) || value < 0.0)
        throw std::runtime_error("EGS5 positron stopping-power table produced invalid stopping power");
    return value;
}

PhysicalDepositLedger::PhysicalDepositLedger(PhysicalDepositContext context,
                                               const Egs5StoppingTable& table)
    : context_(std::move(context)), table_(table) {
    if (context_.phits_version != "3.341")
        throw std::runtime_error("experimental deposit model supports PHITS 3.341 only");
    const auto material = upper(context_.region_material);
    const bool pb_validated = material == "PB" && context_.region_id == 1000 &&
        context_.medium == 2 && context_.dumpall_material_id == 2;
    const bool bi_research = material == "BI" && context_.bi_research_only &&
        context_.region_id > 0 && context_.medium > 0 && context_.dumpall_material_id > 0;
    if (!pb_validated && !bi_research)
        throw std::runtime_error("only validated Pb settings or explicitly unvalidated Bi research are supported");
    if (table_.medium() != context_.medium || table_.material() != material)
        throw std::runtime_error("EGS5 stopping table medium/material mismatch");
    if (context_.electron_cutoff_mev != 0.001)
        throw std::runtime_error("experimental deposit model requires the validated 1 keV electron cutoff");
    if (!context_.unit_weights_only)
        throw std::runtime_error("non-unit particle weights are not supported by this validated model");
}

DepositResult PhysicalDepositLedger::add(const DepositInput& in) {
    if (in.addition_id.empty()) throw std::runtime_error("blank physical tally addition ID");
    if (!seen_ids_.emplace(in.addition_id, true).second)
        throw std::runtime_error("duplicate physical tally addition ID: " + in.addition_id);
    DepositResult result = unresolved(in, "process or record conditions are outside validated scope");
    if (in.region_id != context_.region_id) {
        result.reason = "record region is outside the configured scoring region";
    } else if (!finite_nonnegative(in.transport_weight) || !finite_nonnegative(in.old_weight) ||
               in.transport_weight != 1.0 || in.old_weight != 1.0) {
        result.reason = "validated reconstruction requires transport and old weights equal to 1";
    } else if (in.process == DepositProcess::egs5_electron_continuous &&
               (in.ncol == 13 || in.ncol == 14) && in.jcoll == 14 && in.particle_type == 12 &&
               finite_nonnegative(in.energy_mev) && finite_nonnegative(in.segment_length_cm)) {
        const double value = table_.stopping_power(in.energy_mev, context_.medium) * in.segment_length_cm;
        result = {in.addition_id, in.event_id, in.particle_id, in.region_id, "reconstructed",
                  "EGS5 stopping power at raw Dumpall E multiplied by raw endpoint chord",
                  value, "segment", in.source_line};
    } else if (in.process == DepositProcess::egs5_positron_continuous &&
               (in.ncol == 13 || in.ncol == 14) && in.particle_type == 13 &&
               finite_nonnegative(in.energy_mev) && finite_nonnegative(in.segment_length_cm)) {
        const double value = in.old_weight *
            table_.positron_stopping_power(in.energy_mev, context_.medium) * in.segment_length_cm;
        result = {in.addition_id, in.event_id, in.particle_id, in.region_id, "reconstructed",
                  "PHITS 3.341 EGS5 positron tlv = OLDWT * PDedx(E) * raw endpoint chord",
                  value, "segment", in.source_line};
    } else if (in.process == DepositProcess::egs5_electron_boundary && in.ncol == 10 &&
               in.particle_type == 12 &&
               finite_nonnegative(in.energy_mev) && finite_nonnegative(in.segment_length_cm)) {
        const double value = table_.stopping_power(in.energy_mev, context_.medium) * in.segment_length_cm;
        result = {in.addition_id, in.event_id, in.particle_id, in.region_id, "reconstructed",
                  "PHITS 3.341 EGS5 boundary tlv = ein - ecc = dedx(E) * segment length",
                  value, "segment_to_boundary", in.source_line};
    } else if (in.process == DepositProcess::egs5_electron_cutoff &&
               in.ncol == 11 && in.particle_type == 12 &&
               finite_nonnegative(in.energy_mev) && finite_nonnegative(in.segment_length_cm)) {
        const double value = in.energy_mev +
            table_.stopping_power(context_.electron_cutoff_mev, context_.medium) * in.segment_length_cm;
        result = {in.addition_id, in.event_id, in.particle_id, in.region_id, "reconstructed",
                  "raw residual E plus validated stopping power at electron cutoff times final chord",
                  value, "segment", in.source_line};
    } else if (in.process == DepositProcess::dead_electron_secondary &&
               (in.ncol == 13 || in.ncol == 14) &&
               in.secondary_particle_type == 12 && in.secondary_status < 0 &&
               classify_dead_secondary_score(in.secondary_energy_mev,
                   context_.dead_secondary_score_threshold_mev) == ScoreGateDecision::scored &&
               in.secondary_weight == 1.0 && std::isfinite(in.secondary_xyz_cm[0]) &&
               std::isfinite(in.secondary_xyz_cm[1]) && std::isfinite(in.secondary_xyz_cm[2])) {
        result = {in.addition_id, in.event_id, in.particle_id, in.region_id, "reconstructed",
                  "dead secondary electron QCLUSTS(7) strictly above resolved PHITS eletb(1), unit QCLUSTS(8)",
                  in.secondary_energy_mev, "secondary_point", in.source_line};
    } else if (in.process == DepositProcess::dead_positron_secondary &&
               (in.ncol == 13 || in.ncol == 14) &&
               in.secondary_particle_type == 13 && in.secondary_status < 0 &&
               classify_dead_secondary_score(in.secondary_energy_mev,
                   context_.dead_secondary_score_threshold_mev) == ScoreGateDecision::scored &&
               in.secondary_weight == 1.0 && std::isfinite(in.secondary_xyz_cm[0]) &&
               std::isfinite(in.secondary_xyz_cm[1]) && std::isfinite(in.secondary_xyz_cm[2])) {
        result = {in.addition_id, in.event_id, in.particle_id, in.region_id, "reconstructed",
                  "dead secondary positron QCLUSTS(7) strictly above resolved PHITS eletb(1), unit weight",
                  in.secondary_energy_mev, "secondary_point", in.source_line};
    } else if (in.process == DepositProcess::exact_zero &&
               (in.ncol == 10 || in.ncol == 13 || in.ncol == 14) &&
               finite_nonnegative(in.energy_mev) && finite_nonnegative(in.energy_c_mev) &&
               in.energy_mev == in.energy_c_mev) {
        result = {in.addition_id, in.event_id, in.particle_id, in.region_id, "reconstructed",
                  "raw Dumpall E equals EC exactly; zero contribution has no position", 0.0,
                  "none", in.source_line};
    }
    if (result.energy_mev && (!std::isfinite(*result.energy_mev) || *result.energy_mev < 0.0))
        throw std::runtime_error("physical reconstruction produced invalid negative/non-finite energy");
    rows_.push_back(result);
    return result;
}

const char* to_string(DepositProcess process) noexcept {
    switch (process) {
    case DepositProcess::egs5_electron_continuous: return "egs5_electron_continuous";
    case DepositProcess::egs5_electron_boundary: return "egs5_electron_boundary";
    case DepositProcess::egs5_electron_cutoff: return "egs5_electron_cutoff";
    case DepositProcess::egs5_positron_continuous: return "egs5_positron_continuous";
    case DepositProcess::dead_electron_secondary: return "dead_electron_secondary";
    case DepositProcess::dead_positron_secondary: return "dead_positron_secondary";
    case DepositProcess::exact_zero: return "exact_zero";
    default: return "unsupported";
    }
}

}  // namespace tes_cpp::dump2event
