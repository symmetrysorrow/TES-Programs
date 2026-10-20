#pragma once

#include <map>
#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace tes_cpp::dump2event {

// PHITS 3.341 constructs eletb(1)=min(esmin,1e-6 MeV).  Values at or below
// this conservative upper bound are left unresolved because esmin is not
// serialized in Dumpall and the exact score/no-score decision is unavailable.
inline constexpr double kDeadSecondaryScoreGuardMeV = 1.0e-6;

// Explicitly scoped to the PHITS 3.341 / EGS5 lead validation described in
// PHYSICAL_DEPOSIT_IMPLEMENTATION_PLAN.md. This is an experimental ledger
// engine, not a general PHITS physics implementation.
struct Egs5StoppingRow {
    int medium{};
    std::string material;
    int log_energy_index{};
    double eke0{}, eke1{}, ededx0{}, ededx1{};
    double phits_density{}, egs_density{};
};

class Egs5StoppingTable {
public:
    static Egs5StoppingTable read_csv(const std::string& path, int expected_medium,
                                      const std::string& expected_material);
    double stopping_power(double kinetic_energy_mev, int medium) const;
    const std::string& material() const noexcept { return material_; }
    int medium() const noexcept { return medium_; }
    double phits_density() const noexcept { return phits_density_; }
    double egs_density() const noexcept { return egs_density_; }

private:
    int medium_{};
    std::string material_;
    double phits_density_{}, egs_density_{};
    std::map<int, Egs5StoppingRow> rows_;
};

enum class DepositProcess {
    egs5_electron_continuous,
    egs5_electron_boundary,
    egs5_electron_cutoff,
    dead_electron_secondary,
    exact_zero,
    unsupported
};

struct PhysicalDepositContext {
    std::string phits_version;
    int region_id{};
    std::string region_material;
    int medium{};
    double electron_cutoff_mev{};
    // This validated slice only covers unit particle weights and region 1000.
    bool unit_weights_only{true};
};

struct DepositInput {
    std::string addition_id;
    int event_id{};
    int particle_id{};
    int particle_type{};
    int ncol{};
    int jcoll{};
    DepositProcess process{DepositProcess::unsupported};
    double energy_mev{};       // Dumpall E for the validated record
    double energy_c_mev{};     // Dumpall EC for exact-zero branch
    double transport_weight{1.0};
    double old_weight{1.0};
    double segment_length_cm{}; // Euclidean distance between raw xyz/cxyz
    double secondary_energy_mev{};
    double secondary_weight{1.0};
    int secondary_status{};     // JCLUSTS(4); negative means dead
    int secondary_particle_type{};
    std::array<double, 3> secondary_xyz_cm{};
    int region_id{};
    std::size_t source_line{};
};

struct DepositResult {
    std::string addition_id;
    int event_id{};
    int particle_id{};
    int region_id{};
    std::string status;         // reconstructed or unresolved
    std::string reason;
    std::optional<double> energy_mev;
    std::string position_status; // segment, secondary_point, none, unresolved
    std::size_t source_line{};
};

class PhysicalDepositLedger {
public:
    PhysicalDepositLedger(PhysicalDepositContext context, const Egs5StoppingTable& table);
    DepositResult add(const DepositInput& input);
    const std::vector<DepositResult>& rows() const noexcept { return rows_; }

private:
    PhysicalDepositContext context_;
    const Egs5StoppingTable& table_;
    std::map<std::string, bool> seen_ids_;
    std::vector<DepositResult> rows_;
};

const char* to_string(DepositProcess process) noexcept;

}  // namespace tes_cpp::dump2event
