#include "tes_cpp/physical_deposit.hpp"
#include "tes_cpp/dumpall_record_parser.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace tes_cpp::dump2event;

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

template<class F> void must_throw(F&& call, const char* message) {
    try { call(); } catch (const std::runtime_error&) { return; }
    throw std::runtime_error(message);
}
}

int main() try {
    const std::string record =
        "1\n4\n1 1 1 1\n1 0 14 22 0 0 0 1\n2000 1000 0 0\n"
        "0 0 0 0 0\n1 0 0 0\n1 0 0 0\n1.332 0 -0.8 0 0\n"
        "1.332 0 -0.8 0 0\n0 0 0\n0\n"
        "14\n2 2 12 11 -1 0 -1 1\n0.0005\n2000 1000 0 0\n"
        "0 0 0 0 0\n1 0 0 0\n1 0 0 0\n0.1 0 -0.8 0 0\n"
        "0.1 0 -0.8 0 0\n0 0 0\n0\n0\n"
        "14\n3 2 12 11 -1 0 -1 1\n0.0005\n2000 1000 0 0\n"
        "0 0 0 0 0\n1 0 0 0\n1 0 0 0\n0.1 0 -0.8 0 0\n"
        "0.1 0 -0.8 0 0\n0 0 0\n0\n1\n0 0 14 0\n7\n"
        "0 0 0 12 -1 -1 0 11 0\n0 0 0 0.0005 0.000511 0 0 0.0004 1 0 -0.8 0 0\n"
        "0 0 0\n"
        "15\n4 2 14 22 0 0 0.511 1\n2000 1000 1 1\n"
        "1 2 3 4 5\n6 7 8 9 10\n"
        "0 0 0 0 0\n1 0 0 0\n1 0 0 0\n0.1 0 -0.8 0 0\n"
        "0.1 0 -0.8 0 0\n0 0 0\n0\n3\n2\n";
    std::istringstream parsed_dump(record);
    std::size_t parsed_callbacks = 0;
    const auto parsed = parse_dumpall_records(parsed_dump, [&](const DumpallRecord& item) {
        ++parsed_callbacks;
        if (item.ncol == 4) require(item.nocas == 1 && item.nobch == 1, "NCOL=4 source header");
        if (item.ncol == 14 && item.particle_id == 2)
            require(item.ityp == 12 && item.has_qs && item.nclsts == 0,
                    "NCOL=14 with zero secondaries");
        if (item.ncol == 14 && item.particle_id == 3)
            require(item.nclsts == 1 && item.secondaries.size() == 1 &&
                    item.secondaries[0].jclusts[3] == 12 && item.secondaries[0].jclusts[4] == -1,
                    "NCOL=14 secondary blocks");
        if (item.ncol == 15)
            require(item.particle_id == 4 && !item.has_qs && item.lattice_before.size() == 1 &&
                    item.lattice_after[0][4] == 10 && item.nsurf == 0,
                    "NCOL=15 lattice record");
    });
    require(parsed.records == 4 && parsed_callbacks == 4 && parsed.histories == 1 &&
            parsed.terminal_ncol == 2 && parsed.by_ncol[14] == 2 && parsed.by_ncol[15] == 1,
            "structured Dumpall records and terminators");
    auto wrapped_record = record;
    const std::string unwrapped_particle = "1 1 1 1\n1 0 14 22 0 0 0 1\n";
    const auto wrapped_at = wrapped_record.find(unwrapped_particle);
    require(wrapped_at != std::string::npos, "wrapped-row fixture source found");
    wrapped_record.replace(wrapped_at, unwrapped_particle.size(),
                           "1 1 1 1\n1 0 14 22 0 0\n0 1\n");
    std::istringstream wrapped_dump(wrapped_record);
    const auto wrapped = parse_dumpall_records(wrapped_dump, [](const DumpallRecord&) {});
    require(wrapped.records == parsed.records && wrapped.histories == parsed.histories &&
            wrapped.terminal_ncol == 2,
            "list-directed numeric field rows may wrap across physical lines");
    auto clipped_particle = record;
    clipped_particle.replace(unwrapped_particle.size() ?
        clipped_particle.find(unwrapped_particle) : std::string::npos,
        unwrapped_particle.size(), "1 1 1 1\n1 0 14 22 0 0 0.000\n2000 1000 0 0\n");
    bool clipped_context = false;
    try {
        std::istringstream clipped_dump(clipped_particle);
        (void)parse_dumpall_records(clipped_dump, [](const DumpallRecord&) {},
                                    "gnu-clipped-fixture.dumpall");
    } catch (const std::runtime_error& error) {
        const std::string message = error.what();
        clipped_context = message.find("gnu-clipped-fixture.dumpall") != std::string::npos &&
            message.find("NCOL=4") != std::string::npos &&
            message.find("particle header") != std::string::npos &&
            message.find("not inferred") != std::string::npos;
    }
    require(clipped_context,
            "truncated NCOL=4 particle field fails closed with input and field context");
    must_throw([&] { std::istringstream truncated(record.substr(0, record.size() - 2));
        (void)parse_dumpall_records(truncated, [](const DumpallRecord&) {}); },
        "truncated Dumpall must fail closed");

    const auto path = std::filesystem::temp_directory_path() / "tes_cpp_egs5_test_table.csv";
    {
        std::ofstream table(path);
        table << "medium,medname,lelke,eke0,eke1,ededx0,ededx1,pdedx0,pdedx1,phits_density,egs_density\n";
        table << "1,CU,1,10,1,9,0,0,0,8.9,8.9\n";
        table << "2,PB,1,10,1,9,0,8,0,11.35,11.35\n";
        table << "2,PB,3,10,1,4,0,3,0,11.35,11.35\n";
        table << "2,PB,7,10,1,3,0,2,0,11.35,11.35\n";
        // Synthetic coefficient fixture for Bi code-path tests only.
        // These numbers are NOT HATCH measurements and validate no Bi physics.
        table << "5,BI,1,10,1,9,0,8,0,10,10\\n";
        table << "5,BI,3,10,1,4,0,3,0,10,10\\n";

    }
    const auto table = Egs5StoppingTable::read_csv(path.string(), 2, "Pb");
    require(table.material() == "PB" && table.medium() == 2, "material/medium parse");
    require(std::abs(table.stopping_power(0.1, 2) - 3.0) < 1e-14, "EGS5 table lookup");
    require(std::abs(table.positron_stopping_power(0.1, 2) - 2.0) < 1e-14,
            "EGS5 positron table lookup uses PDedx coefficients");
    require(classify_dead_secondary_score(7.5933446597176157e-7, 1.0e-6) ==
                ScoreGateDecision::not_scored, "below runtime dead-secondary score threshold");
    require(classify_dead_secondary_score(1.0e-6, 1.0e-6) ==
                ScoreGateDecision::not_scored, "score threshold comparison is strict greater-than");
    require(classify_dead_secondary_score(std::nextafter(1.0e-6, 1.0), 1.0e-6) ==
                ScoreGateDecision::scored, "above runtime dead-secondary score threshold");
    require(classify_dead_secondary_score(0.0, std::nullopt) ==
                ScoreGateDecision::unresolved, "unknown score threshold fails closed");
    require(classify_dead_secondary_score(-1.0, 1.0e-6) ==
                ScoreGateDecision::invalid, "negative secondary kinetic energy is invalid");
    PhysicalDepositContext context{"3.341", 1000, "Pb", 2, 0.001, true, 1.0e-6};
    PhysicalDepositLedger ledger(context, table);

    // Bi must require explicit research-only gating, even if a matching
    // synthetic coefficient table is provided. No PHITS oracle has been run.
    const auto bi_table = Egs5StoppingTable::read_csv(path.string(), 5, "Bi");
    must_throw([&] {
        PhysicalDepositLedger no_bi({"3.341", 2007, "Bi", 5, 0.001, true,
                                    1.0e-6, 17, false}, bi_table);
    }, "unvalidated Bi must be rejected in ordinary ledger mode");
    PhysicalDepositLedger bi_research({"3.341", 2007, "Bi", 5, 0.001, true,
                                     1.0e-6, 17, true}, bi_table);
    DepositInput bi_step;
    bi_step.addition_id = "bi-synthetic-1";
    bi_step.event_id = 1; bi_step.particle_id = 1; bi_step.region_id = 2007;
    bi_step.particle_type = 12; bi_step.process = DepositProcess::egs5_electron_continuous;
    bi_step.ncol = 14; bi_step.jcoll = 14; bi_step.energy_mev = 0.1;
    bi_step.segment_length_cm = 0.2;
    const auto bi_result = bi_research.add(bi_step);
    require(bi_result.status == "reconstructed" &&
            std::abs(*bi_result.energy_mev - 0.8) < 1e-14,
            "synthetic Bi coefficients read by medium 5 with independent material ID 17");
    must_throw([&] {
        PhysicalDepositLedger wrong_bi_medium({"3.341", 2007, "Bi", 2, 0.001, true,
                                             1.0e-6, 17, true}, bi_table);
    }, "Bi EGS5 medium must match the supplied HATCH table");


    DepositInput continuous;
    continuous.addition_id = "a1"; continuous.event_id = 3; continuous.particle_id = 2;
    continuous.particle_type = 12; continuous.ncol = 14; continuous.jcoll = 14;
    continuous.process = DepositProcess::egs5_electron_continuous;
    continuous.energy_mev = 0.1; continuous.segment_length_cm = 0.2; continuous.region_id = 1000;
    const auto c = ledger.add(continuous);
    require(c.status == "reconstructed" && std::abs(*c.energy_mev - 0.6) < 1e-14,
            "continuous EGS5 segment loss");
    require(c.position_status == "segment", "transport interval position classification");
    auto repeated_collision = continuous;
    repeated_collision.addition_id = "a1b"; repeated_collision.ncol = 13;
    const auto c13 = ledger.add(repeated_collision);
    require(c13.status == "reconstructed" && std::abs(*c13.energy_mev - 0.6) < 1e-14,
            "NCOL=13 EGS5 electron continuous loss");

    auto positron = continuous;
    positron.addition_id = "p1"; positron.particle_type = 13;
    positron.process = DepositProcess::egs5_positron_continuous;
    positron.energy_mev = positron.energy_c_mev = 0.1;
    positron.old_weight = positron.transport_weight = 1.0;
    const auto p = ledger.add(positron);
    require(p.status == "reconstructed" && std::abs(*p.energy_mev - 0.4) < 1e-14 &&
            p.position_status == "segment", "positron continuous loss uses PDedx and raw segment");

    auto boundary = continuous;
    boundary.addition_id = "a1c"; boundary.ncol = 10; boundary.jcoll = 0;
    boundary.process = DepositProcess::egs5_electron_boundary;
    const auto crossing = ledger.add(boundary);
    require(crossing.status == "reconstructed" && std::abs(*crossing.energy_mev - 0.6) < 1e-14 &&
            crossing.position_status == "segment_to_boundary",
            "NCOL=10 EGS5 loss must use the internal ein-ecc segment formula");

    auto cutoff = continuous;
    cutoff.addition_id = "a2"; cutoff.ncol = 11; cutoff.jcoll = 0;
    cutoff.process = DepositProcess::egs5_electron_cutoff; cutoff.energy_mev = 0.001;
    cutoff.segment_length_cm = 0.2;
    const auto stop = ledger.add(cutoff);
    require(stop.status == "reconstructed" && std::abs(*stop.energy_mev - 0.801) < 1e-14,
            "validated cutoff formula");

    auto dead = continuous;
    dead.addition_id = "a3"; dead.process = DepositProcess::dead_electron_secondary;
    dead.secondary_particle_type = 12; dead.secondary_status = -1;
    dead.secondary_energy_mev = 0.04; dead.secondary_weight = 1.0;
    dead.secondary_xyz_cm = {1.0, 2.0, 3.0};
    const auto d = ledger.add(dead);
    require(d.status == "reconstructed" && *d.energy_mev == 0.04 &&
            d.position_status == "secondary_point", "dead secondary deposit");
    auto below_transport_cutoff = dead;
    below_transport_cutoff.addition_id = "a3b";
    below_transport_cutoff.secondary_energy_mev = 0.0007;
    const auto d_low = ledger.add(below_transport_cutoff);
    require(d_low.status == "reconstructed" && *d_low.energy_mev == 0.0007,
            "dead secondary above score gate is reconstructed even below transport cutoff");
    auto dead_positron = dead;
    dead_positron.addition_id = "p2"; dead_positron.secondary_particle_type = 13;
    dead_positron.process = DepositProcess::dead_positron_secondary;
    dead_positron.secondary_energy_mev = 0.001;
    const auto dp = ledger.add(dead_positron);
    require(dp.status == "reconstructed" && *dp.energy_mev == 0.001 &&
            dp.position_status == "secondary_point",
            "dead secondary positron local kinetic residual, without rest mass");
    auto below_score_guard = dead;
    below_score_guard.addition_id = "a3c";
    below_score_guard.secondary_energy_mev = 5.0e-7;
    require(ledger.add(below_score_guard).status == "unresolved" &&
            classify_dead_secondary_score(5.0e-7, context.dead_secondary_score_threshold_mev) ==
                ScoreGateDecision::not_scored,
            "not-scored candidate is distinct from an actual deposit row");

    auto zero = continuous;
    zero.addition_id = "a4"; zero.process = DepositProcess::exact_zero;
    zero.energy_mev = zero.energy_c_mev = 0.2;
    const auto z = ledger.add(zero);
    require(z.status == "reconstructed" && *z.energy_mev == 0.0 &&
            z.position_status == "none", "strict exact-zero branch");
    zero.addition_id = "a5"; zero.energy_c_mev = std::nextafter(0.2, 1.0);
    require(ledger.add(zero).status == "unresolved", "near-zero must remain unresolved");

    auto unsupported = continuous;
    unsupported.addition_id = "a6"; unsupported.process = DepositProcess::unsupported;
    require(ledger.add(unsupported).status == "unresolved", "unsupported process must fail closed");
    must_throw([&] { ledger.add(continuous); }, "duplicate addition ID must fail");
    must_throw([&] { PhysicalDepositLedger bad({"3.342",1000,"Pb",2,0.001,true,1.0e-6},table); },
               "wrong PHITS version must fail");
    must_throw([&] { PhysicalDepositLedger bad({"3.341",2000,"Pb",2,0.001,true,1.0e-6},table); },
               "wrong region must fail");
    std::filesystem::remove(path);
    std::cout << "physical deposit engine checks passed\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << "physical deposit engine test failure: " << error.what() << '\n';
    return 1;
}
