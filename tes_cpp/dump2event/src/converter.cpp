#include "tes_cpp/converter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <H5Cpp.h>

namespace tes_cpp::dump2event {
namespace {

std::vector<double> parse_numbers(const std::string& line, std::size_t line_number) {
    std::vector<double> values;
    const char* cursor = line.c_str();
    while (*cursor != '\0') {
        while (*cursor == ' ' || *cursor == '\t' || *cursor == '\r') ++cursor;
        if (*cursor == '\0') break;
        char* end = nullptr;
        const double value = std::strtod(cursor, &end);
        if (end == cursor) {
            throw std::runtime_error("invalid numeric token on line " + std::to_string(line_number));
        }
        values.push_back(value);
        cursor = end;
    }
    return values;
}

void require_columns(const std::vector<double>& values, std::size_t count, std::size_t line_number) {
    if (values.size() < count) {
        throw std::runtime_error("expected at least " + std::to_string(count) +
                                 " columns on line " + std::to_string(line_number));
    }
}

void write_array(std::ostream& out, const std::vector<double>& values) {
    out << '[';
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i != 0) out << ", ";
        out << values[i];
    }
    out << ']';
}

void write_event(std::ostream& out, const EventInfo& event, int indent) {
    const std::string pad(indent, ' ');
    const std::string child(indent + 2, ' ');
    out << "{\n";
    out << child << "\"ityp\": " << event.ityp << ",\n";
    const auto field = [&](std::string_view name, const std::vector<double>& values, bool last) {
        out << child << '"' << name << "\": ";
        write_array(out, values);
        out << (last ? "\n" : ",\n");
    };
    field("x", event.x, false);
    field("y", event.y, false);
    field("z", event.z, false);
    field("E", event.energy, false);
    field("x_deposit", event.x_deposit, false);
    field("y_deposit", event.y_deposit, false);
    field("z_deposit", event.z_deposit, false);
    field("E_deposit", event.energy_deposit, true);
    out << pad << '}';
}

}  // namespace

Result read_dump(const std::string& dump_path, const Options& options) {
    std::ifstream file(dump_path, std::ios::binary);
    if (!file) throw std::runtime_error("cannot open dump file: " + dump_path);

    Result result;
    History history;
    double ncol = 1;
    double count = 0;
    double case_number = 0;
    double particle_number = 0;
    double particle_type = 0;
    double deposited_before_collision = 0;
    double energy_at_collision = 0;
    double particle_energy = 0;
    double secondary_energy_sum = 0;
    double secondary_count = 0;
    double secondary_index = 0;
    double secondary_status = 0;
    double collision_type = 0;
    std::vector<double> collision_position(3, 0.0);
    HistorySummary summary;
    bool summary_active = false;
    bool saw_secondary_particle = false;
    bool retain_history = true;
    std::size_t line_number = 0;
    const bool discard_non_full_early = options.full_energy_only || !options.save_all;

    const auto commit_history = [&] {
        if (!summary_active) return;

        // NCOL=13/14 secondary bookkeeping can contain duplicate/dead EGS5
        // entries.  Do not expose the resulting negative local values. Keep
        // their spatial locations, then close the history energy budget using
        // the source energy minus the energy carried out by NCOL=12 leaks.
        const double target_deposit = std::max(
            0.0, options.input_energy - summary.leaked_energy);
        double positive_sum = 0.0;
        for (auto& [_, event] : history) {
            for (double& value : event.energy_deposit) {
                if (value < 0.0) value = 0.0;
                positive_sum += value;
            }
        }
        if (positive_sum > 0.0) {
            const double scale = target_deposit / positive_sum;
            for (auto& [_, event] : history)
                for (double& value : event.energy_deposit) value *= scale;
        }

        double total_deposit = 0.0;
        for (const auto& [_, event] : history) {
            total_deposit += std::accumulate(
                event.energy_deposit.begin(), event.energy_deposit.end(), 0.0);
        }
        const int event_id = static_cast<int>(case_number);
        summary.event_id = event_id;
        summary.total_deposit = total_deposit;
        summary.fully_contained = summary.leaked_particles == 0;
        summary.energy_consistent = std::abs(total_deposit - options.input_energy) <= options.energy_tolerance;
        // A transport history is fully absorbed when no particle reaches a
        // leakage termination. The deposit sum is retained as a diagnostic,
        // but it is not reliable enough for this decision because PHITS may
        // stop low-energy EGS5 particles without writing every residual keV
        // as an explicit deposit record.
        summary.full_energy = summary.fully_contained;
        summary.had_secondary_particle = saw_secondary_particle;
        summary.stored_event = options.save_all || summary.full_energy;
        result.history_summaries.push_back(summary);
        if (summary.stored_event) result.batch.emplace(event_id, history);
        if ((options.save_all || options.collect_full_energy_ids) && summary.full_energy)
            result.full_energy_event_ids.push_back(event_id);
        summary_active = false;
    };

    std::string line;
    while (std::getline(file, line)) {
        ++line_number;
        const auto column = parse_numbers(line, line_number);
        if (column.empty()) continue;

        if (static_cast<int>(ncol) == 1) {
            ncol = 4;
            count = 0;
            case_number = 0;
            particle_number = 0;
            continue;
        }

        const int ncol_int = static_cast<int>(ncol);
        if (ncol_int != 1 && ncol_int != 2 && ncol_int != 3 && ncol_int != 17) {
            if (static_cast<int>(count) == 0) {
                ncol = column[0];
                if (static_cast<int>(ncol) != 4) {
                    ++count;
                } else {
                    commit_history();
                    if (options.max_histories > 0 &&
                        static_cast<int>(result.history_summaries.size()) >= options.max_histories) {
                        break;
                    }
                    history.clear();
                    history[1].ityp = 14;
                    summary = {};
                    summary_active = true;
                    saw_secondary_particle = false;
                    retain_history = true;
                }
            }

            if (static_cast<int>(count) == 1 && static_cast<int>(ncol) == 4) {
                case_number = column[0];
                summary.event_id = static_cast<int>(case_number);
            }
            if (static_cast<int>(count) == 2) {
                require_columns(column, 3, line_number);
                particle_number = column[0];
                particle_type = column[2];
                if (particle_number > 1) saw_secondary_particle = true;
                if (static_cast<int>(particle_type) != 12 && static_cast<int>(particle_type) != 13) ++count;
            }
            if (static_cast<int>(count) == 11) {
                require_columns(column, 2, line_number);
                deposited_before_collision = column[0];
            }
            if (static_cast<int>(count) == 13) {
                require_columns(column, 3, line_number);
                particle_energy = column[0];
                energy_at_collision = column[0];
                collision_position[0] = column[2];
            }
            if (static_cast<int>(count) == 14) {
                require_columns(column, 2, line_number);
                collision_position[1] = column[0];
                collision_position[2] = column[1];
            }
            if (static_cast<int>(count) == 16) {
                if (static_cast<int>(ncol) != 13 && static_cast<int>(ncol) != 14) count = -1;
                const int type = static_cast<int>(particle_type);
                ++summary.particle_records;
                if (static_cast<int>(ncol) == 11) {
                    ++summary.energy_cutoff_particles;
                    if (type == 14) ++summary.energy_cutoff_photons;
                    if (type == 12 || type == 13) ++summary.energy_cutoff_electrons;
                }
                if (static_cast<int>(ncol) == 12) {
                    ++summary.leaked_particles;
                    summary.leaked_energy += std::max(0.0, energy_at_collision);
                    if (type == 14) {
                        ++summary.leaked_photons;
                        if (static_cast<int>(particle_number) == 1) ++summary.leaked_primary_photons;
                        else ++summary.leaked_secondary_photons;
                    }
                    if (discard_non_full_early) {
                        history.clear();
                        retain_history = false;
                    }
                }
                if (static_cast<int>(ncol) == 13 || static_cast<int>(ncol) == 14) {
                    ++summary.reaction_records;
                }
                if (retain_history && (type == 12 || type == 13 || type == 14)) {
                    auto& event = history[static_cast<int>(particle_number)];
                    event.ityp = type;
                    event.x.push_back(collision_position[0]);
                    event.y.push_back(collision_position[1]);
                    event.z.push_back(collision_position[2]);
                    event.energy.push_back(particle_energy);
                    if (static_cast<int>(ncol) == 11) {
                        const double local_deposit = std::max(0.0, energy_at_collision);
                        event.energy_deposit.push_back(local_deposit);
                        event.x_deposit.push_back(collision_position[0]);
                        event.y_deposit.push_back(collision_position[1]);
                        event.z_deposit.push_back(collision_position[2]);
                    }
                }
            }
            if (static_cast<int>(count) == 17) secondary_count = column[0];
            if (static_cast<int>(count) == 18) {
                require_columns(column, 3, line_number);
                collision_type = column[2];
                ncol = 17;
                count = -1;
                secondary_index = 0;
                secondary_energy_sum = 0;
                secondary_status = 0;
                // NCOL=13/14 with no produced particles means that the
                // parent particle is absorbed at this collision.  There is
                // no secondary record from which to trigger the append below,
                // so account for the incoming energy here.
                if (secondary_count <= 0) {
                    const int type = static_cast<int>(particle_type);
                    if (retain_history && (type == 12 || type == 13 || type == 14)) {
                        auto& event = history[static_cast<int>(particle_number)];
                        const double local_deposit = std::max(0.0, deposited_before_collision - energy_at_collision);
                        event.energy_deposit.push_back(local_deposit);
                        event.x_deposit.push_back(collision_position[0]);
                        event.y_deposit.push_back(collision_position[1]);
                        event.z_deposit.push_back(collision_position[2]);
                    }
                    ncol = 13;
                    count = -1;
                    continue;
                }
            }
        }

        if (static_cast<int>(ncol) == 17) {
            if (static_cast<int>(count) == 1) {
                require_columns(column, 5, line_number);
                particle_type = column[3];
                secondary_status = column[4];
            }
            if (static_cast<int>(count) == 5) {
                require_columns(column, 2, line_number);
                particle_energy = column[1];
            }
            if (static_cast<int>(count) == 8) {
                const int type = static_cast<int>(particle_type);
                ++summary.secondary_particles;
                if (type == 14) ++summary.secondary_photons;
                if (type == 12 || type == 13) ++summary.secondary_electrons;
                const double threshold = type == 14
                    ? options.photon_min_energy
                    : (type == 12 || type == 13 ? options.electron_min_energy : 0.0);
                // JCLUSTS(4) is the transport status: 0 is a real particle,
                // while a negative value is dead.  Dead particles do not
                // carry energy out of this collision; their energy belongs
                // in the local deposit.  The same applies to particles below
                // the configured transport cut-off.
                if (secondary_status >= 0 && particle_energy >= threshold)
                    secondary_energy_sum += particle_energy;

                if (static_cast<int>(secondary_index) == static_cast<int>(secondary_count) - 1) {
                    ncol = 13;
                    const int parent_type = static_cast<int>(particle_type);
                    if (retain_history && (parent_type == 12 || parent_type == 13 || parent_type == 14)) {
                        auto& event = history[static_cast<int>(particle_number)];
                        // E is the energy at the preceding event point and
                        // EC is the energy at this event point.  For a
                        // reaction, the generated particles are a transport
                        // bookkeeping record; the local deposition on the
                        // parent track is the non-negative loss E-EC.
                        const double local_deposit = std::max(
                            0.0, deposited_before_collision - energy_at_collision);
                        event.energy_deposit.push_back(local_deposit);
                        event.x_deposit.push_back(collision_position[0]);
                        event.y_deposit.push_back(collision_position[1]);
                        event.z_deposit.push_back(collision_position[2]);
                    }
                } else {
                    ++secondary_index;
                }
                count = -1;
            }
        }
        if (static_cast<int>(ncol) == 3) { ncol = column[0]; count = 0; }
        if (static_cast<int>(ncol) == 2) break;
        ++count;
    }

    // The legacy implementation only committed when it encountered the next
    // event. Committing here preserves the final event as well.
    commit_history();
    return result;
}

void write_event_json(const Batch& batch, const std::string& output_path) {
    std::ofstream out(output_path);
    if (!out) throw std::runtime_error("cannot open JSON output: " + output_path);
    out << std::setprecision(17) << "{\n";
    for (auto outer = batch.begin(); outer != batch.end(); ++outer) {
        out << "  \"" << outer->first << "\": {\n";
        for (auto inner = outer->second.begin(); inner != outer->second.end(); ++inner) {
            out << "    \"" << inner->first << "\": ";
            write_event(out, inner->second, 4);
            out << (std::next(inner) == outer->second.end() ? "\n" : ",\n");
        }
        out << "  }" << (std::next(outer) == batch.end() ? "\n" : ",\n");
    }
    out << "}\n";
}

void write_event_hdf5(const Batch& batch, const std::string& output_path) {
    H5::H5File file(output_path, H5F_ACC_TRUNC);
    file.createAttribute("format", H5::StrType(H5::PredType::C_S1, 15), H5::DataSpace())
        .write(H5::StrType(H5::PredType::C_S1, 15), "tes-dump2event");
    // HDF5 does not create intermediate groups when given a nested path.
    // Create the root container before adding one group per event.
    file.createGroup("/events");
    for (const auto& [event_id, history] : batch) {
        H5::Group event = file.createGroup("/events/" + std::to_string(event_id));
        for (const auto& [particle_id, value] : history) {
            H5::Group particle = event.createGroup(std::to_string(particle_id));
            particle.createAttribute("ityp", H5::PredType::NATIVE_INT, H5::DataSpace())
                .write(H5::PredType::NATIVE_INT, &value.ityp);
            const auto write = [&](const char* name, const std::vector<double>& data) {
                hsize_t size = data.size(); H5::DataSpace space(1, &size);
                H5::DataSet dataset = particle.createDataSet(name, H5::PredType::NATIVE_DOUBLE, space);
                if (!data.empty()) dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
            };
            write("x", value.x); write("y", value.y); write("z", value.z); write("E", value.energy);
            write("x_deposit", value.x_deposit); write("y_deposit", value.y_deposit);
            write("z_deposit", value.z_deposit); write("E_deposit", value.energy_deposit);
        }
    }
}

void write_event_ids(const std::vector<int>& event_ids, const std::string& output_path) {
    std::ofstream out(output_path);
    if (!out) throw std::runtime_error("cannot open event ID output: " + output_path);
    for (const int event_id : event_ids) out << event_id << '\n';
}

void write_history_summary(const std::vector<HistorySummary>& summaries, const std::string& output_path) {
    std::ofstream out(output_path);
    if (!out) throw std::runtime_error("cannot open history summary output: " + output_path);
    out << "event_id,particle_records,had_secondary_particle,secondary_particles,secondary_photons,"
           "secondary_electrons,energy_cutoff_particles,energy_cutoff_photons,energy_cutoff_electrons,"
           "leaked_particles,"
           "leaked_photons,leaked_primary_photons,leaked_secondary_photons,"
           "reaction_records,leaked_energy,total_deposit,"
           "energy_consistent,full_energy,fully_contained,stored_event\n";
    out << std::setprecision(17);
    for (const auto& summary : summaries) {
        out << summary.event_id << ','
            << summary.particle_records << ','
            << (summary.had_secondary_particle ? 1 : 0) << ','
            << summary.secondary_particles << ','
            << summary.secondary_photons << ','
            << summary.secondary_electrons << ','
            << summary.energy_cutoff_particles << ','
            << summary.energy_cutoff_photons << ','
            << summary.energy_cutoff_electrons << ','
            << summary.leaked_particles << ','
            << summary.leaked_photons << ','
            << summary.leaked_primary_photons << ','
            << summary.leaked_secondary_photons << ','
            << summary.reaction_records << ','
            << summary.leaked_energy << ','
            << summary.total_deposit << ','
            << (summary.energy_consistent ? 1 : 0) << ','
            << (summary.full_energy ? 1 : 0) << ','
            << (summary.fully_contained ? 1 : 0) << ','
            << (summary.stored_event ? 1 : 0) << '\n';
    }
}

}  // namespace tes_cpp::dump2event
