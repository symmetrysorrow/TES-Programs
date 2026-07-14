#include "tes_cpp/converter.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string_view>

namespace tes_cpp::dump2json {
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
    constexpr double electron_min_energy = 0.1;
    constexpr double photon_min_energy = 0.001;

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
    double particle_energy = 0;
    double secondary_energy_sum = 0;
    double secondary_count = 0;
    double secondary_index = 0;
    double collision_type = 0;
    std::vector<double> collision_position(3, 0.0);

    const auto commit_history = [&] {
        if (particle_number <= 1) return;
        double total_deposit = 0.0;
        for (const auto& [_, event] : history) {
            total_deposit += std::accumulate(
                event.energy_deposit.begin(), event.energy_deposit.end(), 0.0);
        }
        const bool full_energy = std::abs(total_deposit - options.input_energy) <= options.energy_tolerance;
        const int event_id = static_cast<int>(case_number);
        if (options.save_all || full_energy) result.batch.emplace(event_id, history);
        if (options.save_all && full_energy) result.full_energy_event_ids.push_back(event_id);
    };

    std::string line;
    std::size_t line_number = 0;
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
                    history.clear();
                    history[1].ityp = 14;
                }
            }

            if (static_cast<int>(count) == 1 && static_cast<int>(ncol) == 4) {
                case_number = column[0];
            }
            if (static_cast<int>(count) == 2) {
                require_columns(column, 3, line_number);
                particle_number = column[0];
                particle_type = column[2];
                if (static_cast<int>(particle_type) != 12 && static_cast<int>(particle_type) != 13) ++count;
            }
            if (static_cast<int>(count) == 11) {
                require_columns(column, 2, line_number);
                deposited_before_collision = column[0];
            }
            if (static_cast<int>(count) == 13) {
                require_columns(column, 3, line_number);
                particle_energy = column[0];
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
                if (type == 12 || type == 13 || type == 14) {
                    auto& event = history[static_cast<int>(particle_number)];
                    event.ityp = type;
                    event.x.push_back(collision_position[0]);
                    event.y.push_back(collision_position[1]);
                    event.z.push_back(collision_position[2]);
                    event.energy.push_back(particle_energy);
                    if (static_cast<int>(ncol) == 11) {
                        event.energy_deposit.push_back(deposited_before_collision);
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
                if (static_cast<int>(collision_type) == 14) {
                    ++count;
                    continue;
                }
            }
        }

        if (static_cast<int>(ncol) == 17) {
            if (static_cast<int>(count) == 1) {
                require_columns(column, 4, line_number);
                particle_type = column[3];
            }
            if (static_cast<int>(count) == 5) {
                require_columns(column, 2, line_number);
                particle_energy = column[1];
            }
            if (static_cast<int>(count) == 8) {
                const int type = static_cast<int>(particle_type);
                if (type == 12 || type == 13 || type == 14) {
                    const double threshold = type == 14 ? photon_min_energy : electron_min_energy;
                    if (particle_energy >= threshold) secondary_energy_sum += particle_energy;
                    if (static_cast<int>(secondary_index) == static_cast<int>(secondary_count) - 1) {
                        ncol = 13;
                        auto& event = history[static_cast<int>(particle_number)];
                        event.energy_deposit.push_back(deposited_before_collision - secondary_energy_sum);
                        event.x_deposit.push_back(collision_position[0]);
                        event.y_deposit.push_back(collision_position[1]);
                        event.z_deposit.push_back(collision_position[2]);
                    } else {
                        ++secondary_index;
                    }
                    count = -1;
                }
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

void write_batch_json(const Batch& batch, const std::string& output_path) {
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

void write_event_ids(const std::vector<int>& event_ids, const std::string& output_path) {
    std::ofstream out(output_path);
    if (!out) throw std::runtime_error("cannot open event ID output: " + output_path);
    for (const int event_id : event_ids) out << event_id << '\n';
}

}  // namespace tes_cpp::dump2json
