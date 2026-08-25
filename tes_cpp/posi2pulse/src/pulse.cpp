#include "tes_cpp/pulse.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <atomic>
#include <chrono>
#include <complex>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <utility>
#include <nlohmann/json.hpp>

namespace tes_cpp {
namespace {
struct Input {
    double c_abs, c_tes, g_abs_abs, g_abs_tes, g_tes_bath;
    double resistance, load_resistance, t_c, t_bath, alpha, beta, inductance, exponent;
    double energy, rate, samples;
    int n_abs;
};

Input read_input(const std::string& path) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("cannot open input JSON: " + path);
    nlohmann::json j;
    file >> j;
    return {
        j.at("C_abs").get<double>(), j.at("C_tes").get<double>(),
        j.at("G_abs-abs").get<double>(), j.at("G_abs-tes").get<double>(), j.at("G_tes-bath").get<double>(),
        j.at("R").get<double>(), j.at("R_l").get<double>(), j.at("T_c").get<double>(), j.at("T_bath").get<double>(),
        j.at("alpha").get<double>(), j.at("beta").get<double>(), j.at("L").get<double>(), j.at("n").get<double>(),
        j.at("E").get<double>(), j.at("rate").get<double>(), j.at("samples").get<double>(), j.at("n_abs").get<int>()
    };
}

std::vector<double> linspace(double start, double stop, int count) {
    if (count < 1) return {};
    if (count == 1) return {start};
    std::vector<double> values(count);
    const double step = (stop - start) / static_cast<double>(count - 1);
    for (int i = 0; i < count; ++i) values[i] = start + i * step;
    // Match numpy.linspace's inclusive endpoint exactly.
    values.back() = stop;
    return values;
}

Eigen::MatrixXd make_matrix(const Input& in) {
    const int n = in.n_abs;
    const int size = n + 4;
    const double c_abs = in.c_abs / n;
    const double g_abs_abs = in.g_abs_abs * (n - 1);
    const double current = std::sqrt((in.g_tes_bath * in.t_c * (1 - std::pow(in.t_bath / in.t_c, in.exponent))) /
                                     (in.exponent * in.resistance));
    const double t_el = in.inductance / (in.load_resistance + in.resistance * (1 + in.beta));
    const double loop_gain = (in.alpha * current * current * in.resistance) / (in.g_tes_bath * in.t_c);
    const double t_i = in.c_tes / ((1 - loop_gain) * in.g_tes_bath);

    Eigen::MatrixXd a = Eigen::MatrixXd::Zero(size, size);
    for (int i = 0; i < size; ++i) {
        if (i > 0) a(i, i - 1) = -g_abs_abs / c_abs;
        a(i, i) = 2 * g_abs_abs / c_abs;
        if (i + 1 < size) a(i, i + 1) = -g_abs_abs / c_abs;
    }
    a(0, 0) = 1 / t_el;
    a(0, 1) = loop_gain * in.g_tes_bath / (current * in.inductance);
    a(1, 0) = -current * in.resistance * (2 + in.beta) / in.c_tes;
    a(1, 1) = 1 / t_i + in.g_abs_tes / in.c_tes;
    a(1, 2) = -in.g_abs_tes / in.c_tes;
    a(2, 1) = -in.g_abs_tes / c_abs;
    a(2, 2) = in.g_abs_tes / c_abs + g_abs_abs / c_abs;
    a(2, 3) = -g_abs_abs / c_abs;
    const int n1 = n + 1, n2 = n + 2, n3 = n + 3;
    a(n1, n) = -g_abs_abs / c_abs;
    a(n1, n1) = in.g_abs_tes / c_abs + g_abs_abs / c_abs;
    a(n1, n2) = -in.g_abs_tes / c_abs;
    a(n2, n1) = -in.g_abs_tes / in.c_tes;
    a(n2, n2) = 1 / t_i + in.g_abs_tes / in.c_tes;
    a(n2, n3) = -current * in.resistance * (2 + in.beta) / in.c_tes;
    a(n3, n2) = loop_gain * in.g_tes_bath / (current * in.inductance);
    a(n3, n3) = 1 / t_el;
    return -a;
}

void write_array(std::ostream& out, const std::vector<double>& values) {
    out << '[';
    for (std::size_t i = 0; i < values.size(); ++i) { if (i) out << ','; out << values[i]; }
    out << ']';
}

void write_pulse_channels_json(std::ostream& out, const Pulse& pulse) {
    out << "{\"ch0\":";
    write_array(out, pulse.ch0);
    out << ",\"ch1\":";
    write_array(out, pulse.ch1);
    out << '}';
}

nlohmann::json read_input_document(const std::string& path) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("cannot open input JSON: " + path);
    nlohmann::json document;
    file >> document;
    return document;
}

Pulse make_pulse(
    const Input& in,
    const Eigen::MatrixXcd& eigenvectors,
    const Eigen::VectorXcd& eigenvalues,
    const std::vector<double>& time,
    int position,
    const Eigen::VectorXcd& constants) {
    const int samples = static_cast<int>(time.size());
    Pulse pulse;
    pulse.position = position;
    pulse.time = time;
    pulse.ch0.assign(samples, 0.0);
    pulse.ch1.assign(samples, 0.0);
    for (int mode = 0; mode < constants.size(); ++mode) {
        const std::complex<double> eigenvalue = eigenvalues[mode];
        for (int sample = 0; sample < samples; ++sample) {
            // ``pulse_model.py`` defines its output current with the
            // opposite sign from the state-vector reconstruction.
            const std::complex<double> factor =
                -constants[mode] * std::exp(eigenvalue * time[sample]);
            pulse.ch0[sample] += (factor * eigenvectors(0, mode)).real();
            pulse.ch1[sample] += (factor * eigenvectors(in.n_abs + 3, mode)).real();
        }
    }
    return pulse;
}

struct TemporaryDirectory {
    std::filesystem::path path;
    ~TemporaryDirectory() {
        std::error_code error;
        std::filesystem::remove_all(path, error);
    }
};

TemporaryDirectory make_temporary_directory(const std::string& output_path) {
    const auto base = std::filesystem::path(output_path).filename().string() + ".posi2pulse-";
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    for (int suffix = 0; suffix != 100; ++suffix) {
        const auto candidate = std::filesystem::temp_directory_path() /
            (base + std::to_string(stamp) + "-" + std::to_string(suffix));
        std::error_code error;
        if (std::filesystem::create_directory(candidate, error)) return {candidate};
    }
    throw std::runtime_error("cannot create temporary directory for pulse output");
}
}  // namespace

std::vector<Pulse> generate_pulses(const std::string& input_json_path, const std::vector<int>& positions) {
    const Input in = read_input(input_json_path);
    const int samples = static_cast<int>(in.samples);
    if (in.n_abs < 1 || samples < 1 || in.rate <= 0) throw std::runtime_error("invalid n_abs, samples, or rate");
    for (int position : positions) {
        if (position < 1 || position > in.n_abs) throw std::runtime_error("position must be between 1 and n_abs");
    }

    const Eigen::MatrixXd matrix = make_matrix(in);
    Eigen::EigenSolver<Eigen::MatrixXd> solver(matrix);
    if (solver.info() != Eigen::Success) throw std::runtime_error("eigenvalue decomposition failed");
    Eigen::MatrixXcd vectors = solver.eigenvectors();
    for (int i = 0; i < vectors.cols(); ++i) {
        const auto norm = vectors.col(i).norm();
        if (std::abs(norm) > 0) vectors.col(i) /= norm;
    }
    // Keep the complex eigensystem through reconstruction.  Taking only the
    // real part of the eigenvectors loses the sine component of complex
    // conjugate modes and makes the basis rank-deficient.
    const Eigen::VectorXcd eigenvalues = solver.eigenvalues();
    const auto time = linspace(0, in.samples / in.rate, samples);
    const double c_abs = in.c_abs / in.n_abs;
    constexpr double electron_charge = 1.602e-19;
    std::vector<Pulse> result;
    result.reserve(positions.size());
    for (int position : positions) {
        Eigen::VectorXcd initial = Eigen::VectorXcd::Zero(in.n_abs + 4);
        initial[position + 1] = in.energy * 1e3 * electron_charge / c_abs;
        const Eigen::VectorXcd constants = vectors.colPivHouseholderQr().solve(initial);
        result.push_back(make_pulse(in, vectors, eigenvalues, time, position, constants));
    }
    return result;
}

void write_pulses_json(
    const std::vector<Pulse>& pulses,
    const std::string& input_json_path,
    const std::string& output_path) {
    if (pulses.empty()) throw std::runtime_error("cannot write an empty pulse set");
    std::ofstream out(output_path);
    if (!out) throw std::runtime_error("cannot open pulse output: " + output_path);
    out << std::setprecision(17) << "{\"input\":" << read_input_document(input_json_path).dump()
        << ",\"time\":";
    write_array(out, pulses.front().time);
    out << ",\"pulses\":{";
    for (std::size_t i = 0; i < pulses.size(); ++i) {
        if (i) out << ',';
        out << '"' << pulses[i].position << "\":";
        write_pulse_channels_json(out, pulses[i]);
    }
    out << "}}\n";
}

void generate_pulses_json(
    const std::string& input_json_path,
    const std::vector<int>& positions,
    const std::string& output_path,
    std::size_t thread_count) {
    const Input in = read_input(input_json_path);
    const int samples = static_cast<int>(in.samples);
    if (in.n_abs < 1 || samples < 1 || in.rate <= 0) throw std::runtime_error("invalid n_abs, samples, or rate");
    for (int position : positions) {
        if (position < 1 || position > in.n_abs) throw std::runtime_error("position must be between 1 and n_abs");
    }

    const Eigen::MatrixXd matrix = make_matrix(in);
    Eigen::EigenSolver<Eigen::MatrixXd> solver(matrix);
    if (solver.info() != Eigen::Success) throw std::runtime_error("eigenvalue decomposition failed");
    Eigen::MatrixXcd vectors = solver.eigenvectors();
    for (int i = 0; i < vectors.cols(); ++i) {
        const auto norm = vectors.col(i).norm();
        if (std::abs(norm) > 0) vectors.col(i) /= norm;
    }
    const Eigen::VectorXcd eigenvalues = solver.eigenvalues();
    const auto time = linspace(0, in.samples / in.rate, samples);
    const double c_abs = in.c_abs / in.n_abs;
    constexpr double electron_charge = 1.602e-19;
    const Eigen::ColPivHouseholderQR<Eigen::MatrixXcd> qr(vectors);
    std::vector<Eigen::VectorXcd> constants;
    constants.reserve(positions.size());
    for (int position : positions) {
        Eigen::VectorXcd initial = Eigen::VectorXcd::Zero(in.n_abs + 4);
        initial[position + 1] = in.energy * 1e3 * electron_charge / c_abs;
        constants.push_back(qr.solve(initial));
    }

    TemporaryDirectory temporary = make_temporary_directory(output_path);
    std::vector<std::filesystem::path> fragments;
    fragments.reserve(positions.size());
    for (std::size_t i = 0; i < positions.size(); ++i)
        fragments.push_back(temporary.path / (std::to_string(i) + ".json"));

    const std::size_t hardware_threads = std::thread::hardware_concurrency();
    const std::size_t requested_threads = thread_count == 0 ? (hardware_threads == 0 ? 1 : hardware_threads) : thread_count;
    const std::size_t worker_count = std::min(requested_threads, std::max<std::size_t>(1, positions.size()));
    std::atomic<std::size_t> next{0};
    std::atomic<bool> cancelled{false};
    std::exception_ptr worker_error;
    std::mutex error_mutex;
    auto worker = [&] {
        while (!cancelled.load(std::memory_order_relaxed)) {
            const std::size_t index = next.fetch_add(1, std::memory_order_relaxed);
            if (index >= positions.size()) return;
            try {
                const Pulse pulse = make_pulse(in, vectors, eigenvalues, time, positions[index], constants[index]);
                std::ofstream fragment(fragments[index]);
                if (!fragment) throw std::runtime_error("cannot write temporary pulse fragment");
                fragment << std::setprecision(17);
                write_pulse_channels_json(fragment, pulse);
                if (!fragment) throw std::runtime_error("cannot finish temporary pulse fragment");
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!worker_error) worker_error = std::current_exception();
                cancelled.store(true, std::memory_order_relaxed);
            }
        }
    };
    std::vector<std::thread> workers;
    workers.reserve(worker_count);
    for (std::size_t i = 0; i < worker_count; ++i) workers.emplace_back(worker);
    for (auto& thread : workers) thread.join();
    if (worker_error) std::rethrow_exception(worker_error);

    std::ofstream out(output_path);
    if (!out) throw std::runtime_error("cannot open pulse output: " + output_path);
    out << std::setprecision(17) << "{\"input\":" << read_input_document(input_json_path).dump()
        << ",\"time\":";
    write_array(out, time);
    out << ",\"pulses\":{";
    for (std::size_t i = 0; i < fragments.size(); ++i) {
        if (i) out << ',';
        out << '"' << positions[i] << "\":";
        std::ifstream fragment(fragments[i]);
        if (!fragment) throw std::runtime_error("cannot read temporary pulse fragment");
        out << fragment.rdbuf();
        if (!out) throw std::runtime_error("cannot write pulse output: " + output_path);
    }
    out << "}}\n";
}

}  // namespace tes_cpp
