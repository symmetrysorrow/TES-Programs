
#include <exception>
#include <iostream>
#include <string>
#include <filesystem>
#include <cstdlib>
#include <vector>
#ifdef _WIN32
#include <process.h>
#else
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace {
void usage() {
    std::cerr << "Usage (PHITS 3.341 / Pb EGS5, HDF5 v2 only):\n"
                 "  dump2event <dumpall> <output.h5>\n"
                 "    --physics-config <pb_config.json> [--stopping-table <pb_egs5.csv>]\n"
                 "    [--physics pb-egs5] [--hdf5-version 2]\n"
                 "  The legacy reaction-difference formula and JSON/HDF5 v1 output\n"
                 "  have been removed. Only the validated Pb unit-weight contract\n"
                 "  is supported; invalid or unsupported inputs fail closed.\n";
}

std::string shell_quote(const std::string& value) {
#ifdef _WIN32
    std::string out = "\"";
    for (char c : value) { if (c == '"') out += "\\\""; else out += c; }
    out += '"';
    return out;
#else
    std::string out = "'";
    for (char c : value) { if (c == '\'') out += "'\\''"; else out += c; }
    out += '\'';
    return out;
#endif
}

std::filesystem::path find_physics_driver(const std::filesystem::path& executable) {
    if (const char* configured = std::getenv("TES_CPP_PB_PHYSICS_DRIVER")) {
        const std::filesystem::path path(configured);
        if (std::filesystem::is_regular_file(path)) return path;
        throw std::runtime_error("TES_CPP_PB_PHYSICS_DRIVER does not name a file");
    }
    auto parent = std::filesystem::absolute(executable).parent_path();
    for (int depth = 0; depth < 10 && !parent.empty(); ++depth) {
        const auto repo = parent / "tes_cpp" / "dump2event" / "validation" / "build_physical_hdf5_v2.py";
        if (std::filesystem::is_regular_file(repo)) return repo;
        const auto installed = parent / "share" / "tes_cpp_dump2event" / "build_physical_hdf5_v2.py";
        if (std::filesystem::is_regular_file(installed)) return installed;
        const auto adjacent = parent / "build_physical_hdf5_v2.py";
        if (std::filesystem::is_regular_file(adjacent)) return adjacent;
        parent = parent.parent_path();
    }
    throw std::runtime_error("Pb physics driver not found; set TES_CPP_PB_PHYSICS_DRIVER");
}

int run_pb_physics(int argc, char** argv) {
    if (argc < 3) { usage(); return 2; }
    std::string config, stopping;
    for (int i = 3; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--physics" && i + 1 < argc) {
            if (std::string(argv[++i]) != "pb-egs5") throw std::runtime_error("only --physics pb-egs5 is supported");
        } else if (arg == "--hdf5-version" && i + 1 < argc) {
            if (std::string(argv[++i]) != "2") throw std::runtime_error("Pb physics mode requires --hdf5-version 2");
        } else if (arg == "--physics-config" && i + 1 < argc) config = argv[++i];
        else if (arg == "--stopping-table" && i + 1 < argc) stopping = argv[++i];
        else throw std::runtime_error("unknown or incomplete Pb physics option: " + arg);
    }
    const std::filesystem::path input(argv[1]), output(argv[2]);
    if (!std::filesystem::is_regular_file(input)) throw std::runtime_error("Dumpall input does not exist");
    if (!config.empty() && !std::filesystem::is_regular_file(config))
        throw std::runtime_error("physics config does not exist");
    if (!stopping.empty() && !std::filesystem::is_regular_file(stopping))
        throw std::runtime_error("explicit stopping table does not exist");
    if (output.extension() != ".h5" && output.extension() != ".hdf5")
        throw std::runtime_error("Pb physical mode writes HDF5 v2; output extension must be .h5 or .hdf5");
    const auto exe = std::filesystem::absolute(argv[0]);
    auto reconstructor = exe.parent_path() / "dump2event_physical_reconstruct";
    auto exporter = exe.parent_path() / "dumpall_transport_export";
#ifdef _WIN32
    reconstructor += ".exe"; exporter += ".exe";
#endif
    if (!std::filesystem::is_regular_file(reconstructor) || !std::filesystem::is_regular_file(exporter))
        throw std::runtime_error("physical ledger/export helper binaries must be alongside dump2event");
    const auto driver = find_physics_driver(exe);
    const char* configured_python = std::getenv("TES_CPP_PYTHON");
#ifdef _WIN32
    const std::string python = configured_python ? configured_python : "python";
#else
    const std::string python = configured_python ? configured_python : "python3";
#endif
    std::vector<std::string> args{python, driver.string(), "--dumpall", input.string(),
        "--reconstructor", reconstructor.string(), "--track-exporter", exporter.string(),
        "--converter-executable", exe.string(), "--output", output.string()};
    if (!config.empty()) {
        args.emplace_back("--physics-config");
        args.push_back(config);
    }
    if (!stopping.empty()) {
        args.emplace_back("--stopping-table");
        args.push_back(stopping);
    }
#ifdef _WIN32
    std::vector<std::wstring> wide_args;
    wide_args.reserve(args.size() + 1);
    for (const auto& arg : args) wide_args.push_back(std::filesystem::path(arg).wstring());
    std::vector<const wchar_t*> argv_wide;
    argv_wide.reserve(wide_args.size() + 1);
    for (const auto& arg : wide_args) argv_wide.push_back(arg.c_str());
    argv_wide.push_back(nullptr);
    const intptr_t status = _wspawnvp(_P_WAIT, argv_wide[0], argv_wide.data());
    if (status == -1) throw std::runtime_error("failed to launch Python Pb HDF5 driver");
#else
    std::string command;
    for (const auto& arg : args) { if (!command.empty()) command += ' '; command += shell_quote(arg); }
    const int status = std::system(command.c_str());
#endif
    if (status != 0) {
        std::cerr << "dump2event: experimental Pb reconstruction failed; subprocess status=" << status << '\n';
        return 1;
    }
    return 0;
}
}

int main(int argc, char** argv) {
    if (argc < 3) {
        usage();
        return 2;
    }
    try {
        return run_pb_physics(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << "dump2event: " << error.what() << '\n';
        return 2;
    }
}
