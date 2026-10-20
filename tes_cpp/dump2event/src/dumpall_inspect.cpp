#include "tes_cpp/dumpall_record_parser.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) try {
    if (argc != 2) {
        std::cerr << "Usage: dumpall_record_inspect <ASCII dumpall.dat>\n";
        return 2;
    }
    std::ifstream input(argv[1], std::ios::binary);
    if (!input) throw std::runtime_error("cannot open Dumpall file");
    const auto summary = tes_cpp::dump2event::parse_dumpall_records(input,
        [](const tes_cpp::dump2event::DumpallRecord&) {}, argv[1]);
    std::cout << "records=" << summary.records << " histories=" << summary.histories
              << " terminator=NCOL=" << summary.terminal_ncol << '\n';
    for (std::size_t ncol = 4; ncol < summary.by_ncol.size(); ++ncol)
        if (summary.by_ncol[ncol]) std::cout << "NCOL=" << ncol << ':' << summary.by_ncol[ncol] << '\n';
    return 0;
} catch (const std::exception& error) {
    std::cerr << "dumpall_record_inspect: " << error.what() << '\n';
    return 1;
}
