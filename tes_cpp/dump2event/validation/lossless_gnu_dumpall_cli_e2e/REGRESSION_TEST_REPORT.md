# Regression test report

No source code in the PHITS installation or production converter was changed. On 2026-10-21: `cmake --build build/physdeposit-vs --config Release -j 4` succeeded; `ctest --test-dir build/physdeposit-vs -C Release --output-on-failure` passed 1/1; `python -m pytest tes_cpp -q` reported 44 passed and 10 skipped. The default converter summary-only smoke of 1,000 histories is a historical result, not rerun today. The synthetic GNU Fortran line-copy probe compiled and ran today, but failed exact preservation of an unterminated final line and is not a PHITS integration test.

A fresh full C++/Python regression suite was not run because no implementation change was made; it should be rerun before any future code commit.
