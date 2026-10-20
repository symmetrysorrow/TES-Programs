# Research branch push and CI report

## Pre-push review

- Target branch: `feat/dump2event-physical-validation-20261020`.
- The branch started at remote commit `5ec9264ca090324a6562b358197b4930d1b75c07`.
- Existing implementation work was preserved in this isolated worktree; a
  local backup ref was also created before edits.
- The default `dump2event` mode and existing signed `E_deposit` formula are
  unchanged. The experimental algorithm is a separate executable.
- No PHITS source, PHITS executable, Dumpall, generated build output, or
  stopping-power table is included in the intended commit. The table SHA is
  documented, but the nonpublic table itself stays local.
- No local absolute paths are embedded in source. Runtime inputs are explicit
  CLI paths and manifests hash the supplied files/configuration.
- The staged allowlist contains 24 source, documentation, script and small
  validation-report files. A staged-name scan found no Dumpall, executable,
  DLL/library, PHITS Fortran source, PHITS input/output or stopping table.
  `git diff --cached --check` passed; a credential-pattern scan found no
  credential-like strings. The three small comparison CSVs contain only the
  eight boundary rows, 1,000 per-history totals, and 150 histogram bins.
- MSVC 19.51 build passed for all requested targets; CTest passed 1/1; Python
  tests passed 23/23; default converter regression processed 1,000 histories.
- PHITS same-seed internal comparison matched 201,971/201,971 additions and
  all 150 histogram bins after the NCOL=10 boundary branch was added.
- PHITS source/data exclusions, secrets scan and final staged-file allowlist
  are to be recorded after staging.

## Push and Actions

Commit and push are pending. GitHub CLI is not installed; Git remote read
access works. The commit will be pushed directly with Git after a final staged
diff review. Actions status will be queried for the pushed SHA using the
GitHub integration. PR #6 is not the target branch and remains untouched; this
experimental change is not approved for integration into PR #6.
