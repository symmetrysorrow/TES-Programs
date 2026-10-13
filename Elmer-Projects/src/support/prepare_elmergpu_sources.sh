#!/usr/bin/env bash
# Fetch pinned upstream trees and apply the source patches tracked by this repo.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(realpath "$SCRIPT_DIR/../..")"
PATCH_DIR="$REPO_ROOT/docs"
LOCK="$PATCH_DIR/elmergpu-sources.lock"
SOURCE_ROOT="${1:-/root/src}"

# The lock file is shell-compatible key=value data and is the single source
# of truth for upstream URLs, base commits, and patch paths.
# shellcheck disable=SC1090
source "$LOCK"
ELMER_URL="$elmer_repository"
ELMER_BASE="$elmer_base"
ELMER_PATCH="$REPO_ROOT/$elmer_patch"
HYPRE_URL="$hypre_repository"
HYPRE_BASE="$hypre_base"
HYPRE_PATCH="$REPO_ROOT/$hypre_patch"

echo "$elmer_patch_sha256  $ELMER_PATCH" | sha256sum --check --status || {
  echo "Elmer source patch checksum does not match the lock file" >&2
  exit 2
}
echo "$hypre_patch_sha256  $HYPRE_PATCH" | sha256sum --check --status || {
  echo "HYPRE source patch checksum does not match the lock file" >&2
  exit 2
}

mkdir -p "$SOURCE_ROOT"

checkout_source() {
  local url="$1" base="$2" patch="$3" target="$4"
  if [[ ! -d "$target/.git" ]]; then
    if [[ -e "$target" ]] && [[ -n "$(find "$target" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
      echo "Refusing to reuse non-Git source directory: $target" >&2
      exit 2
    fi
    git clone "$url" "$target"
  fi
  local status head changed_files expected_files
  status="$(git -C "$target" status --porcelain)"
  head="$(git -C "$target" rev-parse HEAD)"
  if [[ "$head" != "$base" ]]; then
    if [[ -n "$status" ]]; then
      echo "Source tree has local changes and is not at locked base $base: $target" >&2
      exit 2
    fi
    git -C "$target" fetch --no-tags origin "$base"
    git -C "$target" checkout --detach "$base"
    status="$(git -C "$target" status --porcelain)"
  fi
  if [[ -n "$status" ]]; then
    if [[ -n "$(git -C "$target" ls-files --others --exclude-standard)" ]]; then
      echo "Source tree has untracked files: $target" >&2
      exit 2
    fi
    if ! git -C "$target" diff --cached --quiet; then
      echo "Source tree has staged changes: $target" >&2
      exit 2
    fi
    # The locked patches are intentionally generated with --unidiff-zero.
    # `git diff --binary HEAD` adds context and therefore has a different byte
    # hash even when the patch is applied exactly. Verify the patch in reverse
    # and require that it accounts for every modified source file.
    if ! git -C "$target" apply --unidiff-zero --reverse --check "$patch" >/dev/null 2>&1; then
      echo "Source tree has edits that do not match its locked patch: $target" >&2
      exit 2
    fi
    changed_files="$(git -C "$target" diff --name-only HEAD | sort)"
    expected_files="$(sed -n 's#^diff --git a/\([^ ]*\) b/.*#\1#p' "$patch" | sort)"
    if [[ "$changed_files" != "$expected_files" ]]; then
      echo "Source tree has modified files outside its locked patch: $target" >&2
      exit 2
    fi
  else
    git -C "$target" fetch --no-tags origin "$base"
    git -C "$target" checkout --detach "$base"
  fi
}

apply_patch_once() {
  local target="$1" patch="$2"
  if git -C "$target" apply --unidiff-zero --reverse --check "$patch" >/dev/null 2>&1; then
    echo "Patch already applied: $(basename "$patch")"
    return
  fi
  git -C "$target" apply --unidiff-zero --check "$patch"
  git -C "$target" apply --unidiff-zero "$patch"
}

checkout_source "$ELMER_URL" "$ELMER_BASE" "$ELMER_PATCH" "$SOURCE_ROOT/elmer"
checkout_source "$HYPRE_URL" "$HYPRE_BASE" "$HYPRE_PATCH" "$SOURCE_ROOT/hypre"
apply_patch_once "$SOURCE_ROOT/elmer" "$ELMER_PATCH"
apply_patch_once "$SOURCE_ROOT/hypre" "$HYPRE_PATCH"

printf 'Elmer base: %s\nElmer patch: %s\nHYPRE base: %s\nHYPRE patch: %s\n' \
  "$ELMER_BASE" "$(sha256sum "$ELMER_PATCH" | cut -d ' ' -f 1)" \
  "$HYPRE_BASE" "$(sha256sum "$HYPRE_PATCH" | cut -d ' ' -f 1)"
