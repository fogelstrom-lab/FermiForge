#!/usr/bin/env bash

set -euo pipefail

readonly script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly project_root="$(cd "${script_directory}/.." && pwd)"
readonly default_branch="${FERMIFORGE_BRANCH:-main}"
readonly remote_name="${FERMIFORGE_REMOTE_NAME:-origin}"
readonly default_remote_url="git@gitlab.com:fogelstrom-lab/FermiForge.git"
readonly maximum_file_bytes="${FERMIFORGE_MAX_FILE_BYTES:-26214400}"

branch="${default_branch}"
remote_url="${FERMIFORGE_REMOTE_URL:-${default_remote_url}}"
commit_message=""
run_tests=true
push_changes=true
dry_run=false
replace_remote=false
assume_yes=false

usage() {
  cat <<'EOF'
Usage: tools/upload_to_git.sh [options]

Safely test, commit, and upload the FermiForge source tree.

Default destination:
  git@gitlab.com:fogelstrom-lab/FermiForge.git

Options:
  -m, --message TEXT       Commit message. A dated message is used by default.
      --remote-url URL     Override the Git remote URL.
      --github             Use git@github.com:fogelstrom-lab/FermiForge.git.
      --gitlab             Use git@gitlab.com:fogelstrom-lab/FermiForge.git.
      --https              Use the GitLab HTTPS URL instead of SSH.
      --replace-remote     Replace an existing origin with --remote-url.
      --branch NAME        Branch to upload (default: main).
      --skip-tests         Do not run the CMake/CTest verification gate.
      --no-push            Commit locally but do not contact or push the remote.
      --dry-run            Run checks and preview files; do not stage or commit.
  -y, --yes                Do not ask for final interactive confirmation.
  -h, --help               Show this help.

Environment overrides:
  FERMIFORGE_REMOTE_URL, FERMIFORGE_BRANCH, FERMIFORGE_REMOTE_NAME,
  FERMIFORGE_FC, FERMIFORGE_BUILD_JOBS, FERMIFORGE_MAX_FILE_BYTES.

Examples:
  tools/upload_to_git.sh -m "Initial FermiForge import"
  tools/upload_to_git.sh --https -m "Initial FermiForge import"
  tools/upload_to_git.sh -m "Add 2D trajectory regression"
EOF
}

fail() {
  printf 'error: %s\n' "$*" >&2
  exit 1
}

while (($# > 0)); do
  case "$1" in
    -m|--message)
      (($# >= 2)) || fail "$1 requires an argument"
      commit_message="$2"
      shift 2
      ;;
    --remote-url)
      (($# >= 2)) || fail "$1 requires an argument"
      remote_url="$2"
      shift 2
      ;;
    --github)
      remote_url="git@github.com:fogelstrom-lab/FermiForge.git"
      shift
      ;;
    --gitlab)
      remote_url="git@gitlab.com:fogelstrom-lab/FermiForge.git"
      shift
      ;;
    --https)
      remote_url="https://gitlab.com/fogelstrom-lab/FermiForge.git"
      shift
      ;;
    --replace-remote)
      replace_remote=true
      shift
      ;;
    --branch)
      (($# >= 2)) || fail "$1 requires an argument"
      branch="$2"
      shift 2
      ;;
    --skip-tests)
      run_tests=false
      shift
      ;;
    --no-push)
      push_changes=false
      shift
      ;;
    --dry-run)
      dry_run=true
      push_changes=false
      shift
      ;;
    -y|--yes)
      assume_yes=true
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      fail "unknown option: $1"
      ;;
  esac
done

cd "${project_root}"
[[ -f CMakeLists.txt && -f README.md && -d src && -d tests ]] || \
  fail "${project_root} does not look like the FermiForge source root"

if [[ ! -d .git ]]; then
  if "${dry_run}"; then
    fail "dry-run needs an initialized repository; run git init -b ${branch} first"
  fi
  git init -b "${branch}"
fi

current_branch="$(git symbolic-ref --quiet --short HEAD 2>/dev/null || true)"
if [[ -z "${current_branch}" ]]; then
  fail "detached HEAD is not supported by the upload script"
fi
if [[ "${current_branch}" != "${branch}" ]]; then
  fail "current branch is ${current_branch}; switch to ${branch} or pass --branch ${current_branch}"
fi

planned_remote_url="${remote_url}"
if "${push_changes}" && git remote get-url "${remote_name}" >/dev/null 2>&1; then
  configured_url="$(git remote get-url "${remote_name}")"
  if [[ "${configured_url}" != "${remote_url}" && "${replace_remote}" != true ]]; then
    fail "${remote_name} is ${configured_url}, not ${remote_url}; use --replace-remote to change it"
  fi
  if [[ "${replace_remote}" != true ]]; then
    planned_remote_url="${configured_url}"
  fi
fi

if "${run_tests}"; then
  compiler="${FERMIFORGE_FC:-}"
  if [[ -z "${compiler}" ]]; then
    for candidate in "$(command -v gfortran 2>/dev/null || true)" \
                     /opt/homebrew/bin/gfortran \
                     /opt/local/bin/gfortran \
                     /usr/local/bin/gfortran; do
      if [[ -n "${candidate}" && -x "${candidate}" ]]; then
        compiler="${candidate}"
        break
      fi
    done
  fi
  [[ -n "${compiler}" ]] || \
    fail "no GNU Fortran compiler found; set FERMIFORGE_FC or use --skip-tests"

  if [[ -z "${DEVELOPER_DIR:-}" && -d /Library/Developer/CommandLineTools ]]; then
    export DEVELOPER_DIR=/Library/Developer/CommandLineTools
  fi

  build_directory="${project_root}/work/git-upload-check"
  strict_flags="-Wall -Wextra -Wimplicit-interface -Wconversion-extra -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace"
  printf 'Running FermiForge verification with %s\n' "${compiler}"
  cmake -S "${project_root}" -B "${build_directory}" \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_Fortran_COMPILER="${compiler}" \
    -DCMAKE_Fortran_FLAGS="${strict_flags}"
  cmake --build "${build_directory}" \
    --parallel "${FERMIFORGE_BUILD_JOBS:-10}"
  ctest --test-dir "${build_directory}" --output-on-failure
fi

printf '\nGit working-tree preview:\n'
git status --short

if "${dry_run}"; then
  printf '\nFiles Git would add or update:\n'
  git add --dry-run -A
  printf '\nDry run complete; nothing was staged, committed, or pushed.\n'
  exit 0
fi

git add -A
# Preserve legacy fixed-form source and numerical reference files byte for byte.
# Apply Git's whitespace-error gate to the maintained source and documentation,
# while excluding imported historical snapshots and benchmark reference data.
git diff --cached --check -- . \
  ':(exclude)incoming/**' \
  ':(exclude)new_src/**' \
  ':(exclude)benchmarks/**/reference/**'

staged_file_is_too_large=false
while IFS= read -r -d '' path; do
  [[ -f "${path}" ]] || continue
  if stat -f '%z' "${path}" >/dev/null 2>&1; then
    file_bytes="$(stat -f '%z' "${path}")"
  else
    file_bytes="$(stat -c '%s' "${path}")"
  fi
  if ((file_bytes > maximum_file_bytes)); then
    printf 'error: staged file exceeds %s bytes: %s (%s bytes)\n' \
      "${maximum_file_bytes}" "${path}" "${file_bytes}" >&2
    staged_file_is_too_large=true
  fi
done < <(git diff --cached --name-only --diff-filter=ACMR -z)
"${staged_file_is_too_large}" && \
  fail "move large run data to external storage or raise FERMIFORGE_MAX_FILE_BYTES deliberately"

printf '\nStaged change summary:\n'
git diff --cached --stat

if "${push_changes}"; then
  printf '\nPlanned destination: %s (%s branch)\n' \
    "${planned_remote_url}" "${branch}"
else
  printf '\nNo remote push will be made.\n'
fi

if ! "${assume_yes}"; then
  [[ -t 0 ]] || fail "interactive confirmation is unavailable; rerun with --yes"
  printf 'Continue with the commit and%s? [y/N] ' \
    "$(if "${push_changes}"; then printf ' push'; fi)"
  read -r reply
  case "${reply}" in
    y|Y|yes|YES|Yes) ;;
    *)
      printf 'Stopped before committing or contacting the remote. Staged files remain available for inspection.\n'
      exit 0
      ;;
  esac
fi

if ! git diff --cached --quiet; then
  [[ -n "$(git config --get user.name || true)" ]] || \
    fail "Git user.name is not configured"
  [[ -n "$(git config --get user.email || true)" ]] || \
    fail "Git user.email is not configured"
  if [[ -z "${commit_message}" ]]; then
    commit_message="FermiForge update $(date '+%Y-%m-%d %H:%M %Z')"
  fi
  git commit -m "${commit_message}"
else
  printf 'No staged changes; no new commit was created.\n'
fi

git rev-parse --verify HEAD >/dev/null 2>&1 || \
  fail "there is no commit to upload"

if ! "${push_changes}"; then
  printf 'Local commit complete; remote upload was skipped.\n'
  exit 0
fi

if git remote get-url "${remote_name}" >/dev/null 2>&1; then
  configured_url="$(git remote get-url "${remote_name}")"
  if [[ -n "${remote_url}" && "${configured_url}" != "${remote_url}" ]]; then
    if "${replace_remote}"; then
      git remote set-url "${remote_name}" "${remote_url}"
      configured_url="${remote_url}"
    else
      fail "${remote_name} is ${configured_url}, not ${remote_url}; use --replace-remote to change it"
    fi
  fi
else
  [[ -n "${remote_url}" ]] || fail \
    "no ${remote_name} remote; pass --github, --gitlab, or --remote-url URL"
  git remote add "${remote_name}" "${remote_url}"
  configured_url="${remote_url}"
fi

printf 'Checking remote %s (%s)\n' "${remote_name}" "${configured_url}"
git ls-remote "${remote_name}" >/dev/null
if git ls-remote --exit-code --heads "${remote_name}" \
     "refs/heads/${branch}" >/dev/null 2>&1; then
  git fetch --no-tags "${remote_name}" "${branch}"
  if ! git merge-base --is-ancestor "${remote_name}/${branch}" HEAD; then
    fail "remote ${branch} is not an ancestor of local HEAD; reconcile it explicitly (the script never force-pushes)"
  fi
fi

git push --set-upstream "${remote_name}" "${branch}"
printf 'Uploaded %s to %s without rewriting remote history.\n' \
  "${branch}" "${configured_url}"
