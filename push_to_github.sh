#!/usr/bin/env bash

set -euo pipefail

readonly script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

exec "${script_directory}/tools/upload_to_git.sh" --github "$@"
