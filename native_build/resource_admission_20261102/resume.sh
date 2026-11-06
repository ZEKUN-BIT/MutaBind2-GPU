#!/usr/bin/env bash
set -euo pipefail
task_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
cd "$task_root"
exec "$task_root/.env/bin/python" -u "$task_root/native_build/resource_admission_20261102/continue_scan.py"
