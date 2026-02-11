#!/usr/bin/env bash
# Run test_old_vs_new.py in the background with nohup.
# Output goes to test_old_vs_new.log in the current directory.
#
# Usage:
#   bash scripts/run_test_old_vs_new.sh [-- extra args for test_old_vs_new.py]
#
# Examples:
#   bash scripts/run_test_old_vs_new.sh
#   bash scripts/run_test_old_vs_new.sh -- --geno /path/to/data --out /path/to/output

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOGFILE="${SCRIPT_DIR}/../test_old_vs_new.log"

echo "Starting test_old_vs_new.py in background..."
echo "Log file: ${LOGFILE}"
echo "Monitor with: tail -f ${LOGFILE}"

nohup python "${SCRIPT_DIR}/test_old_vs_new.py" "$@" > "${LOGFILE}" 2>&1 &

echo "PID: $!"
