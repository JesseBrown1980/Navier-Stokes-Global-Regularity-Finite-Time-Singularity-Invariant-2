#!/usr/bin/env bash
# Full local gate: clippy (warnings denied), all tests, release verifier, receipts.
# Integer-only Rust 1.81.0, offline. Writes logs/ and receipts/VERIFY.hbp only.
set -u
[ -f "$HOME/.cargo/env" ] && . "$HOME/.cargo/env"
cd "$(dirname "$0")/.."
TC="${TOOLCHAIN:-+1.81.0}"
export CARGO_CACHE_RUSTC_INFO=0
mkdir -p logs receipts
rc=0

cargo $TC clippy --offline --release --all-targets -- -D warnings >logs/clippy.log 2>&1 || rc=1
echo "clippy_exit=$rc"

cargo $TC test --offline --release --all-targets >logs/tests.log 2>&1 || rc=1
grep -E '^(running|test result)' logs/tests.log

cargo $TC run --offline --release --bin verify >receipts/VERIFY.hbp 2>logs/verify.stderr.log || rc=1
cat receipts/VERIFY.hbp
sha256sum receipts/VERIFY.hbp | awk '{print $1"  VERIFY.hbp"}' >receipts/VERIFY.hbp.sha256

echo "GATE_EXIT=$rc"
exit $rc
