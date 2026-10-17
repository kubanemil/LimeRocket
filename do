#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

examples=(examples/*.yaml)

case "${1:-help}" in
  test)
    shift
    uv run pytest "$@"
    ;;
  check)
    uv run ruff check .
    uv run ruff format --check .
    uv run mypy
    ;;
  fmt)
    uv run ruff check --fix .
    uv run ruff format .
    ;;
  run)
    shift
    uv run limerocket "$@"
    ;;
  examples)
    for spec in "${examples[@]}"; do
      name="$(basename "$spec" .yaml)"
      uv run limerocket design "$spec" -o "examples/out/$name"
    done
    ;;
  help|*)
    cat <<'USAGE'
usage: ./do <task>

  test [args]   run pytest (extra args go to pytest)
  check         ruff lint, ruff format check, mypy
  fmt           apply ruff fixes and formatting
  run [args]    run the limerocket CLI, e.g. ./do run design examples/lime-1.yaml -o out
  examples      regenerate the example reports and exports in examples/out
USAGE
    ;;
esac
