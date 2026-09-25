#!/usr/bin/env bash
set -euo pipefail

app_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
output_root="$app_root/out"
host_binary="$output_root/host/ProteinInMembrane.Host.dll"
worker_python="$output_root/python/bin/python"

if [[ ! -f "$host_binary" || ! -x "$worker_python" ]]; then
    echo "The local product is not built. Run scripts/build-local.sh after prerequisites are supplied." >&2
    exit 1
fi

if [[ -z "${PIM_PPM_EXECUTABLE:-}" || ! -x "${PIM_PPM_EXECUTABLE:-}" ]]; then
    echo "PPM 2.0 is unavailable; placement actions will remain unavailable." >&2
fi
if [[ -z "${PIM_PACKMOL_EXECUTABLE:-}" || ! -x "${PIM_PACKMOL_EXECUTABLE:-}" ]]; then
    echo "Packmol is unavailable; explicit-system construction will remain unavailable." >&2
fi

port="${PIM_PORT:-4185}"
if [[ ! "$port" =~ ^[0-9]{1,5}$ ]] || (( port < 1024 || port > 65535 )); then
    echo "PIM_PORT must be an unprivileged TCP port number." >&2
    exit 1
fi

workspace_root="${PIM_WORKSPACE_ROOT:-$output_root/workspace}"
mkdir -p "$workspace_root"
policy_catalogue="${PIM_POLICY_CATALOGUE:-$app_root/config/policies/protein-slice1.json}"
if [[ ! -f "$policy_catalogue" ]]; then
    echo "No local policy catalogue is available; policy-bound scientific actions will remain unavailable." >&2
fi

export PIM_WORKER_PYTHON="$worker_python"
export PIM_OWNER_SOURCE_ROOT="$app_root/src/ProteinInMembrane.Host"
export PIM_WORKSPACE_ROOT="$workspace_root"
export PIM_POLICY_CATALOGUE="$policy_catalogue"
export ASPNETCORE_URLS="http://127.0.0.1:$port"

echo "Opening local protein-in-membrane workspace at $ASPNETCORE_URLS"
exec dotnet "$host_binary"
