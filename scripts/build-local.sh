#!/usr/bin/env bash
set -euo pipefail

app_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
browser_root="$app_root/browser"
owner_source_root="$app_root/src/ProteinInMembrane.Host"
host_project="$app_root/src/ProteinInMembrane.Host/ProteinInMembrane.Host.csproj"
architecture_project="$app_root/tools/ArchitectureCheck/ArchitectureCheck.csproj"
output_root="$app_root/out"

for command_name in node npm dotnet python3.11; do
    if ! command -v "$command_name" >/dev/null 2>&1; then
        echo "Required build tool is unavailable: $command_name" >&2
        exit 1
    fi
done

dotnet run --project "$architecture_project" --configuration Release -- "$app_root/src/ProteinInMembrane.Host"

(
    cd "$browser_root"
    if [[ ! -f package-lock.json ]]; then
        npm install --package-lock-only --ignore-scripts
        echo "Generated package-lock.json. Review and commit it with the first verified build." >&2
    fi
    npm ci
    npm run build
)

dotnet publish "$host_project" --configuration Release --output "$output_root/host"

python3.11 -m venv "$output_root/python"
"$output_root/python/bin/python" -m pip install --requirement "$owner_source_root/ProteinInMembraneSystem/worker/requirements.txt"
"$output_root/python/bin/python" -m pip freeze --all > "$output_root/python-resolved-requirements.txt"

if [[ ! -f "$output_root/host/wwwroot/index.html" ]]; then
    echo "The browser bundle was not included in the published host." >&2
    exit 1
fi

echo "Local build prepared at $output_root"
