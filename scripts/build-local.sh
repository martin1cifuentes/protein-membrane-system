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
worker_python="$output_root/python/bin/python"
worker_requirements="$owner_source_root/ProteinInMembraneSystem/worker/requirements.txt"
if ! "$worker_python" - "$worker_requirements" <<'PY'
import importlib.metadata as metadata
import json
from pathlib import Path
import sys

for raw in Path(sys.argv[1]).read_text().splitlines():
    requirement = raw.strip()
    if not requirement or requirement.startswith("#"):
        continue
    try:
        if " @ git+" in requirement:
            name, source = requirement.split(" @ git+", 1)
            url, commit = source.rsplit("@", 1)
            installed = json.loads(metadata.distribution(name).read_text("direct_url.json") or "{}")
            if installed.get("url") != url or installed.get("vcs_info", {}).get("commit_id") != commit:
                sys.exit(1)
        else:
            name, version = requirement.split("==", 1)
            if metadata.version(name) != version:
                sys.exit(1)
    except (metadata.PackageNotFoundError, ValueError, json.JSONDecodeError):
        sys.exit(1)
PY
then
    "$worker_python" -m pip install --requirement "$worker_requirements"
fi
"$worker_python" -m pip freeze --all > "$output_root/python-resolved-requirements.txt"

ambertools_root="${PIM_AMBERTOOLS_HOME:-$output_root/ambertools26}"
ambertools_manifest="$app_root/config/providers/ambertools-environment.yml"
if [[ ! -x "$ambertools_root/bin/packmol-memgen" ]]; then
    if command -v micromamba >/dev/null 2>&1; then
        if ! micromamba create --yes --prefix "$ambertools_root" --file "$ambertools_manifest" --strict-channel-priority; then
            echo "Memgen runtime could not be prepared; its route will be reported unavailable." >&2
        fi
    else
        echo "Memgen runtime unavailable: micromamba is needed to prepare $ambertools_manifest" >&2
    fi
fi
if [[ -x "$ambertools_root/bin/python" ]]; then
    if ! "$ambertools_root/bin/python" - "$ambertools_root" <<'PY'
import importlib.metadata as metadata
from pathlib import Path
import sys

root = Path(sys.argv[1])
if metadata.version("packmol-memgen") != "2026.3.25":
    raise SystemExit("Memgen runtime version differs from pinned 2026.3.25")
if not all((root / "bin" / name).is_file() for name in
           ("packmol-memgen", "packmol", "tleap", "sander", "ambpdb")):
    raise SystemExit("Memgen runtime lacks a required CPU executable")
PY
    then
        echo "Installed Memgen runtime does not match the selected route; other routes remain independent." >&2
    fi
    if command -v micromamba >/dev/null 2>&1; then
        micromamba list --prefix "$ambertools_root" --json > "$output_root/ambertools26-resolved-packages.json"
    fi
fi

if [[ ! -f "$output_root/host/wwwroot/index.html" ]]; then
    echo "The browser bundle was not included in the published host." >&2
    exit 1
fi

echo "Local build prepared at $output_root"
