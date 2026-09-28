# Identified POPC starting patch

`POPC-Lipid21-Zenodo-64x64.gro` is adapted from Patrick Frankel,
[*Amber Lipid21 bilayer simulations*](https://doi.org/10.5281/zenodo.14776136),
Zenodo v1 (2025). The dataset's [DataCite rights record](https://api.datacite.org/dois/10.5281/zenodo.14776136)
declares [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The source
GRO digest is `5ec34d5a49f208af18b20eaeb97276bd5737e3c8bf6ef265f20fde3718b3e702`.
The dataset describes 128-lipid GROMACS Lipid21 systems and 500 ns simulation
trajectories with their first 100 ns removed. The exact time represented by
this GRO file was not independently established. Its solvated periodic
coordinates and observed organization are used as a starting patch, not as an
equilibrium claim for a later protein-containing system.

`POPC-Lipid21-Zenodo-64x64.cif` is an **adaptation** for OpenMM 8.6
`Modeller.addMembrane`. The script
[`scripts/map-zenodo-popc-patch.py`](../../../scripts/map-zenodo-popc-patch.py)
applies the exact atom permutation in `POPC-Lipid21-Zenodo-64x64-map.json`,
joining each source lipid's three linked Amber residues into the selected
one-residue Lipid21 POPC template. It preserves every coordinate, all 128
lipids, 5,120 waters and the orthorhombic cell; it does not use GRO velocities.
The mapped CIF digest is `464ecd96499b9574b8d35e60dbf61b185e487a10296cd09aabdde4db8b2c7c37`.
The script reproduced those exact bytes. The source and adapted files are
included with this attribution; no original author-repository files are
redistributed here.

The selected POPC molecular graph and stereochemistry were checked for every
source lipid, and the mapped patch was checked against the exact selected
Lipid21/OpenMM templates and periodic cell. Parameter matching, construction,
same-candidate minimization, stage assessment and export are separate runtime
checks. A passing build does not establish biological suitability, membrane
equilibration or the deferred optional-equilibration capability.
