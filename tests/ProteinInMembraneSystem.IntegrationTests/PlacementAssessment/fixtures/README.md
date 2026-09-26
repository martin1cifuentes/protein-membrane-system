# Prepared 6QWR PPM crossing fixture

`6qwr-prepared.pdb` is the first model of RCSB [6QWR](https://www.rcsb.org/structure/6qwr),
prepared through this product's accepted Slice 1 local host and OpenMM 8.6 path.
The downloaded RCSB PDB SHA-256 was
`f4c1503a60321c0cfe513e8211e43e20e2199aa5f71f22429e0e0ae8c97b0779`.
The actor approved the existing HID state for HIS 108 and the missing terminal
OXT on GLU 211; no other chemical or structural change was approved.
The resulting 3,205-atom fixture SHA-256 is
`80e2dcade32491cd49f86750ea101b050899d900ae72190f5b97f79b2052abc6`.
The companion preparation bond-graph SHA-256 was
`c3ad823d50c7088490c998e0a094d09e04c11ca4f84563e132525bc57ab37036`.

This fixed fixture tests the provider crossing and all-atom coordinate
correspondence. Fresh protein preparation is tested separately. Different
OpenMM runs can differ by a few thousandths of an angstrom in terminal H2/H3,
so this fixture's byte identity is not a general policy witness for all future
preparations of the same source.

The local PPM executable and its GPL-licensed `res.lib` are installed
separately under ignored `out/ppm2/`; neither is part of this fixture.
