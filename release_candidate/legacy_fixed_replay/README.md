# Minimal fixed-prediction replay

All 22 required input files are included under source-specific conditions. Seventeen are unchanged originals; five are exact-token column projections. Original enriched research tables are preserved outside the public tree.

Use the safe runner from the repository root, after installing the root dependency lock:

```sh
python tools/run_fixed_replay.py --bundle release_candidate/legacy_fixed_replay --output-dir liver_replay_results --timeout-seconds 1200
```

Expected result: 238 checks at absolute tolerance 2e-12, zero new random draws and zero new model queries. Old output is excluded and a fresh working directory is used. A direct `python replay.py` writes inside this bundle and requires a writable copy; it is not the preferred entry. Neither command is itself an OS sandbox.

See [data access and retained grouping](../../DATA_ACCESS_AND_REPLAY.md), [22 input decisions](../../sources/REPLAY_INPUT_DECISIONS.csv), [field inventory](field_inventory.csv), [data dictionary](data_dictionary.md), [provenance](provenance.csv), [group membership](group_membership.csv) and [validation](../../VALIDATION.md). `provenance.csv` gives both current projected bytes and original bytes. Historical code and training are separate from this fixed-score layer.
