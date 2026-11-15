# Preparing the complete research archive

`prepare_v12_archives.py` creates unpublished code and data ZIP files, checks
their CRCs and records their SHA-256 digests. It does not deposit data, register
a DOI or declare a release gate passed. Use a clean, reviewed Git snapshot.
A local Git bundle can transfer that exact commit to the execution server
without storing GitHub credentials there.

The default study directories are `evidence`, `case_studies`, `manuscript`,
`formal_execution` and `performance_workflows`. Large formal results may live
outside this local study layout. Include each of those directories explicitly:

```sh
python scripts/prepare_v12_archives.py \
  --repository /project/source/RustBEARS \
  --study-directory /project/study \
  --gate-status /project/study/release_gate_status.json \
  --research-source formal_suite=/project/results/formal_suite \
  --research-source workflow_benchmark=/project/results/workflow_benchmark \
  --research-source matched_benchmark=/project/results/matched_benchmark \
  --output-directory /project/archives/frozen_snapshot
```

Use unique archive directory names. Missing explicit sources, unsafe names and
an output directory inside an included source are rejected. Files are streamed
into ZIP64 entries and hashed from the bytes actually archived; detected source
changes abort preparation. Run only after the selected computations and figure
exports have stopped writing. Symlinks, Python caches and writing-strategy notes
are excluded. Original input licenses and attribution remain authoritative.

Keep the raw results and large archives on the execution server. A checksum
inventory or a local summary alone is not the full research-data archive.
The archive manifest records every included file and all selected source roots;
check that it contains the full formal suite, both benchmark types, input trees,
all retained failure records and final figure source data. Historical development
archives remain unchanged. Prepare a new immutable snapshot after the final
manuscript, figures and independent review are complete.
