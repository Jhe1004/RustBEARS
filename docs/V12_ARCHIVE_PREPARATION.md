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

The research ZIP now lists every omitted file in `ARCHIVE_MANIFEST.json` and
`archive_index.json` with its original path, intended archive path, byte size,
SHA-256 and exclusion reason. Both records include the actual omitted file count
and total size. Inputs are never removed or changed. A symlink is not followed;
its recorded hash and size describe its link-target text, as explicitly labelled.

Default exclusions cover compiled/distribution suffixes (including wheels and
shared libraries), extensionless native executable headers, installation-marked R packages,
Python package trees and conda/virtualenv runtimes. Source-package directories
without installation markers remain included. Formal NPZ data, arbitrary event
`.bin` files, RDS/JSON checkpoints, complete failure records, tables, logs and
source manifests remain included unless they actually reside inside an installed
runtime/package. Do not mix authoritative research outputs into runtime trees.
The existing private `qa_toolchain/` script copies are withheld until their
redistribution license is reviewed; the known `figure_qa/private_toolchain.tar.gz`
bundle is excluded under the same rule. The adjacent QA inventories/results remain.
Unrelated scientific tar/ZIP archives, including a same-named archive outside
that private-tool context, are retained.
Private writing notes, Python caches and platform metadata are also recorded as
exclusions instead of silently omitted.

These rules apply to selected research roots, including explicit server result
roots, and leave the exact Git code archive unchanged. They are a bounded file
classification policy, not a recursive inspection of arbitrary ZIP/tar contents
or a license grant. Original dataset/source archives are retained; review their
own license and contents before deposit. Keep reproduction/build instructions
and corresponding reviewed source with the final release; excluded binaries can
be rebuilt or obtained from their recorded sources. No private tool is silently
relicensed by this archive operation.
