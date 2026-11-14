# Security Policy

## Supported Versions

| Version | Supported |
| ------- | --------- |
| 1.1.x   | Yes       |

## Reporting a Vulnerability

If you discover a memory-safety defect across the C-ABI boundary (`bgb-ffi`, `r_package/RustBEARS/src/init.c`, or `rustbears/_ffi.py`), arbitrary code execution risk in input parsers, or any other security vulnerability:

1. Please **do not** post proof-of-concept crash payloads publicly before a fix is available.
2. Report the issue privately via GitHub Security Advisories on <https://github.com/Jhe1004/RustBEARS/security/advisories> or email the maintainer (`j.he930724@gmail.com`) with reproduction steps and platform details.
3. We aim to acknowledge reports within 72 hours and ship a patched release with regression tests promptly.
