# svp code style — house rules

These rules are non-negotiable on this branch. The user reads every line and reviews with Codex; brevity and directness beat "clean architecture."

## Hard rules

1. **C++17.** No newer standard features unless a hard dependency forces it.
2. **No OOP.** No classes except (a) POD structs and (b) trivial RAII wrappers around C resources we cannot otherwise free. No inheritance, no virtual, no member functions on POD types, no ctors/dtors beyond the RAII case.
3. **No user-written templates.** STL templates (`std::vector<T>`, `std::unordered_map<K,V>`, `std::string_view`) are fine — that's the price of using the STL.
4. **No fancy C++.** No lambdas with non-trivial captures; no inline function pointers; no clever `using` aliases; no namespaces beyond the implicit global one; no `auto` for non-iterator return types; no `constexpr` cleverness; no `std::variant`/`std::optional` if a sentinel does the job.
5. **Fail loudly, fail fast.** No `try`/`catch`. No log-and-continue. Any unexpected condition: `fprintf(stderr, "[FATAL] ...\n"); std::abort();`. `assert()` is fine in tight inner code. Do **not** add fallback paths "just in case."
6. **No defensive coding.** Trust internal callers — only validate at the CLI boundary. Callee functions assume their preconditions hold.
7. **Avoid if-else trees.** Prefer early return; arithmetic / table lookup; `std::min`/`std::max`; the `cond ? a : b` ternary for trivial selections. If you find yourself writing more than two `else if` chained, stop and rethink.
8. **Function size ≤ 30 LOC.** Single-purpose. Split if it grows.
9. **Total binary target ≤ 700 LOC** across all `.cpp` we write (excluding vendored/system carries: `unordered_dense.h`, `xxhash.h`, and the btllib bloom headers).
10. **Comments only when the WHY is non-obvious.** Never restate WHAT — well-named identifiers do that. No comments referencing "added for issue X", "used by Y", "removed Z" — that belongs in PR descriptions / git log.

## What's allowed

- POD structs with `operator==` if needed for hashtable keys (only when STL forces it).
- File-scope `static` helpers.
- `std::string_view` for non-owning slices; never store owning strings when a view will do.
- Buffer reuse via `clear()` + `reserve()` to avoid repeated allocs.
- Minimal `#include` set per file. No catch-all umbrella headers.

## Coordinate system

- All positions are **0-based half-open** (matches lcptools' `core.start/end`).
- Bubble TSV coordinates are flat positions into the **concatenated input buffer**, not into any FASTA record. Mapping back to record coordinates is downstream.

## File organisation

- One `svp.h` for shared structs + decls. No other cross-module headers.
- Each `.cpp` owns its forward decls or pulls them from `svp.h`.
- No header-only modules we wrote ourselves. (Carry-overs like `zstr.hpp` are exceptions.)

## Build

- One `Makefile`. No CMake.
- Single binary `svp`. No installer, no test harness, no CI config.
- `make clean && make` must succeed end-to-end with only `lib/lcptools` built and zlib present.

## When in doubt

Three concrete short lines is better than one clever abstracted line. The user will read this code and edit it — optimise for *the next reader*, not for elegance.
