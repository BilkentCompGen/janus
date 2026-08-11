# svp code style — house rules

Non-negotiable on this branch. Every line is read and reviewed; brevity and directness beat "clean
architecture".

## Hard rules

1. **C++17.** No newer standard features unless a hard dependency forces it.
2. **No OOP.** No classes except (a) POD structs and (b) trivial RAII wrappers around C resources we cannot otherwise free. No inheritance, no virtual, no member functions on POD types, no ctors/dtors beyond the RAII case.
3. **No user-written templates.** STL templates (`std::vector<T>`, `std::string_view`) are fine — that's the price of using the STL.
4. **No fancy C++.** No inline function pointers; no clever `using` aliases; no `auto` for non-iterator return types; no `constexpr` cleverness; no `std::variant`/`std::optional` if a sentinel does the job.
5. **Fail loudly, fail fast.** No `try`/`catch`. No log-and-continue. Any unexpected condition: `fprintf(stderr, "[FATAL] ...\n"); std::abort();`. Do **not** add fallback paths "just in case".
6. **No defensive coding.** Trust internal callers — only validate at the CLI boundary. Callee functions assume their preconditions hold.
7. **Avoid if-else trees.** Prefer early return; arithmetic / table lookup; `std::min`/`std::max`; the `cond ? a : b` ternary for trivial selections. More than two chained `else if` means stop and rethink.
8. **Function size ≤ 30 LOC.** Single-purpose. Split if it grows.
9. **Every stage stays small enough to read in one sitting.** The whole binary is ~2 kLOC of our own code, excluding the vendored `unordered_dense.h`. A stage that cannot be explained in its own header comment is too big.
10. **Comments only when the WHY is non-obvious.** Never restate WHAT — well-named identifiers do that. A comment that says what a measurement was and what it cost is worth keeping; "added for issue X" is not.

## What's allowed

- POD structs with `operator==` if needed for hashtable keys (only when STL forces it).
- File-scope `static` helpers, and one anonymous namespace per `.cpp` for the stage's own internals.
- `std::string_view` for non-owning slices; never store owning strings when a view will do.
- Buffer reuse via `clear()` + `reserve()` to avoid repeated allocs.
- Minimal `#include` set per file. No catch-all umbrella headers.

## Coordinates

- Read/contig positions are **0-based half-open** (matches lcptools' `core.start/end`, end exclusive).
- VCF positions are 1-based from the moment a junction is built, so the writer converts nothing.

## File organisation

- `svp.h` holds what more than one stage needs: the extraction decl, the shared PODs, the reads
  interface, the checkpoint IO helpers. Nothing else goes in it.
- Each stage is one `.cpp` plus a header that carries its one-line entry point and the explanation of
  what the stage is for.

## Build

- One `Makefile`. No CMake, no CI, no test harness. A single binary, `svp`.
- `make clean && make` must succeed with only `lib/lcptools` built (`make install-deps-lcptools`) and
  system htslib present.
