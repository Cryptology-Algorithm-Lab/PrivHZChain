# Privacy-Preserving Monitoring of Controlled Substance Supply Chains Using Hierarchical Blockchain and ZKP

A reference implementation of the hidden accumulator and the zero-knowledge proof systems for set relations (membership, non-membership, set-split, and aggregated membership) that power PrivHZChain, together with the Entry / Transfer / Exit / Summary transactions and the (non-)ownership proofs used for investigation.

## Requirements

- A C++17 compiler: GCC ≥ 10 or Clang (older GCC fails to build mcl's AVX-512 MSM unit; see [Building with Docker](#building-with-docker))
- CMake ≥ 3.16
- OpenMP — the setup, witness-generation, and some MSM loops in `zkacc` (manager-side) are parallelized with `#pragma omp parallel for`
- [herumi/mcl](https://github.com/herumi/mcl) for the BLS12-381 curve

### Obtaining mcl

The build pulls mcl in via `add_subdirectory(mcl)`, so you only need the source tree at `./mcl` — no separate mcl build step is required (it is compiled as part of this project). From the project root:

```bash
git clone https://github.com/herumi/mcl
```

No GMP is needed: the BLS12-381 path in mcl uses its own assembly backend.

## Building

A standard out-of-source CMake build:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

This produces the benchmark executable at `build/test`.

The project's own sources compile with `-O2` and **without** `-march=native` by default, so numbers stay comparable across machines. To let the compiler use the host's full instruction set (e.g. AVX-512 or IFMA), enable the native option:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPRIVHZCHAIN_NATIVE=ON
cmake --build build -j
```

## Building with Docker

If you would rather not install a toolchain or clone mcl by hand, build everything in a container. The provided `Dockerfile` is based on Ubuntu 24.04 (GCC 13), clones mcl at build time, and compiles both mcl and this project — nothing extra needs to live in your tree.

From the project root (with `Dockerfile` and `.dockerignore` present):

```bash
docker build -t privhzchain .
```

This produces an image with the benchmark built at `/app/build/test`. Run it with the same flags as the native binary:

```bash
# run a benchmark mode
docker run --rm privhzchain /app/build/test --mode local

# keep the CSV output on the host
docker run --rm -v "$PWD/out":/out -w /out privhzchain /app/build/test --mode all
```

Ubuntu 24.04's GCC 13 also sidesteps a build failure seen on older toolchains: GCC < 10 (e.g. the GCC 9 shipped with Ubuntu 20.04) lacks an intrinsic used by mcl's AVX-512 IFMA MSM unit (`msm_avx.cpp`), so mcl does not compile there.

## Running the benchmarks

```bash
./build/test [--mode local|server|all] [--verbose|--no-verbose]
```

- `--mode all` (default) runs every benchmark.
- `--mode local` runs the participant-/peer-side workloads.
- `--mode server` runs the manager-side workloads.
- `--verbose` echoes each measurement as it is taken; `--no-verbose` prints only the final summary table.

Each run prints an aligned summary table (mean / median / p95 / min / max per metric) and writes the same data to `bench_<mode>.csv`. Times are in microseconds; proof sizes are in bytes.

## Benchmark catalog

Benchmarks are grouped by which entity performs the work in the protocol. Generation and verification can be done by different entities (e.g. a participant generates a transaction, a peer verifies it); the grouping reflects protocol-level roles, not a single machine.

### `local` — Permitted Participants (generation) and Blockchain Peers (verification)

| Test            | What it measures                                                                               |
| --------------- | ---------------------------------------------------------------------------------------------- |
| `ZKMP_test`     | ZK membership proof (`ΠMP`) — prove and verify                                                 |
| `ZKNMP_test`    | ZK non-membership proof (`ΠNMP`) — prove and verify                                            |
| `ZKSP_test`     | ZK set-split proof (`ΠSP`, the disjoint-union / quantity-conservation core) — prove and verify |
| `OwnPf_test`    | Ownership proof for investigation — prove and verify                                           |
| `NonOwnPf_test` | Non-ownership proof for dynamic suspension — prove and verify                                  |
| `TxGenTest`     | Entry / Transfer / Exit transactions — full TxCom + TxPrv and verify                           |

### `server` — Global Manager (setup) and Local Managers (summary aggregation)

| Test            | What it measures                                                                        |
| --------------- | --------------------------------------------------------------------------------------- |
| `AZKMP_test`    | Aggregated membership proof (`ΠAMP`) — prove and verify (the Summary Tx building block) |
| `TxSummaryTest` | Summary transaction generation (local manager) and verification                         |
| `SetupTest`     | Global-manager setup: public-parameter / accumulator / ID-commitment / witness generation |

## Repository layout

```
include/        public headers
  poly.h          polynomial arithmetic over Fr (NTT-backed multiply, division, xGCD, ...)
  utils.h         characteristic-polynomial construction helpers
  zkacc.h         hidden accumulator + ZK (non-)membership / set-split / aggregated proofs
  ownpf.h         (non-)ownership proofs and (non-)membership witness generation
  txs.h           Entry / Transfer / Exit / Summary transactions
  test.h          benchmark entry points
src/            implementations matching the headers above (main() lives in test.cpp)
CMakeLists.txt  build configuration
```