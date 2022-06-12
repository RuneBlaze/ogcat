ogcat
=================

Fast and scalable phylogenomic utilities :cat:.

## Installation

### Prebuilt binaries

See releases. The `musl` binary for Linux should be the most compatible, at the possible
expense of some speed.

(I keep a separate "better" copy on the computing cluster at UIUC, so let me know if you need it
on the cluster.)

### Compiling from scratch

Install the [Rust toolchain](https://www.rust-lang.org/tools/install) and then compile the binary:

```shell
RUSTFLAGS="-C target-feature=+avx" cargo build --release
```

## Manual

See [wiki](https://github.com/RuneBlaze/ogcat/wiki).

## Acknowledgements

`ogcat` internally uses code translated from [TreeSwift](https://github.com/niemasd/TreeSwift), and its sum-of-pairs calculation is based on that from [FastSP](https://github.com/smirarab/FastSP).

`ogcat` stands for "orange cat" and does not relate to `cat` the unix utility.