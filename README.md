## libsais suffix array driver

### building

You must have a C++17 compliant compiler, meson, CMake and ninja installed.  From the top-level directory issue:

```bash
$ meson setup builddir
```

then

```bash
$ cd builddir
$ ninja
```

The resulting executable should be called `driver`.

### running 

The driver takes in the input file using the `-f` option, the output location with the `-o` option
the input type must be provided (options are `dna`, `integer` and `text`, though currently the only
difference between `dna` and `text` is that `dna` will attempt to parse a `FASTA` file with one record
while `text` will attempt to read a single string from the fiile. The `-t` options sets the number of 
threads. With 1 thread, the driver will use the serial algorithm and with > 1 it will use the openmp variants.

* The format of `dna` input is a FASTA file; only the first record will be read.
* The format of `text` input is a file containing just a single string.
* The format of `integer` input is a binary file starting with 2 u64 integers follwed by an array
  * The first u64 integer `N` is the length of the input
  * The second u64 `M` is the maximum character in the input
  * If either `N` or `M` is >= `i32::MAX`, then the following array is an array of `N` `u64` integers,
  otherwise it is an array of `N` `u32` integers.

The algorithm will be dispatched based on the input size (and in the integer case, the maximum token).

```bash
libsais driver


./builddir/driver [OPTIONS]


OPTIONS:
  -h,     --help              Print this help message and exit
  -f,     --file TEXT         input filename
  -o,     --output TEXT REQUIRED
                              output filename
          --input-type ENUM:value in {dna->,integer->,text->} OR {,,} REQUIRED
                              input type
  -t,     --threads UINT      number of threads

```
