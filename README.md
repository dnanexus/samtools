[dnanexus / samtools](https://github.com/dnanexus/samtools)
===================

This is a fork of [samtools](http://samtools.sourceforge.net/) maintained (without warranty) by [DNAnexus](https://www.dnanexus.com/). See [COPYING](https://github.com/dnanexus/samtools/blob/dnanexus/COPYING).

### Compiling

This version of samtools incorporates [RocksDB](http://rocksdb.org/) as a [git submodule](http://git-scm.com/docs/git-submodule). RocksDB has several build dependencies; see [facebook/rocksdb/INSTALL.md](https://github.com/facebook/rocksdb/blob/master/INSTALL.md). Additionally, a development installation of [jemalloc](http://www.canonware.com/jemalloc/) is required (available in most Linux package managers).

Then,
```{bash}
git clone https://github.com/dnanexus/htslib.git
git clone https://github.com/dnanexus/samtools.git
make -C samtools
```

This will produce the executable `samtools/samtools`, which will probably assume several shared library dependencies (try `ldd samtools/samtools`).

### Using rocksort

This section has specific instructions for using the `rocksort` subcommand; for a more general introduction, see our blog post. The basic invocation of `samtools rocksort` is similar to `samtools sort`:

```
Usage:   samtools rocksort [options] <in.bam> <out.prefix>

Options: -@ INT    number of sorting and compression threads [1]
         -m INT    max memory per thread; suffix K/M/G/T recognized [768M]
```

We recommend using at least four threads, with compute cores to run them. The memory requirements should be similar to `samtools sort` with the same `-@` and `-m` settings. 

To enable background compactions, supply the option `-s`.

```
         -s INT    plan background compactions assuming this uncompressed
                   BAM data size; suffix K/M/G/T recognized [off]
```

See below for guidance on setting the value of this option.

Additional options:

```
         -o        final output to stdout
         -f        use <out.prefix> as full file name instead of prefix
         -n        sort by read name
         -k        keep RocksDB instead of deleting it when done
         -l INT    compression level, from 0 to 9 [-1]
         -u INT    unsort: shuffle the BAM using given random seed
```

Like `samtools sort`, the input BAM can be streamed through standard input by supplying `-` instead of `<in.bam>`. In fact, if the input is a compressed BAM file, it's usually faster to decompress with [pigz](http://zlib.net/pigz/) like this:

```
cat <in.bam> | pigz -dc | samtools rocksort [options] - <out.prefix>
```

#### Estimating data size (`-s`) for background compactions

**Note:** background compactions are generally useful only when the dataset to be sorted is *many* times larger than provisioned memory. See our blog post introducing rocksort for more details.

To plan efficient background compactions, rocksort needs a rough estimate of the total *uncompressed* size of the BAM data. A rough rule of thumb is to quadruple the expected size of the final BAM file. For example, if you expect to produce a 125 GiB final BAM (roughly a deep human WGS), a size estimate of 500 GiB would work pretty well. Note that sorted BAMs are substantially smaller than unsorted BAMs, since they're more compressible.

Alternatively, here's a bottom-up formula for the size of a single BAM alignment block ([source](http://genome.sph.umich.edu/wiki/SAM)), which you can multiply by the expected number of read alignments:

```
Block Size = 8*4 + ReadNameLength(including null) + CigarLength*4 + (ReadLength+1)/2 + ReadLength + TagLength
```

The estimate needs not be fantastically accurate; +/- 20% or so is fine. If in doubt, overestimate. Once all the data are loaded, `samtools rocksort` will log some feedback about the accuracy of the hint to standard error. 


#### Setting the scratch directory

By default, temporary files are written into the same directory as the final output BAM, similar to `samtools sort`. You can override this by setting the `TMPDIR` environment variable to something else. The scratch directory should be high-performance storage, and should have at least twice the expected size of the final BAM (4-5X if using background compactions) in free space.
