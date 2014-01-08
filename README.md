[dnanexus / samtools](https://github.com/dnanexus/samtools)
===================

This is a fork of [samtools](http://samtools.sourceforge.net/) maintained (without warranty) by [DNAnexus](https://www.dnanexus.com/). See [COPYING](https://github.com/dnanexus/samtools/blob/dnanexus/COPYING).

### Compiling

This version of samtools incorporates [RocksDB](http://rocksdb.org/) as a [git submodule](http://git-scm.com/docs/git-submodule). Unfortunately RocksDB has several build dependencies; see [facebook/rocksdb/INSTALL.md](https://github.com/facebook/rocksdb/blob/master/INSTALL.md). Additionally, a development installation of [jemalloc](http://www.canonware.com/jemalloc/) is required (available in most Linux package managers).

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
         -s INT    hint as to total uncompressed BAM data size; suffix K/M/G/T recognized [512G]
```

We recommend using at least four threads, with compute cores to run them. The memory requirements should be similar to `samtools sort` with the same `-@` and `-m` settings. See below for guidance on setting `-s`.

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

#### Estimating data size (`-s`)

To plan background compactions efficiently, rocksort needs a rough estimate of the total *uncompressed* size of the BAM data. The default assumption is 512GiB, which corresponds to a deep human WGS, 125GiB compressed and sorted final product BAM as of this writing (early 2014). If your data is much smaller or larger than this, it will probably be beneficial to adjust the hint. The 4:1 ratio is a reasonable assumption, or alternatively, here's a bottom-up formula for the size of a single BAM alignment block ([source](http://genome.sph.umich.edu/wiki/SAM)), which you can multiply by the expected number of read alignments:

```
Block Size = 8*4 + ReadNameLength(including null) + CigarLength*4 + (ReadLength+1)/2 + ReadLength + TagLength
```

Once all the data are loaded, `samtools rocksort` will log some feedback about the accuracy of the hint to standard error. The estimate need not be fantastically accurate; +/- 20% or so is fine.

#### Disabling background compaction

With limited hardware configurations, the background compaction performed by `samtools rocksort` may be counterproductive. You can effectively disable background compaction by setting the data size estimate absurdly large, e.g. `-s 999999T`. The result will usually still be modestly faster than `samtools sort`, except for datasets that fit entirely in the provisioned RAM.

#### Setting the scratch directory

By default, temporary files are written into the same directory as the final output BAM, similar to `samtools sort`. You can override this by setting the `TMPDIR` environment variable to something else. The scratch directory should be high-performance storage, and should have at least the uncompressed data size (4-5X the expected size of the final BAM) in free space.
