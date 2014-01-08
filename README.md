dnanexus / samtools
===================

This is a fork of [samtools](http://samtools.sourceforge.net/) maintained (without warranty) by [DNAnexus](https://www.dnanexus.com/). See [COPYING](https://github.com/dnanexus/samtools/blob/dnanexus/COPYING).

### Compiling

This repo incorporates [RocksDB](http://rocksdb.org/) as a submodule. Unfortunately RocksDB has several build dependencies; see [facebook/rocksdb/INSTALL.md](https://github.com/facebook/rocksdb/blob/master/INSTALL.md). Additionally, a development installation of [jemalloc](http://www.canonware.com/jemalloc/) is required (available in most Linux package managers).

Then,
```{bash}
git clone https://github.com/dnanexus/htslib.git
git clone https://github.com/dnanexus/samtools.git
make -C samtools
```

This will produce the executable `samtools/samtools`, which will probably assume several shared library dependencies (try `ldd samtools/samtools`).

### How to use rocksort

This section has specific instructions for using the `rocksort` subcommand; for a more general introduction, see our blog post.

#### Basic invocation

#### Estimating data size

#### Disabling background compaction

With limited hardware configurations, the background compaction performed by `samtools rocksort` may be counterproductive. You can effectively disable background compaction by setting the data size estimate absurdly large, e.g. `-s 999999T`. The result will usually still be modestly faster than `samtools sort`, except for datasets that fit entirely in the provisioned RAM.

#### Setting the scratch directory

By default, temporary files are written into the same directory as the final output BAM, similar to `samtools sort`. You can override this by setting the `TMPDIR` environment variable to something else. The scratch directory should be high-performance storage, and should have at least 3-4 times the expected size of the final BAM in free space.
