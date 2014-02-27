#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "bam.h"
#include <rocksdb/c.h>

#define SORT_KEY_SHUFFLE -1
#define SORT_KEY_POS 0
#define SORT_KEY_QNAME 1

unsigned int shuffling_seed;
/* simple 64-bit hashing for shuffling
   adapted from https://gist.github.com/badboy/6267743 and bamshuf.c */
static inline uint64_t hash64shift(uint64_t key)
{
  key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  return key;
}
static inline uint64_t hash64(const char *key, size_t len, uint64_t h) {
	h = hash64shift(h);
	while(len>0) h = (h << 5) - h + key[--len];
	return hash64shift(h);
}

/* copied from bam_sort.c */
static int strnum_cmp(const char *_a, const char *_b)
{
	const unsigned char *a = (const unsigned char*)_a, *b = (const unsigned char*)_b;
	const unsigned char *pa = a, *pb = b;
	while (*pa && *pb) {
		if (isdigit(*pa) && isdigit(*pb)) {
			while (*pa == '0') ++pa;
			while (*pb == '0') ++pb;
			while (isdigit(*pa) && isdigit(*pb) && *pa == *pb) ++pa, ++pb;
			if (isdigit(*pa) && isdigit(*pb)) {
				int i = 0;
				while (isdigit(pa[i]) && isdigit(pb[i])) ++i;
				return isdigit(pa[i])? 1 : isdigit(pb[i])? -1 : (int)*pa - (int)*pb;
			} else if (isdigit(*pa)) return 1;
			else if (isdigit(*pb)) return -1;
			else if (pa - a != pb - b) return pa - a < pb - b? 1 : -1;
		} else {
			if (*pa != *pb) return (int)*pa - (int)*pb;
			++pa; ++pb;
		}
	}
	return *pa? 1 : *pb? -1 : 0;
}

/* copied from bam_sort.c */
static int change_SO(bam_header_t *h, const char *so)
{
	char *p, *q, *beg = 0, *end = 0, *newtext;
	if (h->l_text > 3) {
		if (strncmp(h->text, "@HD", 3) == 0) {
			if ((p = strchr(h->text, '\n')) == 0) return -1;
			*p = '\0';
			if ((q = strstr(h->text, "\tSO:")) != 0) {
				*p = '\n'; // change back
				if (strncmp(q + 4, so, p - q - 4) != 0) {
					beg = q;
					for (q += 4; *q != '\n' && *q != '\t'; ++q);
					end = q;
				} else return 0; // no need to change
			} else beg = end = p, *p = '\n';
		}
	}
	if (beg == 0) { // no @HD
		h->l_text += strlen(so) + 15;
		newtext = malloc(h->l_text + 1);
		sprintf(newtext, "@HD\tVN:1.3\tSO:%s\n", so);
		strcat(newtext, h->text);
	} else { // has @HD but different or no SO
		h->l_text = (beg - h->text) + (4 + strlen(so)) + (h->text + h->l_text - end);
		newtext = malloc(h->l_text + 1);
		strncpy(newtext, h->text, beg - h->text);
		sprintf(newtext + (beg - h->text), "\tSO:%s", so);
		strcat(newtext, end);
	}
	free(h->text);
	h->text = newtext;
	return 0;
}

/* compare function for variable-length, little-endian sequence numbers
   we append to RocksDB key names to ensure uniqueness and sort stability */
static inline int compare_sequence_numbers(const char *seqnoa, size_t alen, const char *seqnob, size_t blen) {
	uint64_t ai = 0, bi = 0;
	if (alen == blen) {
		for(; alen; alen--) {
			ai <<= 8; bi <<= 8;
			ai |= (uint8_t) seqnoa[alen-1];
			bi |= (uint8_t) seqnob[alen-1];
		}
		if (ai < bi) {
			return -1;
		} else if (ai > bi) {
			return 1;
		}
		return 0;
	} else if (alen < blen) {
		return -1;
	}
	return 1;
}

/* RocksDB comparator for genome positions encoded as uint64_t's (as in
   bam1_lt in bam_sort.c) */
static int pos_rocksdb_comparator(void *ctx, const char *a, size_t alen, const char *b, size_t blen) {
	register uint64_t ai = *((uint64_t*) a);
	register uint64_t bi = *((uint64_t*) b);
	if (ai < bi) {
		return -1;
	} else if (ai > bi) {
		return 1;
	}
	/* rarely - break tie with sequence numbers */
	return compare_sequence_numbers(a+sizeof(uint64_t),alen-sizeof(uint64_t),
		                            b+sizeof(uint64_t),blen-sizeof(uint64_t));
}
static const char* pos_rocksdb_comparator_name(void *ctx) {
	return "samtools_rocksort_pos";
}

/* RocksDB comparator for qnames. The key encoding consists of the null-
   terminated qname string followed immediately by flag:uint32_t. */
static int qname_rocksdb_comparator(void *ctx, const char *a, size_t alen, const char *b, size_t blen) {
	register int t;
	uint32_t aflag, bflag;
	size_t qnamealen, qnameblen;
	if ((t = strnum_cmp(a,b)) != 0) {
		return t;
	}
	qnamealen = strlen(a);
	qnameblen = strlen(b);
	aflag = *(uint32_t*)(a + qnamealen + 1);
	bflag = *(uint32_t*)(b + qnameblen + 1);
	if (aflag < bflag) {
		return -1;
	} else if (aflag > bflag) {
		return 1;
	}
	/* rarely - break tie with sequence numbers */
	return compare_sequence_numbers(a+qnamealen+1+sizeof(uint32_t), alen-qnamealen-1-sizeof(uint32_t),
		                            b+qnameblen+1+sizeof(uint32_t), blen-qnameblen-1-sizeof(uint32_t));
}
static const char* qname_rocksdb_comparator_name(void *ctx) {
	return "samtools_rocksort_qname";
}

/* No-op RocksDB comparator destructor */
static void nop_rocksdb_comparator_destructor(void *c) {
}

/* given a bam1_t, formulate the appropriate key for insertion into RocksDB */
static size_t formulate_key(const int sort_key, bam1_t *b, uint64_t seqno, char **keybufp, size_t *keybufsz) {
	char *keybuf = *keybufp;
	size_t keylen = 0;
	assert(*keybufsz >= 2*sizeof(uint64_t));
	switch (sort_key) {
		case SORT_KEY_POS:
			*(uint64_t*)keybuf = (((uint64_t)b->core.tid)<<32|(b->core.pos+1)<<1|bam1_strand(b));
			keylen = sizeof(uint64_t);
			break;
		case SORT_KEY_QNAME:
			/* null-terminated qname followed by flags byte */
			keylen = b->core.l_qname + sizeof(uint32_t);
			if (keylen+sizeof(uint64_t) > *keybufsz) {
				*keybufsz = keylen+sizeof(uint64_t);
				kroundup32(*keybufsz);
				keybuf = realloc(keybuf, *keybufsz);
				if (!keybuf) return 0;
			}
			memcpy(keybuf, bam1_qname(b), b->core.l_qname);
			*(uint32_t*)(keybuf+b->core.l_qname) = b->core.flag;
			break;
		case SORT_KEY_SHUFFLE:
			/* hash the qname with salt from the stdlib PRNG.
			   NB the period of the PRNG alone would be much too small to
			      avoid a lot of collisions in any reasonably large BAM */
			*(uint64_t*)keybuf = hash64(bam1_qname(b), b->core.l_qname, rand_r(&shuffling_seed));
			keylen = sizeof(uint64_t);
			break;
		default: assert(0);
	}
	/* append the variable-length, little-endian sequence number to the key,
	   ensuring uniqueness & sort stability */
	while (seqno > 0) {
		*(uint8_t*)(keybuf+keylen) = (uint8_t) (seqno & 0xFF);
		seqno >>= 8;
		keylen++;
	}
	*keybufp = keybuf;
	return keylen;
}

/* Load contents of BAM fp into RocksDB rdb, keyed by leftmost genome position
   (encoded as a uint64_t) or by qname (if sort_key==SORT_KEY_QNAME) */
static int bam_to_rocksdb(bamFile fp, rocksdb_t *rdb, const int sort_key, unsigned long long *count, unsigned long long *actual_data_size) {
	int ret = 0;
	bam1_t *b = 0;
	size_t buflen = 1024;
	rocksdb_writeoptions_t *rdbwropts = 0;
	char *rdberr = 0;
	char *key = 0;
	size_t keybufsz = 0;
	size_t keylen = 0;

	key = calloc(2, sizeof(uint64_t));
	if (!key) {
		ret = -4;
		goto cleanup;
	}
	keybufsz = 2*sizeof(uint64_t);

	*count=0;
	*actual_data_size=0;

	rdbwropts = rocksdb_writeoptions_create();
	b = (bam1_t*) malloc(sizeof(bam1_t)+buflen);
	if (!rdbwropts || !b) {
		ret = -4;
		goto cleanup;
	}
	/* disable RocksDB write-ahead log since durability (to power loss, etc.)
	   is unnecessary in this context */
	rocksdb_writeoptions_disable_WAL(rdbwropts, 1);
	memset(b, 0, sizeof(bam1_t));
	
	/* for each BAM record */
	while ((ret = bam_read1(fp,b)) >= 0) {
		/* formulate key for insertion into RocksDB */
		if (!(keylen = formulate_key(sort_key, b, ++(*count), &key, &keybufsz))) {
			ret = -4;
			goto cleanup;
		}
		
		/* Copy b->data into the buffer directly following the bam1_t at b,
		   so that we can then stick the whole thing into RocksDB. It would
		   be nice to avoid this copying, but bam_read1 may call realloc on
		   b->data so unfortunately we can't just put it anywhere we want. */
		if (buflen < b->l_data) {
			// expand the buffer
			buflen = b->l_data;
			kroundup32(buflen);
			b = (bam1_t*) realloc(b, sizeof(bam1_t)+buflen);
			if (!b) {
				ret = -4;
				goto cleanup;
			}
		}
		memcpy(b+1, b->data, b->l_data);

		/* insert this key & value into RocksDB */
		rocksdb_put(rdb, rdbwropts, key, keylen, (char*) b, sizeof(bam1_t) + b->l_data, &rdberr);
		if (rdberr) {
			fprintf(stderr, "[bam_rocksort_core] %s\n", rdberr);
			ret = -1;
			goto cleanup;
		}
		*actual_data_size += keylen + sizeof(bam1_t) + b->l_data;
	}
	if (ret != -1)
		fprintf(stderr, "[bam_rocksort_core] truncated file. Continue anyway.\n");
	ret = 0;

cleanup:
	if (rdbwropts) rocksdb_writeoptions_destroy(rdbwropts);
	if (b) {
		if (b->data) free(b->data);
		free(b);
	}
	if (keybufsz) {
		free(key);
	}
	if (rdberr) {
		free(rdberr);
	}
	return ret;
}

/* Traverse the RocksDB and output a BAM file */
static int rocksdb_to_bam(rocksdb_t *rdb, const bam_header_t *header, const char *fnout, const int n_threads, const int level, unsigned long reopen, unsigned long long *count) {
	int ret = 0;
	rocksdb_readoptions_t *rdbrdopts = 0;
	rocksdb_iterator_t *rdbiter = 0, *rdbiter2 = 0;
	char *rdberr = 0;
	bamFile fp = 0;
	bam1_t *b;
	size_t vsz = 0;
	char mode[8];
	const char *reopen_key = 0, *reopened_key = 0;
	size_t reopen_key_len = 0, reopened_key_len = 0;

	*count = 0;

	if (!(rdbrdopts = rocksdb_readoptions_create())) {
		ret = -4;
		goto cleanup;
	}
	rocksdb_readoptions_set_verify_checksums(rdbrdopts, 0);
	rocksdb_readoptions_set_prefetch(rdbrdopts, 1);
	if (!(rdbiter = rocksdb_create_iterator(rdb, rdbrdopts))) {
		ret = -4;
		goto cleanup;
	}

	strcpy(mode, "w");
	if (level >= 0) sprintf(mode + 1, "%d", level < 9? level : 9);
	if (!(fp = strcmp(fnout, "-") ? bam_open(fnout, mode) : bam_dopen(fileno(stdout), mode))) {
		fprintf(stderr, "[bam_rocksort_core] fail to create the output file %s\n",fnout);
		ret = -1;
		goto cleanup;
	}
	bam_header_write(fp, header);
	if (n_threads > 1) bgzf_mt(fp, n_threads, 256);

	/* For each record in RocksDB */
	rocksdb_iter_seek_to_first(rdbiter);
	while (rocksdb_iter_valid(rdbiter)) {
        /* Extract bam1_t, knowing that the RocksDB value consists of a bam1_t
		   with a garbage data pointer, immediately followed by the data */
		b = (bam1_t*) rocksdb_iter_value(rdbiter, &vsz);
		b->data = (uint8_t*)(b+1);
		b->m_data = vsz - sizeof(bam1_t);
		/* Write to output BAM */
		bam_write1(fp, b);
		(*count)++;

		/* if requested, reopen the iterator at the current key, to start
		   taking advantage of any background compactions that've finished
		   in the meantime. */
		if (reopen > 0 && (*count) == reopen) {
			reopen_key = rocksdb_iter_key(rdbiter, &reopen_key_len);
			if (!(rdbiter2 = rocksdb_create_iterator(rdb, rdbrdopts))) {
				ret = -4;
				goto cleanup;
			}
			rocksdb_iter_seek(rdbiter2, reopen_key, reopen_key_len);
			if (!rocksdb_iter_valid(rdbiter2)) {
				ret = -1;
				break;
			}
			/* sanity check */
			reopened_key = rocksdb_iter_key(rdbiter2, &reopened_key_len);
			if (reopen_key_len != reopened_key_len || memcmp(reopen_key, reopened_key, reopen_key_len)) {
				fprintf(stderr, "[bam_rocksort_core] BUG: RocksDB iterator reopened at different key!\n");
				ret = -1;
				break;
			}
			rocksdb_iter_destroy(rdbiter);
			rdbiter = rdbiter2;
			rdbiter2 = 0;
		}

		rocksdb_iter_next(rdbiter);
	}

	rocksdb_iter_get_error(rdbiter, &rdberr);
	if (rdberr) {
		fprintf(stderr, "[bam_rocksort_core] %s\n", rdberr);
		ret = -1;
		goto cleanup;
	}

	if (fp->errcode) {
		fprintf(stderr, "[bam_rocksort_core] error writing final BAM: %x\n", fp->errcode);
		ret = -1;
		goto cleanup;
	}

cleanup:
	if (fp) bam_close(fp);
	if (rdbiter) rocksdb_iter_destroy(rdbiter);
	if (rdbiter2) rocksdb_iter_destroy(rdbiter2);
	if (rdberr) free(rdberr);
	if (rdbrdopts) rocksdb_readoptions_destroy(rdbrdopts);

	return ret;
}

/* If TMPDIR environment variable is defined, create a temp directory there.
   Otherwise create a directory name based on prefix.*/
char* choose_rocksdb_path(const char *prefix) {
	char *tmpdir = getenv("TMPDIR");
	char *ans = 0;
	const char *template_const = "/samtools_rocksort_XXXXXX";
	const char *suffix = ".rocksort";

	if (tmpdir && *tmpdir) {
		ans = malloc(strlen(tmpdir)+strlen(template_const)+1);
		strcpy(ans, tmpdir);
		strcat(ans, template_const);
		if (!mkdtemp(ans)) {
			free(ans);
			return NULL;
		}
	} else {
		ans = malloc(strlen(prefix)+strlen(suffix)+1);
		strcpy(ans, prefix);
		strcat(ans, suffix);
	}

	return ans;
}

/*!
  @abstract Sort an unsorted BAM file based on the chromosome order
  and the leftmost position of an alignment

  @param  sort_key one of SORT_KEY_POS, SORT_KEY_QNAME, or SORT_KEY_SHUFFLE
  @param  fn       name of the file to be sorted
  @param  prefix   prefix of the temporary files (prefix.NNNN.bam are written)
  @param  fnout    name of the final output file to be written
  @param  max_mem  approxiate maximum memory (very inaccurate)
  @return 0 for successful sorting, negative on errors

  @discussion It may create multiple temporary subalignment files
  and then merge them by calling bam_merge_core(). This function is
  NOT thread safe.
 */
int bam_rocksort_core_ext(const int sort_key, const char *fn, const char *prefix, const char *fnout, size_t max_mem_per_thread, size_t data_size_hint, int n_threads, const int level, const int keep_db) {
	int ret = 0;
	bamFile fp = 0;
	bam_header_t *header = 0;
	rocksdb_env_t *rdbenv = 0;
	rocksdb_t *rdb = 0;
	rocksdb_options_t *rdbopts = 0;
	rocksdb_comparator_t *rdbcomp = 0;
	rocksdb_universal_compaction_options_t *rdbucopts = 0;
	rocksdb_cache_t *rdbcache = 0;
	char *rdberr = 0;
	char *rdbpath = 0;
	unsigned long long count1 = 0, count2 = 0, actual_data_size = 0;
	size_t bytes_per_file = 0;
	unsigned int max_background_compactions = 0, max_background_flushes = 0;
	unsigned int files_per_compaction = 0;
	const size_t MB = ((size_t)1)<<20;
	const size_t GB = ((size_t)1)<<30;
	size_t max_mem;

	if (n_threads < 1) n_threads=1;
	max_mem = max_mem_per_thread * n_threads;
	if (max_mem < 512*MB) max_mem = 512*MB;

	if (!(fp = strcmp(fn, "-")? bam_open(fn, "r") : bam_dopen(fileno(stdin), "r"))) {
		fprintf(stderr, "[bam_rocksort_core] fail to open file %s\n", fn);
		ret = -1;
		goto cleanup;
	}
	if (!(header = bam_header_read(fp))) {
		fprintf(stderr, "[bam_rocksort_core] fail to parse file %s\n", fn);
		ret = -1;
		goto cleanup;
	}
	switch (sort_key) {
		case SORT_KEY_POS:
			change_SO(header, "coordinate");
			break;
		case SORT_KEY_QNAME:
			change_SO(header, "queryname");
			break;
		case SORT_KEY_SHUFFLE:
			change_SO(header, "unsorted");
			break;
		default:
			assert(0);
	}

	/* Configure & open RocksDB */
	if (sort_key == SORT_KEY_QNAME) {
		rdbcomp = rocksdb_comparator_create(0, nop_rocksdb_comparator_destructor, qname_rocksdb_comparator, qname_rocksdb_comparator_name);
	} else {
		/* the position comparator is used for both SORT_KEY_POS and SORT_KEY_SHUFFLE */
		rdbcomp = rocksdb_comparator_create(0, nop_rocksdb_comparator_destructor, pos_rocksdb_comparator, pos_rocksdb_comparator_name);
	}
	rdbenv = rocksdb_create_default_env();
	rdbucopts = rocksdb_universal_compaction_options_create();
	rdbopts = rocksdb_options_create();
	rdbcache = rocksdb_cache_create_lru(512*MB);
	if (!rdbcomp || !rdbopts || !rdbenv || !rdbucopts || !rdbcache) {
		ret = -4;
		goto cleanup;
	}

	/* ROCKSDB CONFIGURATION: ESSENTIALS */
	rocksdb_options_set_create_if_missing(rdbopts, 1);
	rocksdb_options_set_error_if_exists(rdbopts, 1);
	rocksdb_options_set_comparator(rdbopts, rdbcomp);
	rocksdb_options_set_num_levels(rdbopts, 1);
	rocksdb_options_set_compression(rdbopts, rocksdb_snappy_compression);
	rocksdb_options_set_cache(rdbopts, rdbcache);
	rocksdb_options_set_max_open_files(rdbopts, 256);
	/* turn off some durability and throttling features (unnecessary here) */
	rocksdb_options_set_disable_data_sync(rdbopts, 1);
	rocksdb_options_set_paranoid_checks(rdbopts, 0);
	rocksdb_options_set_disable_seek_compaction(rdbopts, 1);
	rocksdb_options_set_level0_slowdown_writes_trigger(rdbopts, GB);
	rocksdb_options_set_level0_stop_writes_trigger(rdbopts, GB);

	/* CONCURRENCY */
	/* up to two concurrent memtable flushes (see below) */
	rocksdb_env_set_background_threads(rdbenv, n_threads);
	max_background_flushes = n_threads > 1 ? 2 : 1;
	rocksdb_env_set_high_priority_background_threads(rdbenv, max_background_flushes);
	rocksdb_options_set_max_background_flushes(rdbopts, max_background_flushes);
	max_background_compactions = n_threads-1;
	if (max_background_compactions > 3) {
		/* set a ceiling on concurrent background compactions, since each one
		   seems to use up substantial memory */
		max_background_compactions = 3;
	}
	if (data_size_hint > 0 && max_background_compactions > 0) {
		rocksdb_options_set_max_background_compactions(rdbopts, max_background_compactions);
	} else {
		rocksdb_options_set_disable_auto_compactions(rdbopts, 1);
	}
	rocksdb_options_set_env(rdbopts, rdbenv);
	/* note: for n_threads <= 4, we may have provisioned n_threads+1
	   background threads. in practice however, each is likely to spend a
	   lot of time either idle or I/O-bound. */

	/* TUNE IN-MEMORY SORTING & COMPRESSION */

	/* use vectors instead of skip lists for memtables - better for bulk
	   loading */
	rocksdb_options_set_memtable_vector_rep(rdbopts);

	/* Keep up to 3 memtables in memory (one loading, two flushing). If
	   background compaction is enabled, set aside some memory for them by
	   making each memtable 1/4 instead of 1/3 of the total provisioned.
       With two flushing memtables we hope that RocksDB takes no more than
       2X as long to sort, Snappy-compress and write each buffer as it
       does for us to initially read, parse and insert it. */
	bytes_per_file = max_mem / (data_size_hint > 0 ? 4 : 3);
	rocksdb_options_set_write_buffer_size(rdbopts, bytes_per_file);
	rocksdb_options_set_max_write_buffer_number(rdbopts, 3);
	/* compress in 2MB blocks */
	rocksdb_options_set_block_size(rdbopts, 2*MB);

	/* TUNE ON-DISK COMPACTION  */
	/* use 'universal' compaction which has lower write amplification and
	   higher read amplification than level compaction - this is the right
	   tradeoff given that we're ultimately just going to do one sequential
	   scan over the DB and then delete it. */
	if (data_size_hint > 0) {
		rocksdb_options_set_compaction_style(rdbopts, rocksdb_universal_compaction);

		/* if there will be T total buffers, merge them in batches of sqrt(T) */
		files_per_compaction = ((unsigned int)sqrt(((float)data_size_hint)/bytes_per_file)) + 1;
		rocksdb_options_set_level0_file_num_compaction_trigger(rdbopts, 2*files_per_compaction);
		rocksdb_universal_compaction_options_set_min_merge_width(rdbucopts, files_per_compaction);
		rocksdb_universal_compaction_options_set_max_merge_width(rdbucopts, 2*files_per_compaction);
		rocksdb_universal_compaction_options_set_stop_style(rdbucopts, rocksdb_similar_size_compaction_stop_style);
		rocksdb_universal_compaction_options_set_size_ratio(rdbucopts, 50);
		/* disable compaction for size amplification - useless since we won't be
		   overwriting or deleting anything */
		rocksdb_universal_compaction_options_set_max_size_amplification_percent(rdbucopts, 1<<30);
		rocksdb_options_set_universal_compaction_options(rdbopts, rdbucopts);
	}
	/* END OF ROCKSDB CONFIGURATION */

	rdbpath = choose_rocksdb_path(prefix);
	if (!rdbpath) {
		fprintf(stderr, "[bam_rocksort] Failed creating temporary directory; check given prefix and TMPDIR environment variable.\n");
		ret = -1;
		goto cleanup;
	}
	if (!(rdb = rocksdb_open(rdbopts, rdbpath, &rdberr))) {
		fprintf(stderr, "[bam_rocksort] %s\n", rdberr ? rdberr : "failed creating RocksDB");
		ret = -4;
		goto cleanup;
	}

	/* Load input BAM into RocksDB */ 
	fprintf(stderr, "[bam_rocksort_core] Sorting in %s/ (you can change this with the TMPDIR environment variable)...\n", rdbpath);
	if ((ret = bam_to_rocksdb(fp, rdb, sort_key, &count1, &actual_data_size)) != 0) {
		goto cleanup;
	}

	/* Provide feedback on data_size_hint */
	if (data_size_hint > 0) {
		if (((float)actual_data_size)/data_size_hint > 1.2) {
			fprintf(stderr, "[bam_rocksort_core] WARNING: actual uncompressed data size (%llu bytes) was well above the assumption (%llu bytes); consider providing an accurate size hint (-s) to optimize sort performance\n", actual_data_size, (unsigned long long) data_size_hint);
		} else if (((float)actual_data_size)/data_size_hint < 0.8) {
			fprintf(stderr, "[bam_rocksort_core] WARNING: actual uncompressed data size (%llu bytes) was well below the assumption (%llu bytes); consider providing an accurate size hint (-s) to optimize sort performance\n", actual_data_size, (unsigned long long) data_size_hint);
		}  else {
			fprintf(stderr, "[bam_rocksort_core] actual uncompressed data size (%llu bytes) was pretty close to the assumption (%llu bytes)\n", actual_data_size, (unsigned long long) data_size_hint);
		}
	}

	/* Export sorted BAM from RocksDB */
	fprintf(stderr, "[bam_rocksort_core] Writing %llu records to %s...\n", count1, fnout);
	if ((ret = rocksdb_to_bam(rdb, header, fnout, n_threads, level, count1/3, &count2)) != 0) {
		goto cleanup;
	}

	/* sanity check */
	if (count1 != count2) {
		fprintf(stderr, "[bam_rocksort_core] BUG: read %llu, wrote %llu records!\n", count1, count2);
		ret = -1;
		goto cleanup;
	}

	fprintf(stderr, "[bam_rocksort_core] OK\n");

cleanup:
	if(fp) bam_close(fp);
	if(header) bam_header_destroy(header);
	if(rdb) {
		rocksdb_close(rdb);
		if (!keep_db) rocksdb_destroy_db(rdbopts, rdbpath, &rdberr);
	}
	if(rdbcomp) rocksdb_comparator_destroy(rdbcomp);
	if(rdbopts) rocksdb_options_destroy(rdbopts);
	if(rdbenv) rocksdb_env_destroy(rdbenv);
	if(rdbucopts) rocksdb_universal_compaction_options_destroy(rdbucopts);
	if(rdbcache) rocksdb_cache_destroy(rdbcache);
	if(rdberr) free(rdberr);
	if(rdbpath) free(rdbpath);

	if (ret != 0) {
		fprintf(stderr, "[bam_rocksort_core] Exit status %d\n", ret);
	}
	return ret;
}

int bam_rocksort_core(int sort_key, const char *fn, const char *prefix, size_t max_mem_per_thread, size_t data_size_hint)
{
	int ret;
	char *fnout = calloc(strlen(prefix) + 4 + 1, 1);
	sprintf(fnout, "%s.bam", prefix);
	ret = bam_rocksort_core_ext(sort_key, fn, prefix, fnout, max_mem_per_thread, data_size_hint, 0, -1, 0);
	free(fnout);
	return ret;
}

int bam_rocksort(int argc, char *argv[])
{
	size_t max_mem = 768<<20; // 768MB
	size_t size_hint = 0;
	int c, sort_key = SORT_KEY_POS, is_stdout = 0, ret = 0,
	    n_threads = 0, level = -1, full_path = 0, keep_db = 0;
	char *fnout;
	shuffling_seed = 0;
	while ((c = getopt(argc, argv, "fnoks:m:@:l:u:")) >= 0) {
		switch (c) {
		case 'f': full_path = 1; break;
		case 'o': is_stdout = 1; break;
		case 'n': sort_key = SORT_KEY_QNAME; break;
		case 'k': keep_db=1; break;
		case 'm': {
				char *q;
				max_mem = strtol(optarg, &q, 0);
				if (*q == 'k' || *q == 'K') max_mem <<= 10;
				else if (*q == 'm' || *q == 'M') max_mem <<= 20;
				else if (*q == 'g' || *q == 'G') max_mem <<= 30;
				else if (*q == 't' || *q == 'T') max_mem <<= 40;
				break;
			}
		case 's': {
				char *q2;
				size_hint = strtol(optarg, &q2, 0);
				if (*q2 == 'k' || *q2 == 'K') size_hint <<= 10;
				else if (*q2 == 'm' || *q2 == 'M') size_hint <<= 20;
				else if (*q2 == 'g' || *q2 == 'G') size_hint <<= 30;
				else if (*q2 == 't' || *q2 == 'T') size_hint <<= 40;
				break;
			}
		case '@': n_threads = atoi(optarg); break;
		case 'l': level = atoi(optarg); break;
		case 'u':
			shuffling_seed = (unsigned int) atoi(optarg);
			sort_key = SORT_KEY_SHUFFLE;
			break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   samtools rocksort [options] <in.bam|-> <out.prefix>\n\n");
		fprintf(stderr, "    -@ INT    number of sorting and compression threads [1]\n");
		fprintf(stderr, "    -m INT    max memory per thread; suffix K/M/G/T recognized [768M]\n");
		fprintf(stderr, "    -s INT    plan background compactions assuming this uncompressed\n");
		fprintf(stderr, "              BAM data size; suffix K/M/G/T recognized [off]\n\n");
		fprintf(stderr, "    -f        use <out.prefix> as full file name instead of prefix\n");
		fprintf(stderr, "    -o        final output to stdout\n");
		fprintf(stderr, "    -l INT    compression level, from 0 to 9 [-1]\n");
		fprintf(stderr, "    -k        keep RocksDB instead of deleting it when done\n\n");
		fprintf(stderr, "    -n        sort by read name\n");
		fprintf(stderr, "    -u INT    unsort: shuffle the BAM using given random seed\n");
		fprintf(stderr, "\n");
		return 1;
	}

	if (is_stdout) fnout = strdup("-");
	else {
		fnout = calloc(strlen(argv[optind+1]) + 4 + 1, 1);
		sprintf(fnout, "%s%s", argv[optind+1], full_path? "" : ".bam");
	}

	if (bam_rocksort_core_ext(sort_key, argv[optind], argv[optind+1], fnout, max_mem, size_hint, n_threads, level, keep_db) < 0) ret = 1;
	free(fnout);
	return ret;
}
