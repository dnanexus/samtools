#!/bin/bash -e

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ] ; do SOURCE="$(readlink "$SOURCE")"; done
srcdir="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

htslibdir="${srcdir}/../../../htslib"
samtoolsdir="${srcdir}/../.."
samtools="${samtoolsdir}/samtools"

if [ -z "$TMPDIR" ]; then
	export TMPDIR="/tmp"
fi

# Fetch a test BAM file
if [ -z "$TEST_BAM" ]; then
	TEST_BAM="${TMPDIR}/samtools_rocksort_test.bam"
	if ! [ -f "$TEST_BAM" ]; then
		echo "Fetching a test BAM file (25MB) from the internet..."
		wget -O "$TEST_BAM" "https://test.galaxyproject.org/library_common/download_dataset_from_folder?library_id=04c38ea196069a5d&cntrller=library&use_panels=False&id=4a29c5f588caf7fb"
	fi
fi

function sam_unique_lines {
	echo $($samtools view $1 | sort | uniq | wc -l)
}

# count alignments
expected_sam_count=$(sam_unique_lines $TEST_BAM)
echo "SAM unique lines: $expected_sam_count"

# Shuffle it
testdir=$(mktemp -d --tmpdir samtools_rocksort_test.XXXXXX)
echo "Temporary test directory: $testdir"
cd "$testdir"
shuffled_bam="${testdir}/samtools_rocksort_test_shuffled.bam"
echo "Shuffling $TEST_BAM..."
$samtools rocksort -u 0 -l 1 $TEST_BAM "${testdir}/samtools_rocksort_test_shuffled"

# Compile and test bamsorted.c (FIXME ugly)
cc -o bamsorted -O2 ${srcdir}/bamsorted.c -I"$htslibdir" -I"$samtoolsdir" "$htslibdir/libhts.a" "$samtoolsdir/libbam.a" "$htslibdir/faidx.o" -lz -lpthread
echo -n "Positive control: "
./bamsorted "$TEST_BAM"
echo -n "Negative control: "
code=1
./bamsorted "$shuffled_bam" || code=$?
if [ "$code" -eq "0" ]; then
	echo "should have failed: bamsorted $shuffled_bam"
	exit 1
fi

echo "Running samtools rocksort..."
time $samtools rocksort -l 1 -@ 2 -m 16M "$shuffled_bam" "${testdir}/samtools_rocksort_test_sorted"
sorted_bam="${testdir}/samtools_rocksort_test_sorted.bam"

echo -n "Verifying product: "
result_sam_count=$(sam_unique_lines "$sorted_bam")
if [ "$result_sam_count" != "$expected_sam_count" ]; then
	echo "Mismatch of SAM unique lines: expected $expected_sam_count, got $result_sam_count"
	exit 1
fi
./bamsorted "$sorted_bam"

echo "Sort by read name..."
time $samtools rocksort -n -l 1 -@ 2 -m 16M "$shuffled_bam" "${testdir}/samtools_rocksort_test_sorted"
sorted_bam="${testdir}/samtools_rocksort_test_sorted_by_name.bam"

# Clean up
rm -rf $testdir
echo "Test passed. Cleaned up $testdir"
