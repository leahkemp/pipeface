#!/usr/bin/env python3

"""
Compute chunk boundaries for parallel Jasmine merging.

Usage:
    python compute_chunks.py \
        --gaps gap_regions.sorted.bed.gz \
        --genome genome_sizes.txt \
        --positions variant_positions.bed \
        --n-chunks 5 \
        --min-gap 10000 \
        --output chunk_boundaries.bed

"""

import argparse
import gzip
import sys
from bisect import bisect_left

# minimum number of distinct variants in a chunk across all samples
MIN_VARIANTS = 50

def parse_args():
    p = argparse.ArgumentParser(description="Compute genomic chunk boundaries from gap regions.")
    p.add_argument("--gaps",      required=True, help="Gap regions BED (gzipped)")
    p.add_argument("--genome",    required=True, help="Chromosome sizes file (2-column TSV)")
    p.add_argument("--positions", required=True, help="Per-variant positions BED (CHROM\\tSTART\\t...)")
    p.add_argument("--n-chunks",  required=True, type=int, help="Number of chunks per chromosome")
    p.add_argument("--min-gap",   required=True, type=int, help="Minimum gap size in bp")
    p.add_argument("--output",    required=True, help="Output BED file")
    return p.parse_args()


def load_genome_sizes(path):
    sizes = {}
    with open(path) as f:
        for line in f:
            chrom, size = line.strip().split("\t")[:2]
            sizes[chrom] = int(size)
    return sizes


def load_gaps(path, min_gap):
    """Load gaps per chromosome, filtering by minimum size."""
    gaps = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            chrom, start, end = line.strip().split("\t")[:3]
            start, end = int(start), int(end)
            if end - start >= min_gap:
                gaps.setdefault(chrom, []).append((start, end))
    return gaps


def load_positions(path):
    """Load distinct variant start positions per chromosome"""
    positions = {}
    with open(path) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 2:
                continue
            chrom, start = cols[0], int(cols[1])
            positions.setdefault(chrom, set()).add(start)
    return {chrom: sorted(starts) for chrom, starts in positions.items()}


def find_nearest_gap(pos, gaps):
    """Find the gap midpoint nearest to pos."""
    best = None
    best_dist = float("inf")
    for start, end in gaps:
        mid = (start + end) // 2
        dist = abs(mid - pos)
        if dist < best_dist:
            best_dist = dist
            best = mid
    return best


def compute_boundaries(chrom, size, gaps, n_chunks):
    """
    Compute chunk boundaries for a chromosome by snapping
    evenly spaced target positions to the nearest valid gap.
    """
    if n_chunks == 1 or not gaps:
        return [(0, size)]

    targets = [int(size * i / n_chunks) for i in range(1, n_chunks)]
    boundaries = []
    prev = 0

    for target in targets:
        snap = find_nearest_gap(target, gaps)
        if snap is None or snap <= prev:
            # no valid gap near this target, skip this boundary
            continue
        boundaries.append(snap)
        prev = snap

    # build chunks from boundaries
    chunks = []
    starts = [0] + boundaries
    ends = boundaries + [size]
    for s, e in zip(starts, ends):
        chunks.append((s, e))

    return chunks


def merge_small_chunks(chunks, positions, min_variants):
    """
    Within a chromosome, merge adjacent chunks until each holds
    >= min_variants variants. If the trailing chunk falls short,
    fold it back into the previous chunk.
    """
    def count_in(start, end):
        return bisect_left(positions, end) - bisect_left(positions, start)

    if not chunks:
        return []

    merged = []
    cur_start, cur_end = chunks[0]
    cur_count = count_in(cur_start, cur_end)

    for nxt_start, nxt_end in chunks[1:]:
        if cur_count < min_variants:
            cur_end = nxt_end
            cur_count += count_in(nxt_start, nxt_end)
        else:
            merged.append((cur_start, cur_end, cur_count))
            cur_start, cur_end = nxt_start, nxt_end
            cur_count = count_in(cur_start, cur_end)

    if cur_count < min_variants and merged:
        # trailing chunk too small — fold back into the previously emitted chunk
        prev_start, _, prev_count = merged[-1]
        merged[-1] = (prev_start, cur_end, prev_count + cur_count)
    else:
        merged.append((cur_start, cur_end, cur_count))

    return merged


def pool_across_chromosomes(per_chrom_chunks, min_variants):
    """Pool small chunks across chromosomes until each holds
    >= min_variants variants. If the trailing chunk falls short,
    fold it back into the previous chunk.
    """

    groups = []
    pool_buffer = []
    pool_count = 0

    for chrom, start, end, count in per_chrom_chunks:
        if count >= min_variants:
            if pool_buffer and pool_count >= min_variants:
                groups.append(pool_buffer)
                pool_buffer = []
                pool_count = 0
                groups.append([(chrom, start, end)])
            elif pool_buffer:
                pool_buffer.append((chrom, start, end))
                groups.append(pool_buffer)
                pool_buffer = []
                pool_count = 0
            else:
                groups.append([(chrom, start, end)])
        else:
            pool_buffer.append((chrom, start, end))
            pool_count += count
            if pool_count >= min_variants:
                groups.append(pool_buffer)
                pool_buffer = []
                pool_count = 0

    if pool_buffer:
        if groups:
            groups[-1].extend(pool_buffer)
        else:
            groups.append(pool_buffer)

    result = []
    per_chrom_counter = {}
    pool_counter = 0
    for group in groups:
        chroms_in_group = {c for c, _, _ in group}
        if len(chroms_in_group) == 1:
            chrom = next(iter(chroms_in_group))
            per_chrom_counter[chrom] = per_chrom_counter.get(chrom, 0) + 1
            name = f"chunk_{chrom}_{per_chrom_counter[chrom]}"
        else:
            pool_counter += 1
            name = f"pool_{pool_counter}"
        for c, s, e in group:
            result.append((c, s, e, name))

    return result


def main():
    args = parse_args()

    genome    = load_genome_sizes(args.genome)
    gaps      = load_gaps(args.gaps, args.min_gap)
    positions = load_positions(args.positions)

    per_chrom_chunks = []
    for chrom, size in genome.items():
        chrom_gaps = gaps.get(chrom, [])
        if not chrom_gaps:
            print(f"[warn] no gaps found for {chrom} after filtering, writing as single chunk", file=sys.stderr)
        chunks = compute_boundaries(chrom, size, chrom_gaps, args.n_chunks)
        merged = merge_small_chunks(chunks, positions.get(chrom, []), MIN_VARIANTS)
        for start, end, count in merged:
            per_chrom_chunks.append((chrom, start, end, count))

    named_chunks = pool_across_chromosomes(per_chrom_chunks, MIN_VARIANTS)

    with open(args.output, "w") as out:
        for chrom, start, end, name in named_chunks:
            out.write(f"{chrom}\t{start}\t{end}\t{name}\n")

    print(f"[done] wrote {len(named_chunks)} chunks to {args.output}", file=sys.stderr)

if __name__ == "__main__":
    main()
