#!/usr/bin/env python3
"""Cluster a TE library from an all-vs-all MMseqs2 search, the way CD-HIT-EST does, and
write CD-HIT-EST's output format.

CD_HIT_EST clusters large libraries with MMseqs2 instead of CD-HIT (see that module).
`mmseqs cluster` has no strand option and, on TE libraries, misses most pairs whose
consensi lie on opposite strands (rice: 26% of CD-HIT's opposite-strand pairs found, 54%
of same-strand ones). So the module runs a strand-aware search of the library against
itself, and this script does CD-HIT's greedy step on the hits:

  sequences are taken longest first (ties: input order); each joins the most identical
  earlier representative it aligns to with identity >= --identity over at least
  --coverage of its own length, on either strand; otherwise it becomes a representative.

Identity is MMseqs2's fident (identical columns / alignment length), the analogue of
cd-hit-est -G 0; coverage is of the member, the shorter sequence (-aS).

Everything downstream reads CD-HIT's files, so this script writes them:

  <out>          representatives, in input order (as cd-hit-est writes them)
  <out>.clstr    one block per cluster, as cd-hit-est -p 1 -d 0 writes it:
                   0   2535nt, >rep#LTR/Gypsy... *
                   1   1088nt, >member#LTR/Gypsy... at 1:1088:120:1207/+/91.20%
                 member start:end on the member, then on the representative, the
                 strand and the identity, from the alignment the member was assigned by.
                 DEDUP_REFINE reads the representative span.

Input: the clustered FASTA and the hit table (convertalis --format-output
query,target,fident,qstart,qend,qlen,tstart,tend,tlen). Each hit is read in both
directions, so a pair the search reported only one way still counts.

The hit table names each sequence by its record number in the FASTA, s1, s2, ...: the
module searches a copy with headers renamed that way, because MMseqs2 does not export
header words as given (its accession parser drops a 'consensus_' prefix and reads
'lcl|', 'ref|', ... as database prefixes), so >consensus_1 would come back as 1. The
original names are restored here, so the output carries them unchanged.
"""
import argparse
import sys


def read_fasta(path):
    seqs, order, name, buf = {}, [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(buf)
                name = line[1:].split()[0]
                order.append(name)
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = ''.join(buf)
    return seqs, order


def candidates(path, name, pos, length, identity, coverage):
    """member -> {rep: (ident, mstart, mend, rstart, rend, strand)}, the best alignment of
    each member on each sequence that comes before it in the greedy order. `name` maps the
    hit table's record numbers (s1, s2, ...) to the FASTA names."""
    best = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if f[0] == f[1]:
                continue
            fid = float(f[2])
            if fid < identity:
                continue
            if f[0] not in name or f[1] not in name:
                sys.exit(f'mmseqs_to_clstr: hit {f[0]} -> {f[1]} is not a record number (s1..s{len(name)}) of the FASTA')
            q, t = name[f[0]], name[f[1]]
            qs, qe, ts, te = int(f[3]), int(f[4]), int(f[6]), int(f[7])
            # MMseqs2 reports a minus-strand hit with one span reversed
            strand = '-' if (qs > qe) != (ts > te) else '+'
            # the later sequence in greedy order is the member, the earlier the candidate
            if pos[q] > pos[t]:
                m, r, ms, me, rs, re_ = q, t, qs, qe, ts, te
            else:
                m, r, ms, me, rs, re_ = t, q, ts, te, qs, qe
            if abs(me - ms) + 1 < coverage * length[m]:
                continue
            hit = (fid, min(ms, me), max(ms, me), min(rs, re_), max(rs, re_), strand)
            cur = best.setdefault(m, {}).get(r)
            if cur is None or hit[0] > cur[0]:
                best[m][r] = hit
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fasta', required=True, help='the library that was searched against itself')
    ap.add_argument('--hits', required=True,
                    help='mmseqs convertalis: query, target, fident, qstart, qend, qlen, tstart, tend, tlen; '
                         'query and target are FASTA record numbers, s1, s2, ...')
    ap.add_argument('--identity', type=float, default=0.8, help='minimum identity (fraction)')
    ap.add_argument('--coverage', type=float, default=0.8, help='minimum coverage of the member (fraction)')
    ap.add_argument('--out', required=True, help='representative FASTA; <out>.clstr is written too')
    args = ap.parse_args()

    seqs, order = read_fasta(args.fasta)
    if len(seqs) != len(order):
        sys.exit(f'mmseqs_to_clstr: {args.fasta} repeats a sequence name')
    rank = {n: i for i, n in enumerate(order)}
    length = {n: len(s) for n, s in seqs.items()}
    greedy = sorted(order, key=lambda n: (-length[n], rank[n]))
    pos = {n: i for i, n in enumerate(greedy)}
    name = {f's{i}': n for i, n in enumerate(order, 1)}
    cands = candidates(args.hits, name, pos, length, args.identity, args.coverage)

    members, aln = {}, {}                          # rep -> [member, ...]; member -> alignment
    for s in greedy:
        ok = [(hit[0], -rank[r], r) for r, hit in cands.get(s, {}).items() if r in members]
        if ok:
            rep = max(ok)[2]
            members[rep].append(s)
            aln[s] = cands[s][rep]
        else:
            members[s] = []

    reps = sorted(members, key=lambda r: rank[r])
    with open(args.out, 'w') as fa, open(args.out + '.clstr', 'w') as cl:
        for c, rep in enumerate(reps):
            seq = seqs[rep]
            fa.write(f'>{rep}\n')
            for i in range(0, len(seq), 60):
                fa.write(seq[i:i + 60] + '\n')
            cl.write(f'>Cluster {c}\n')
            cl.write(f'0\t{len(seq)}nt, >{rep}... *\n')
            # cd-hit lists members by length, longest first
            for k, mem in enumerate(sorted(members[rep], key=lambda m: (-length[m], rank[m])), 1):
                fid, ms, me, rs, re_, strand = aln[mem]
                cl.write(f'{k}\t{length[mem]}nt, >{mem}... at {ms}:{me}:{rs}:{re_}/{strand}/{fid * 100:.2f}%\n')
    print(f'mmseqs_to_clstr: {len(order)} sequences -> {len(reps)} clusters', file=sys.stderr)


if __name__ == '__main__':
    main()
