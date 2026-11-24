"""bin/mmseqs_to_clstr.py: CD-HIT's greedy rule on an all-vs-all MMseqs2 hit table, on
both strands, written in CD-HIT's .clstr format."""
import random
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from conftest import BIN, headers

MODULE = Path(__file__).resolve().parents[2] / 'modules' / 'local' / 'cd_hit_est' / 'main.nf'
LONG, MID, SHORT, OTHER = 'A' * 1000, 'C' * 900, 'G' * 500, 'T' * 950
# the module searches record numbers (s1, s2, ...), not names: FASTA order below
ID = {'short': 's1', 'long': 's2', 'mid': 's3', 'other': 's4'}


def setup(tmp_path, hits):
    (tmp_path / 'lib.fa').write_text(
        f'>short#DNA\n{SHORT}\n>long#LTR/Gypsy\n{LONG}\n>mid#LTR/Gypsy\n{MID}\n>other#LINE\n{OTHER}\n')
    # query target fident qstart qend qlen tstart tend tlen
    (tmp_path / 'hits.tsv').write_text(
        ''.join('\t'.join(map(str, (ID[h[0]], ID[h[1]]) + h[2:])) + '\n' for h in hits))


def run(run_script, tmp_path, *extra):
    run_script('mmseqs_to_clstr', '--fasta', tmp_path / 'lib.fa', '--hits', tmp_path / 'hits.tsv',
               '--out', tmp_path / 'dedup.fa', *extra)
    return headers(tmp_path / 'dedup.fa'), (tmp_path / 'dedup.fa.clstr').read_text()


def test_minus_strand_member_joins_with_its_coordinates(run_script, tmp_path):
    # mid on the minus strand of long, reported only as long -> mid
    setup(tmp_path, [('long', 'mid', 0.9, 51, 950, 1000, 900, 1, 900)])
    reps, clstr = run(run_script, tmp_path)
    assert reps == ['short#DNA', 'long#LTR/Gypsy', 'other#LINE']
    assert '900nt, >mid#LTR/Gypsy... at 1:900:51:950/-/90.00%' in clstr


def test_best_earlier_representative_wins(run_script, tmp_path):
    # short hits long (85%) and other (95%); other is a representative before it, so it wins
    setup(tmp_path, [('short', 'long', 0.85, 1, 500, 500, 1, 500, 1000),
                     ('short', 'other', 0.95, 1, 500, 500, 101, 600, 950)])
    reps, clstr = run(run_script, tmp_path)
    assert 'short#DNA' not in reps
    block = clstr.split('>Cluster')[3]
    assert '>other#LINE... *' in block and '>short#DNA... at 1:500:101:600/+/95.00%' in block


def test_member_of_a_member_stays_a_representative_if_unlinked(run_script, tmp_path):
    # mid joins long; short aligns only to mid (not a representative) so it stands alone
    setup(tmp_path, [('mid', 'long', 0.9, 1, 900, 900, 1, 900, 1000),
                     ('short', 'mid', 0.99, 1, 500, 500, 1, 500, 900)])
    reps, _ = run(run_script, tmp_path)
    assert reps == ['short#DNA', 'long#LTR/Gypsy', 'other#LINE']


def test_identity_and_coverage_thresholds(run_script, tmp_path):
    setup(tmp_path, [('mid', 'long', 0.79, 1, 900, 900, 1, 900, 1000),   # identity too low
                     ('short', 'long', 0.9, 1, 350, 500, 1, 350, 1000)])  # covers 70% of short
    reps, _ = run(run_script, tmp_path)
    assert len(reps) == 4
    reps, _ = run(run_script, tmp_path, '--identity', '0.75', '--coverage', '0.7')
    assert reps == ['long#LTR/Gypsy', 'other#LINE']


def test_a_hit_on_an_unknown_record_number_stops(tmp_path):
    setup(tmp_path, [])
    (tmp_path / 'hits.tsv').write_text('s1\ts9\t0.9\t1\t500\t500\t1\t500\t1000\n')
    r = subprocess.run([sys.executable, str(BIN / 'mmseqs_to_clstr.py'), '--fasta', tmp_path / 'lib.fa',
                        '--hits', tmp_path / 'hits.tsv', '--out', tmp_path / 'dedup.fa'],
                       capture_output=True, text=True)
    assert r.returncode != 0 and 'record number' in r.stderr


# --- the real MMseqs2 boundary: names MMseqs2 would rewrite on export -------------------

def mutate(seq, rate, rng):
    return ''.join(rng.choice('ACGT'.replace(b, '')) if rng.random() < rate else b for b in seq)


def revcomp(seq):
    return seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def rename_script():
    """The module's renaming command, cut out from between its markers."""
    block = re.search(r'# mmseqs-rename-begin\n(.*?)\n\s*# mmseqs-rename-end', MODULE.read_text(), re.S)
    return (block.group(1).replace('\\\\', '\\').replace('\\$', '$')
            .replace('${prefix}.combined_library.fa', 'lib.fa'))


@pytest.mark.skipif(shutil.which('mmseqs') is None or shutil.which('awk') is None,
                    reason='needs mmseqs and awk')
def test_real_mmseqs_search_keeps_names_it_would_rewrite(tmp_path):
    rng = random.Random(7)
    x = ''.join(rng.choice('ACGT') for _ in range(1000))
    y = ''.join(rng.choice('ACGT') for _ in range(1000))
    # MMseqs2 exports 'consensus_1' as '1', which is also a real name here, and parses 'lcl|'
    lib = {'consensus_1#DNA/hAT': x,
           '1#LINE/L1': y,
           'lcl|fam#DNA/hAT': revcomp(mutate(x[50:950], 0.05, rng)),   # minus strand of consensus_1
           'plain#LINE/L1': mutate(y[:900], 0.05, rng)}               # plus strand of 1
    (tmp_path / 'lib.fa').write_text(''.join(f'>{k}\n{v}\n' for k, v in lib.items()))
    sh = lambda cmd: subprocess.run(['bash', '-c', cmd], cwd=tmp_path, check=True, capture_output=True)
    sh(rename_script())
    sh('mmseqs createdb mm_input.fa mm_db --dbtype 2 -v 1')
    # the module's search settings
    sh('mmseqs search mm_db mm_db mm_res mm_tmp --search-type 3 --strand 2 -a --max-seqs 2000 -e 10 '
       '--min-seq-id 0.8 --threads 2 -v 1')
    sh('mmseqs convertalis mm_db mm_db mm_res mm_hits.tsv '
       '--format-output query,target,fident,qstart,qend,qlen,tstart,tend,tlen -v 1')
    subprocess.run([sys.executable, str(BIN / 'mmseqs_to_clstr.py'), '--fasta', 'lib.fa',
                    '--hits', 'mm_hits.tsv', '--out', 'dedup.fa'], cwd=tmp_path, check=True)
    assert headers(tmp_path / 'dedup.fa') == ['consensus_1#DNA/hAT', '1#LINE/L1']
    blocks = (tmp_path / 'dedup.fa.clstr').read_text().split('>Cluster')[1:]
    assert '>consensus_1#DNA/hAT... *' in blocks[0] and '>lcl|fam#DNA/hAT... at ' in blocks[0]
    assert '/-/' in blocks[0]
    assert '>1#LINE/L1... *' in blocks[1] and '>plain#LINE/L1... at ' in blocks[1] and '/+/' in blocks[1]
