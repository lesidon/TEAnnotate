"""bin/rip_scan.py: parsing of RepeatMasker .align headers into TE copies."""
from conftest import load_script

rip = load_script('rip_scan')

# RepeatMasker 4.x: raw alignment ID, then the copy ID that joins a copy's fragments
PLUS = '1000 12.50 0.00 0.00 chr1 101 108 (892) famA#LTR/Gypsy 1 8 (92) c_b1s251i0 5'
MINUS = '900 10.00 0.00 0.00 chr1 201 208 (792) C famB#DNA/hAT (92) 8 1 m_b1s252i3 6'


def test_overlap_flag_does_not_change_copy_identity_or_family():
    assert rip.parse_header(PLUS + ' *') == rip.parse_header(PLUS) == \
        ('chr1', 101, 108, 'famA', 'LTR/Gypsy', '5', 12.5)


def test_minus_strand_header():
    assert rip.parse_header(MINUS) == ('chr1', 201, 208, 'famB', 'DNA/hAT', '6', 10.0)
    assert rip.parse_header(MINUS + ' *')[5] == '6'


def test_older_header_with_only_the_copy_id():
    assert rip.parse_header('1000 12.50 0.00 0.00 chr1 101 108 (892) famA#LTR/Gypsy 1 8 (92) 5')[5] == '5'


def test_unexpected_header_is_rejected_not_guessed():
    assert rip.parse_header('1000 12.50 0.00 0.00 chr1 101 108 (892) famA#LTR/Gypsy 1 8 (92)') is None
    assert rip.parse_header(PLUS + ' extra') is None
    assert rip.parse_header('1000 12.50 0.00 0.00 chr1 101 108 (892) famA#LTR/Gypsy 1 8 (92) c_b1s251i0') is None


def write_align(path, headers):
    aln = '\n  chr1          1 ACGTACGT 8\n\n  famA#LTR/Gypsy 1 ACGTACGT 8\n\n'
    path.write_text(''.join(h + '\n' + aln for h in headers))


def test_two_starred_ids_on_one_contig_stay_two_copies(run_script, tmp_path):
    other = '950 9.00 0.00 0.00 chr1 301 308 (692) famC#LINE/L1 1 8 (92) 7 *'
    write_align(tmp_path / 'g.align', [PLUS + ' *', other])
    (tmp_path / 'g.fa').write_text('>chr1\n' + 'ACGT' * 250 + '\n')
    run_script('rip_scan', '--genome', tmp_path / 'g.fa', '--align', tmp_path / 'g.align',
               '--prefix', tmp_path / 'g', '--min-copy-len', '1')
    rows = [l.split('\t') for l in (tmp_path / 'g.rip_copies.tsv').read_text().splitlines()[1:]]
    assert sorted((r[3], r[4]) for r in rows) == [('5', 'famA'), ('7', 'famC')]


def test_fragments_of_one_copy_are_joined(run_script, tmp_path):
    frag2 = '800 11.00 0.00 0.00 chr1 121 128 (872) famA#LTR/Gypsy 9 16 (84) 5'
    write_align(tmp_path / 'g.align', [PLUS, frag2])
    (tmp_path / 'g.fa').write_text('>chr1\n' + 'ACGT' * 250 + '\n')
    run_script('rip_scan', '--genome', tmp_path / 'g.fa', '--align', tmp_path / 'g.align',
               '--prefix', tmp_path / 'g', '--min-copy-len', '1')
    rows = [l.split('\t') for l in (tmp_path / 'g.rip_copies.tsv').read_text().splitlines()[1:]]
    assert len(rows) == 1 and rows[0][1:4] == ['101', '128', '5'] and rows[0][6] == '2'
