"""modules/local/edta/main.nf: the filter that removes the rice sequences EDTA substituted
for a module that found nothing. The test runs the module's own awk program, cut out from
between its rice-filter markers, so it checks the code the pipeline runs."""
import re
import shutil
import subprocess
from pathlib import Path

import pytest

MODULE = Path(__file__).resolve().parents[2] / 'modules' / 'local' / 'edta' / 'main.nf'
LIB = '>fam1#LTR/Gypsy\nACGTACGT\nACGT\n>rice1#SINE/tRNA\nGGGG\nGG\n>fam2#DNA/hAT\nTTTT\n'

pytestmark = pytest.mark.skipif(shutil.which('awk') is None, reason='needs awk')


def filter_script():
    """The shell lines between the markers, with Groovy's escapes undone."""
    block = re.search(r'# rice-filter-begin\n(.*?)\n\s*# rice-filter-end', MODULE.read_text(), re.S)
    return block.group(1).replace('\\$', '$').replace('\\\\', '\\')


def run(tmp_path, ids):
    (tmp_path / 'lib.fa').write_text(LIB)
    (tmp_path / 'rice_removed.ids').write_text(''.join(f'{i}\n' for i in ids))
    subprocess.run(['bash', '-c', 'LIB=lib.fa\n' + filter_script()], cwd=tmp_path, check=True)
    return (tmp_path / 'edta_norice.fa').read_text()


def test_empty_removal_list_keeps_the_whole_library(tmp_path):
    assert run(tmp_path, []) == LIB


def test_one_removal_keeps_every_other_record_and_its_sequence_lines(tmp_path):
    assert run(tmp_path, ['rice1#SINE/tRNA']) == '>fam1#LTR/Gypsy\nACGTACGT\nACGT\n>fam2#DNA/hAT\nTTTT\n'


def test_all_records_removed_gives_an_empty_library(tmp_path):
    assert run(tmp_path, ['fam1#LTR/Gypsy', 'rice1#SINE/tRNA', 'fam2#DNA/hAT']) == ''


def test_missing_removal_list_keeps_the_whole_library(tmp_path):
    (tmp_path / 'lib.fa').write_text(LIB)
    subprocess.run(['bash', '-c', 'LIB=lib.fa\n' + filter_script()], cwd=tmp_path, check=True)
    assert (tmp_path / 'edta_norice.fa').read_text() == LIB
