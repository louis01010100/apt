from apt import utils
from pathlib import Path
import filecmp


def test_merge_static_column_file(tmp_path):
    fixture_dir = Path(
        __file__).parents[0] / 'fixture_merge_static_column_file'

    default_file = fixture_dir / 'default.tsv'
    mod_file = fixture_dir / 'mod.tsv'
    expected_file = fixture_dir / 'merged.tsv'
    expected_filtered_file = fixture_dir / 'merged_filtered.tsv'

    actual_file = tmp_path / 'merged.tsv'

    utils.merge_static_column_file(
        default_file,
        mod_file,
        actual_file,
    )

    assert filecmp.cmp(expected_file, actual_file, shallow=False)

    actual_filtered_file = tmp_path / 'merged_filtered.tsv'

    utils.merge_static_column_file(
        default_file,
        mod_file,
        actual_filtered_file,
        target_probesets={'AX-200', 'AX-300', 'AX-400'},
    )
    assert filecmp.cmp(
        expected_filtered_file,
        actual_filtered_file,
        shallow=False,
    )
