from apt import utils
from pathlib import Path
import filecmp

#           defualt mod merged	merged_target	merged_target_improved
# AX-100    x           x       x               x
# AX-200    x       x   x       x               x
# AX-300    x       x   x       x               xx
# AX-400    x       x   x
# AX-500            x   x       x               x


def test_merge_static_column_file(tmp_path):
    fixture_dir = Path(
        __file__).parents[0] / 'fixture_merge_static_column_file'

    default_file = fixture_dir / 'default.tsv'
    mod_file = fixture_dir / 'mod.tsv'
    expected_file = fixture_dir / 'merged.tsv'
    expected_target_file = fixture_dir / 'merged_target.tsv'
    expected_target_improved_file = fixture_dir / 'merged_target_improved.tsv'

    actual_file = tmp_path / 'merged.tsv'

    utils.merge_static_column_file(
        default_file,
        mod_file,
        actual_file,
    )

    assert filecmp.cmp(expected_file, actual_file, shallow=False)

    actual_target_file = tmp_path / 'merged_target.tsv'

    utils.merge_static_column_file(
        default_file,
        mod_file,
        actual_target_file,
        target_probesets={'AX-100', 'AX-200', 'AX-300', 'AX-500'},
    )
    assert filecmp.cmp(
        expected_target_file,
        actual_target_file,
        shallow=False,
    )

    actual_target_improved_file = tmp_path / 'merged_target_improved.tsv'

    utils.merge_static_column_file(
        default_file,
        mod_file,
        actual_target_improved_file,
        improved_probesets={'AX-300'},
        target_probesets={'AX-100', 'AX-200', 'AX-300', 'AX-500'},
    )
    assert filecmp.cmp(
        expected_target_improved_file,
        actual_target_improved_file,
        shallow=False,
    )
