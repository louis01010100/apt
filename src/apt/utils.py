import gzip
import json
import logging
import pickle
import re
import shutil
import sys
from pathlib import Path

import humanfriendly
import numpy as np
import pandas as pd

from . import config

PROBESET_ID_PTN = re.compile(r'^((?:AX|AFFX-SP|AFFX-NP)-\d+).*$')


def find_apt_cmds(apt_cmds, bin_dirs):

    def find_cmd_path(cmd, bin_dir_list):
        for bin_dir in bin_dir_list:
            cmd_path = shutil.which(cmd, path=bin_dir)
            if cmd_path:
                return cmd_path

        return None

    apt_cmd_paths = {}
    for cmd in apt_cmds:
        cmd_path = find_cmd_path(cmd, bin_dirs)
        if not cmd_path:
            error(f'Error: cannot find needed APT program {cmd}, '
                  f'set the environment variable AXIOM_TOP to its location,'
                  f'or specify the APT bin folder using --apt-bin-dir,')
        apt_cmd_paths[cmd] = cmd_path
    return apt_cmd_paths


def euclidean_distance(v1, v2):
    return np.sqrt(np.dot(v1, v1) + np.dot(v2, v2))


def cosine(v1, v2):

    x = np.dot(v1, v2)
    y = np.sqrt(np.dot(v1, v1)) * np.sqrt(np.dot(v2, v2))

    return x / y


def init_logging(log_file: Path):

    h0 = logging.StreamHandler(sys.stderr)
    h0.setLevel(logging.INFO)
    h1 = logging.FileHandler(log_file, 'w')
    h1.setLevel(logging.DEBUG)

    logging.basicConfig(
        level=logging.INFO,
        handlers=[h0, h1],
        format='[%(levelname)s] %(asctime)s\t%(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )


def export_probesets(probesets, filepath):
    with filepath.open('wt') as fd:
        fd.write('probeset_id\n')
        for x in probesets:
            fd.write(f'{x}\n')


def random_sample(category, probesets, n, output_dir):
    random_ps_file = output_dir / f'{category}.{n}.ps'

    sample = probesets.sample(min(n, len(probesets)))
    sample.to_csv(
        random_ps_file,
        header=True,
        index=False,
    )

    return random_ps_file


def progress(msg):
    print(msg, file=sys.stderr)


def error(message):
    logging.error(message)
    sys.exit(message)


def save(obj, filepath):
    with gzip.open(filepath, 'wb') as fd:
        pickle.dump(obj, fd)


def load(filepath):
    with gzip.open(filepath, 'rb') as fd:
        return pickle.load(fd)


def get_plate_format(samples):

    x_list = []
    y_list = []

    for _, row in samples.iterrows():
        x = ord(row['well'][0]) - ord('A') + 1
        y = int(row['well'][1:])
        x_list.append(x)
        y_list.append(y)

    if max(x_list) > 8 or max(y_list) > 12:
        return {'height': 16, 'width': 24}
    else:
        return {'height': 8, 'width': 12}


def wellToXY(well, plate_format):
    if len(well) != 3:
        error(f"Error: illegal well specified {well}")
    x = int(well[1:3]) - 1
    y = (plate_format['height'] - 1) - (ord(well[0]) - ord('A'))
    return (x, y)


def get_point_size(n_samples):
    return f'ps {max(0.25, n_samples * -0.00025 + 1):.3f}'


def eval_sample_quality(
    samples,
    dqc_threshold,
    qccr_threshold,
    avg_qccr_threshold,
    min_samples_for_plate_qccr,
):

    qc_result_list = []
    qc_note_list = []
    idx_list = []

    for idx, row in samples.iterrows():
        idx_list.append(idx)

        if row['dqc'] < dqc_threshold:
            qc_result_list.append('FAIL')
            qc_note_list.append(f'DQC < {dqc_threshold}')
            continue

        if row['qccr'] < qccr_threshold:
            qc_result_list.append('FAIL')
            qc_note_list.append(f'QCCR < {qccr_threshold}%')
            continue

        if row['n_samples_qccr_ok'] < min_samples_for_plate_qccr:
            qc_result_list.append('PASS')
            qc_note_list.append('.')
            continue

        if row['avg_qccr'] < avg_qccr_threshold:
            qc_result_list.append('FAIL')
            qc_note_list.append(
                f'Average QCCR of passing samples < {avg_qccr_threshold}%')
            continue

        qc_result_list.append('PASS')
        qc_note_list.append('.')

    samples = samples.join(
        pd.DataFrame(
            {
                'qc_result': qc_result_list,
                'qc_note': qc_note_list,
            },
            index=idx_list,
        ))

    return samples


def call_rate_from_report_file(axiom_gt1_report_file):

    call_rates = pd.read_csv(
        axiom_gt1_report_file,
        comment='#',
        header=0,
        sep='\t',
    ).rename(columns={'cel_files': 'cel_name'})[['cel_name', 'call_rate']]

    return call_rates


def call_rate_from_calls_file(axiom_gt1_calls_file):
    calls = pd.read_csv(
        axiom_gt1_calls_file,
        comment='#',
        header=0,
        sep='\t',
        index_col=0,
    )

    samples = []
    call_rates = []

    for sample, callcodes in calls.items():
        n_total = len(callcodes)
        n_calls = len(callcodes[lambda x: x != -1])

        cr = n_calls / n_total * 100

        samples.append(sample)
        call_rates.append(cr)

    return pd.DataFrame({'cel_name': samples, 'call_rate': call_rates})


def find_files(root_dir_path, file_pattern, missing_ok=False):

    filepaths = [x.resolve() for x in root_dir_path.glob(file_pattern)]

    if len(filepaths) > 0:
        return filepaths

    if missing_ok:
        return []

    raise Exception(
        f'{root_dir_path} does not contain file with pattern {file_pattern}')


def find_file(root_dir_path, file_pattern, missing_ok=False):
    found_files = find_files(root_dir_path, file_pattern, missing_ok)

    if len(found_files) > 1:
        found_file = found_files[0]
        logging.warning(
            f'{root_dir_path} contains more than one file with pattern {file_pattern }. {found_file} is used'
        )
    elif len(found_files) == 1:
        found_file = found_files[0]
    else:
        if missing_ok:
            found_file = None
        else:
            raise Exception(
                f'{root_dir_path} does not contain file with pattern {file_pattern}'
            )

    return found_file


def export_cels(cels: pd.DataFrame, filepath):
    cels.sort_values(by=['cel_order'])[['cel_path']].to_csv(
        filepath,
        header=['cel_files'],
        index=False,
    )


def df2md(df, n=10, index=False, tablefmt='simple'):
    print(df.head(n).to_markdown(index=index, tablefmt=tablefmt))


def df2tsv(
    df,
    filepath,
    header=True,
    index=False,
    sep='\t',
    float_format=None,
    compress=False,
):

    if compress:
        filepath = filepath.with_suffix('.tsv.gz')
        with gzip.open(filepath, 'wt') as fh:
            df.to_csv(
                fh,
                header=header,
                index=index,
                sep=sep,
                float_format=float_format,
            )
    else:
        df.to_csv(
            filepath,
            header=header,
            index=index,
            sep=sep,
            float_format=float_format,
        )


def tsv2df(filepath):
    return pd.read_csv(
        filepath,
        comment='#',
        header=0,
        sep='\t',
        dtype='str',
    )


def get_array_name(lib_dir: Path):

    ax_package_file = find_file(lib_dir, '*.ax_package')

    with ax_package_file.open('r') as fh:
        ax_package = json.load(fh)

    return ax_package['array_name']


def log_start(
    analysis_type: str,
    version: str,
    workspace_dir: Path,
    lib_dir: Path = None,
    apt_bin_dir: Path = None,
):
    banner = f'axbp {analysis_type} v{version}'
    logging.info('#' * config.BANNER_WIDTH)
    logging.info(banner.center(config.BANNER_WIDTH))
    logging.info('#' * config.BANNER_WIDTH)
    if lib_dir:
        logging.info(f"lib_dir: {lib_dir}")
    if apt_bin_dir:
        logging.info(f"apt_bin_dir: {apt_bin_dir}")
    logging.info(f"workspace_dir: {workspace_dir}")


def log_stop(analysis_type, start_time, stop_time):
    logging.info('')
    banner = f'axbp {analysis_type} finished'
    logging.info('#' * config.BANNER_WIDTH)
    logging.info(banner.center(config.BANNER_WIDTH))
    logging.info('#' * config.BANNER_WIDTH)
    logging.info(f"start time: {start_time}")
    logging.info(f"stop time: {stop_time}")
    logging.info(
        f"duration: {humanfriendly.format_timespan(stop_time - start_time)}")


def cscores_iter(filepath):
    col2idx = {}
    with filepath.open('rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            samples = line.strip().split('\t')[1:]
            break

        for line in fh:
            line = line.strip()

            probeset_id, scores_str = line.split('\t', 1)

            yield {
                'probeset_id': probeset_id,
                'cscores': scores_str,
                'sample': samples,
            }


def posteriors_iter(filepath):
    col2idx = {}
    with filepath.open('rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            break

        current_probeset_id = None
        params = {}
        for line in fh:
            line = line.strip()

            psid_ploidy, params_str = line.split('\t', 1)
            items = psid_ploidy.split(':')
            probeset_id = items[0]
            if len(items) == 1:
                ploidy = '2'
            else:
                assert len(items) == 2
                ploidy = items[1]

            if not current_probeset_id:
                current_probeset_id = probeset_id
                params[ploidy] = params_str
            elif current_probeset_id == probeset_id:
                params[ploidy] = params_str
            else:
                assert current_probeset_id != probeset_id

                yield {
                    'probeset_id': current_probeset_id,
                    'params': params,
                }
                params = {}
                current_probeset_id = probeset_id
                params[ploidy] = params_str
        if params:
            yield {
                'probeset_id': current_probeset_id,
                'params': params,
            }


def summaries_iter(filepath):
    col2idx = {}
    with filepath.open('rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            samples = line.strip().split('\t')[1:]
            break

        current_probeset_id = None
        bag = {}
        for line in fh:
            line = line.strip()
            if line.startswith('AFFX-NP-'):
                continue

            psid_allele, signal_str = line.split('\t', 1)

            probeset_id = psid_allele[0:-2]
            allele = psid_allele[-1:]

            if not current_probeset_id:
                current_probeset_id = probeset_id
                bag[allele] = signal_str
            elif current_probeset_id == probeset_id:
                bag[allele] = signal_str
            else:
                assert current_probeset_id != probeset_id

                yield {
                    'probeset_id': current_probeset_id,
                    'signals': bag,
                    'samples': samples,
                }
                bag = {}
                current_probeset_id = probeset_id
                bag[allele] = signal_str
        if bag:
            yield {
                'probeset_id': current_probeset_id,
                'signals': bag,
                'samples': samples,
            }


def subset_file(input_file: Path, output_file: Path, probeset_ids: set):

    with input_file.open('rt') as ifh, output_file.open('wt') as ofh:
        for line in ifh:
            if line.startswith('#'):
                ofh.write(line)
                continue
            ofh.write(line)    # write header
            break

        for line in ifh:
            id_ = line.split('\t', 1)[0]

            match = PROBESET_ID_PTN.match(id_)

            if not match:
                raise Exception(line[0:78])

            probeset_id = match.group(1)
            if probeset_id in probeset_ids:
                ofh.write(line)


def _create_index(data_file):
    index = dict()
    with data_file.open('rt') as fh:
        header_found = False
        pos = None
        while True:
            pos = fh.tell()
            line = fh.readline()

            if not line:
                break

            if line.startswith('#'):
                continue

            if not header_found:
                header_found = True
                continue

            key = line.split('\t', 1)[0]
            index[key] = pos

            pos = fh.tell()

    return index


def merge_dynamic_column_file(
        default_file: Path,
        modified_file: Path,
        merged_file: Path,
        improved_probesets: set = set(),
        target_probesets: set = set(),
):
    default = pd.read_csv(
        default_file,
        comment='#',
        header=0,
        sep='\t',
        dtype='str',
    )

    modified = pd.read_csv(
        modified_file,
        comment='#',
        header=0,
        sep='\t',
        dtype='str',
    )

    default = default[lambda x: ~x['probeset_id'].isin(improved_probesets)]
    modified = modified[lambda x: x['probeset_id'].isin(improved_probesets)]

    data = pd.concat([default, modified])

    if target_probests:
        data = data[lambda x: x['probeset_id'].isin(target_probesets)]

    data.to_csv(merged_file, header=True, index=False, sep='\t')


def merge_static_column_file(
        default_file: Path,
        modified_file: Path,
        merged_file: Path,
        improved_probesets: set = set(),
        target_probesets: set = set(),
):

    index = _create_index(modified_file)

    with default_file.open('rt') as ifd, modified_file.open(
            'rt') as xfd, merged_file.open('wt') as ofd:

        for line in ifd:
            if line.startswith('#'):
                ofd.write(line)
                continue
            ofd.write(line)
            break

        for line in ifd:
            line = line.strip()
            key = line.split('\t', 1)[0]

            match = PROBESET_ID_PTN.match(key)
            assert match

            psid = match.group(1)

            pos = index.pop(key, None)

            if target_probesets and psid not in target_probesets:
                continue

            if pos:
                xfd.seek(pos)
                v = xfd.readline().strip()
            else:
                v = None

            if improved_probesets:
                if psid in improved_probesets:
                    assert v
                    ofd.write(v)
                else:
                    ofd.write(line)
            elif v:
                ofd.write(v)
            else:
                ofd.write(line)
            ofd.write('\n')

        for key in index:

            match = PROBESET_ID_PTN.match(key)
            assert match

            psid = match.group(1)

            if target_probesets and psid not in target_probesets:
                continue

            pos = index[key]
            xfd.seek(pos)
            v = xfd.readline()
            ofd.write(v)
