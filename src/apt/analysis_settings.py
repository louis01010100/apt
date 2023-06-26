import json
from enum import Enum
from pathlib import Path

Key = Enum('Key', [
    'SNP_PRIORS_FILENAME',
    'ARGS_FILENAME',
    'PROBESET_IDS_FILENAME',
    'RARE_HET_ADJUSTMENT',
    'PS2SNP_FILENAME',
    'GENOTYPE_FREQ_FILENAME',
    'PSCT_FILENAME',
    'ANNOTDB_FILENAME',
    'CN_ARGS_FILENAME',
    'CNREF_FILENAME',
    'CN_PRIORS_FILENAME',
])


class AnalysisSettings():

    def __init__(self, filepath: Path):
        with filepath.open('rt') as fh:
            settings = json.load(fh)

        self._values = settings

    @property
    def sample_qc(self):
        data = dict()
        for step in self._values['steps']:
            if step['name'] != 'sample_qc':
                continue

            analysis_files = step['analysis_files']

            analysis_files = step['analysis_files']

            data[Key.ARGS_FILENAME] = get_arg(
                analysis_files,
                'arg-file',
                'file_name',
            )

            data[Key.SNP_PRIORS_FILENAME] = get_arg(
                analysis_files,
                'snp-priors-input-file',
                'file_name',
            )
            data[Key.PROBESET_IDS_FILENAME] = get_arg(
                analysis_files,
                'probeset-ids',
                'file_name',
            )

        return data

    @property
    def genotype(self):
        data = dict()
        for step in self._values['steps']:
            if step['name'] != 'probeset_genotyping':
                continue

            analysis_files = step['analysis_files']

            data[Key.ARGS_FILENAME] = get_arg(
                analysis_files,
                'arg-file',
                'file_name',
            )

            data[Key.SNP_PRIORS_FILENAME] = get_arg(
                analysis_files,
                'snp-priors-input-file',
                'file_name',
            )
            data[Key.PROBESET_IDS_FILENAME] = get_arg(
                analysis_files,
                'probeset-ids',
                'file_name',
            )
            data[Key.PS2SNP_FILENAME] = get_arg(
                analysis_files,
                'ps2snp-file',
                'file_name',
            )
            data[Key.GENOTYPE_FREQ_FILENAME] = get_arg(
                analysis_files,
                'genotype-freq-file',
                'file_name',
            )

            data[Key.PSCT_FILENAME] = get_arg(
                analysis_files,
                'psct-file',
                'file_name',
            )
            data[Key.RARE_HET_ADJUSTMENT] = get_arg(
                analysis_files,
                'do-rare-het-adjustment',
                'enabled',
            )

            data[Key.CN_ARGS_FILENAME] = get_arg(
                analysis_files,
                'cn-arg-file',
                'file_name',
            )

            data[Key.CNREF_FILENAME] = get_arg(
                analysis_files,
                'reference-file',
                'file_name',
            )
            data[Key.CN_PRIORS_FILENAME] = get_arg(
                analysis_files,
                'cn-priors-file',
                'file_name',
            )

        return data

    @property
    def cnref(self):
        data = dict()
        for step in self._values['steps']:
            if step['name'] != 'copynumber_ref_file':
                continue

            analysis_files = step['analysis_files']

            data[Key.ARGS_FILENAME] = get_arg(
                analysis_files,
                'arg-file',
                'file_name',
            )

            data[Key.ANNOTDB_FILENAME] = get_arg(
                analysis_files,
                'apt-axiom-cnv-node:annotation-file',
                'file_name',
            )

        return data

    @property
    def cnvmix(self):
        data = dict()
        for step in self._values['steps']:
            if step['name'] == 'copy_number_fixed':
                analysis_files = step['analysis_files']

                data[Key.ARGS_FILENAME] = get_arg(
                    analysis_files,
                    'arg-file',
                    'file_name',
                )
                data[Key.ARGS_FILENAME] = get_arg(
                    analysis_files,
                    'cn_genotyping',
                    'file_name',
                )

        return data

    @property
    def array_type(self):
        return self._values['array_type']

    @property
    def is_default(self):
        return self._values['is_default']

    @property
    def version(self):
        return self._values['version']


def get_arg(analysis_files: dict, name: str, key: str):
    for analysis_file in analysis_files:
        if analysis_file['name'] == name:
            return analysis_file[key]

    return None
