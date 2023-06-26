import json
import logging
from pathlib import Path

from . import config, utils

PS_METRICS_THRESHOLDS = [
    'y-restrict', 'min-genotype-freq-samples', 'clustermin'
]


class Thresholds():

    def __init__(self, ax_thresholds_file):

        bag = dict()

        with ax_thresholds_file.open('rt') as fh:
            ax_thresholds_json = json.load(fh)
            for category in ax_thresholds_json['categories']:
                thresholds = category['thresholds']
                category_name = category['name']
                if category_name not in bag:
                    bag[category_name] = {}

                for threshold in thresholds:
                    threshold_name = threshold['name']
                    current_value = threshold['current_value']
                    bag[category_name][threshold_name] = current_value

        self._values = bag

    @property
    def dqc(self):
        return float(self._values['sample_qc']['axiom_dishqc_DQC'])

    @property
    def qccr(self):
        return float(self._values['sample_qc']['qc_call_rate'])

    @property
    def avg_qccr(self):
        return float(self._values['sample_qc']['plate_qc_averagecallrate'])

    @property
    def plate_pass_rate(self):
        return float(
            self._values['sample_qc']['plate_qc_percentsamplespassed'])

    @property
    def mapd(self):
        return float(self._values['cn_sample_qc']['mapd-max'])

    @property
    def waviness_sd(self):
        return float(self._values['cn_sample_qc']['waviness-sd-max'])

    # don't perform plate QCCR (step 6) if less than 20 samples
    @property
    def min_samples_for_plate_qccr(self):
        return 20

    @property
    def cnqcc(self):
        bag = []
        for k, v in self._values['cn_sample_qc'].items():
            bag.append(_arg2option(k, v) + '\n')
        return ''.join(bag)

    @property
    def cnqc(self):
        bag = []
        for k, v in self._values['cn_sample_qc'].items():
            if k in ['mapdc-max', 'waviness-sdc-max']:
                continue
            bag.append(_arg2option(k, v) + '\n')
        return ''.join(bag)

    @property
    def genotype_p_value(self):
        return self._values['snp_qc']['genotype-p-value-cutoff']

    @property
    def ps_classification(self):
        bag = []

        thresholds = dict()

        thresholds.update(self._values['snp_qc'])
        thresholds.update(self._values['MA_SNP_qc'])
        thresholds.update(self._values['ps_supplemental'])

        for k, v in thresholds.items():
            if k in PS_METRICS_THRESHOLDS:
                continue
            bag.append(_arg2option(k, v) + '\n')

        return ''.join(bag)

    @property
    def ps_metrics(self):
        bag = []

        thresholds = dict()

        thresholds.update(self._values['snp_qc'])
        thresholds.update(self._values['MA_SNP_qc'])
        thresholds.update(self._values['ps_supplemental'])

        for k, v in thresholds.items():
            if k in PS_METRICS_THRESHOLDS:
                bag.append(_arg2option(k, v) + '\n')

        return ''.join(bag)


def load_thresholds(
    libpath: Path,
    species: str,
    het_so_cutoff: float = None,
):

    filepaths = [x for x in libpath.glob('*.ax_thresholds')]

    if len(filepaths) == 0:
        if species == 'human':
            filepath = config.RESOURCES_DIR / 'Human.v5.ax_thresholds'
        elif species == 'diploid':
            filepath = config.RESOURCES_DIR / 'Diploid.v5.ax_thresholds'
        elif species == 'polyplid':
            filepath = config.RESOURCES_DIR / 'Polyploid.v5.ax_thresholds'
        else:
            utils.error(f'Unknown species: {species}')
    elif len(filepaths) == 1:
        filepath = filepaths[0]
    else:
        filepath = filepaths[0]
        logging.warning(
            f'{libpath} contains more than one file with pattern *.ax_thresholds; '
            f'{filepath} is used')
    thresholds = Thresholds(filepath)
    if het_so_cutoff:
        thresholds._values['snp_qc']['het-so-cutoff'] = str(het_so_cutoff)

    return thresholds


def _arg2option(k, v):

    if isinstance(v, str):
        v = v.replace(' ', '')
    return f'\t--{k}\t{v}'
