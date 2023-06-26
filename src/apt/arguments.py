import json
import xml
from enum import Enum
from pathlib import Path

from . import utils
from .thresholds import Thresholds

Key = Enum('Key', [
    'PROBESET_IDS_FILE',
    'SNP_PRIORS_FILE',
    'SNP_PARAMS_FILE',
    'GENDER_THRESHOLDS',
    'COPY_NUMBER_CALL_CODES',
    'PSCT_FILE',
    'CN_MODELS_FILE',
    'TRANSLATION_FILE',
    'METABOLIZER_FILE',
    'GENOTYPE_P_VALUE',
    'ALLELE_TRANSLATION',
    'CNAG',
    'RARE_HET_ADJUSTMENT',
    'SKIP_SNPOLISHER',
    'PNORM',
    'GENOTYPE_FREQ_FILE',
    'KNOWN_SUBOPTIMAL_PROBESETS_FILE',
    'RECALL_PROBESETS_FILE',
])


class SampleQcArguments():

    def __init__(
        self,
        lib_dir: Path,
        thresholds: Thresholds,
        geno_qc_args_file: Path = None,
        step1_args_file: Path = None,
    ):

        self._lib_dir = lib_dir
        self._thresholds = thresholds
        self._step1_args = dict()
        if not step1_args_file:
            step1_args_file = utils.find_file(lib_dir, '*Step1*.xml')
        if not geno_qc_args_file:
            geno_qc_args_file = utils.find_file(lib_dir,
                                                '*apt-geno-qc.AxiomQC1.xml')

        self._step1_args_file = step1_args_file
        self._geno_qc_args_file = geno_qc_args_file
        self._step1_args_dom = xml.dom.minidom.parse(
            str(self._step1_args_file))

    @property
    def lib_dir(self):
        return self._lib_dir

    @property
    def gender_thresholds(self):
        if Key.GENDER_THRESHOLDS in self._step1_args:
            return self._step1_args[Key.GENDER_THRESHOLDS]

        female_threshold = _get_arg(
            self._step1_args_dom,
            'igender-female-threshold',
            'raw-gender-node',
            converter=float,
        )
        male_threshold = _get_arg(
            self._step1_args_dom,
            'igender-male-threshold',
            'raw-gender-node',
            converter=float,
        )

        return {'male': male_threshold, 'female': female_threshold}

    @gender_thresholds.setter
    def gender_thresholds(self, male_threshold: float,
                          female_threshold: float):
        self._step1_args[Key.GENDER_THRESHOLD] = {
            'male': male_threshold,
            'female': female_threshold
        }

    @property
    def dqc_threshold(self):
        return self._thresholds.dqc

    @property
    def qccr_threshold(self):
        return self._thresholds.qccr

    @property
    def avg_qccr_threshold(self):
        return self._thresholds.avg_qccr

    @property
    def plate_pass_rate_threshold(self):
        return self._thresholds.plate_pass_rate

    @property
    def min_samples_for_plate_qccr_threshold(self):
        return self._thresholds.min_samples_for_plate_qccr

    @property
    def geno_qc_args_file(self):
        return self._geno_qc_args_file

    @property
    def step1_args_file(self):
        return self._step1_args_file

    @property
    def probeset_ids_file(self):
        if Key.PROBESET_IDS_FILE in self._step1_args:
            return self._step1_args[Key.PROBESET_IDS_FILE]

        return self._lib_dir / _get_arg(
            self._step1_args_dom,
            'probeset-ids',
            'libary-file-node',
        )

    @probeset_ids_file.setter
    def probeset_ids_file(self, filepath):
        self._step1_args[Key.PROBESET_IDS_FILE] = filepath

    @property
    def geno_qc_args_file(self):
        return self._geno_qc_args_file

    @property
    def step1_args_file(self):
        return self._step1_args_file

    @property
    def snp_priors_file(self):
        if Key.SNP_PRIORS_FILE in self._step1_args:
            return self._step1_args[Key.SNP_PRIORS_FILE]

        return self._lib_dir / _get_arg(
            self._step1_args_dom,
            'snp-priors-input-file',
            'genotyping-node',
        )

    @snp_priors_file.setter
    def snp_priors_file(self, filepath):
        self._step1_args[Key.SNP_PRIORS_FILE] = filepath


class SnvArguments():

    def __init__(
        self,
        lib_dir: Path,
        thresholds: Thresholds,
        args_file: Path = None,
    ):

        if not args_file:
            args_file = utils.find_file(lib_dir, '*Step2*.xml')

        self._args_file = args_file
        self._args_dom = xml.dom.minidom.parse(str(self._args_file))
        self._args = dict()
        self._lib_dir = lib_dir
        self._thresholds = thresholds

    @property
    def lib_dir(self):
        return self._lib_dir

    @property
    def psct_file(self):
        if Key.PSCT_FILE in self._args:
            return self._args[Key.PSCT_FILE]

        return utils.find_file(self._lib_dir, '*psct', missing_ok=True)

    @psct_file.setter
    def psct_file(self, psct_file: Path):
        self._args[Key.PSCT_FILE] = psct_file

    @property
    def ps2snp_file(self):
        if self.multi_alleles:
            return utils.find_file(
                self._lib_dir,
                '*.ps2multisnp_map.ps',
            )
        else:
            return utils.find_file(self._lib_dir, '*.ps2snp_map.ps')

    @property
    def genotype_p_value(self):
        if Key.GENOTYPE_P_VALUE in self._args:
            return self._args[Key.GENOTYPE_P_VALUE]
        return self._thresholds.genotype_p_value

    @genotype_p_value.setter
    def genotype_p_value(self, p_value):
        self._args[Key.GENOTYPE_P_VALUE] = p_value

    @property
    def ps_metrics_thresholds(self):
        return self._thresholds.ps_metrics

    @property
    def ps_classification_thresholds(self):
        return self._thresholds.ps_classification

    @property
    def probeset_ids_file(self):
        if Key.PROBESET_IDS_FILE in self._args:
            return self._args[Key.PROBESET_IDS_FILE]
        return self._lib_dir / _get_arg(self._args_dom, 'probeset-ids',
                                        'library-file-node')

    @probeset_ids_file.setter
    def probeset_ids_file(self, filepath: Path):
        self._args[Key.PROBESET_IDS_FILE] = filepath

    @property
    def snp_priors_file(self):
        if Key.SNP_PRIORS_FILE in self._args:
            return self._args[Key.SNP_PRIORS_FILE]
        return self._lib_dir / _get_arg(
            self._args_dom, 'snp-priors-input-file', 'genotyping-node')

    @snp_priors_file.setter
    def snp_priors_file(self, filepath):
        self._args[Key.SNP_PRIORS_FILE] = filepath

    @property
    def snp_params_file(self):
        if Key.SNP_PARAMS_FILE in self._args:
            return self._args[Key.SNP_PARAMS_FILE]
        return self._lib_dir / _get_arg(
            self._args_dom, 'snp-specific-param-file', 'library-file-node')

    @snp_params_file.setter
    def snp_params_file(self, filepath):
        self._args[Key.SNP_PARAMS_FILE] = filepath

    @property
    def special_snps_file(self):
        return self._lib_dir / _get_arg(self._args_dom, 'special-snps',
                                        'library-file-node')

    @property
    def copynumber_call_codes(self):
        if Key.COPY_NUMBER_CALL_CODES in self._args:
            return self._args[Key.COPY_NUMBER_CALL_CODES]

        return _get_arg(
            self._args_dom,
            'use-copynumber-call-codes',
            'probeset-summarize-genotype-node',
            converter=lambda x: x.lower() == 'true',
            missing_ok=True,
        )

    @copynumber_call_codes.setter
    def copynumber_call_codes(self, value: bool):
        self._args[Key.COPY_NUMBER_CALL_CODES] = value

    @property
    def cnag(self):
        if Key.CNAG in self._args:
            return self._args[Key.CNAG]

    @cnag.setter
    def cnag(self, value: bool):
        self._args[Key.CNAG] = value
        if value:
            self._args[Key.COPY_NUMBER_CALL_CODES] = value

    @property
    def pnorm(self):
        if Key.PNORM in self._args:
            return self._args[Key.PNORM]

    @pnorm.setter
    def pnorm(self, value: bool):
        self._args[Key.PNORM] = value

    @property
    def multi_priors_input_file(self):
        return self.lib_dir / _get_arg(
            self._args_dom,
            'multi-priors-input-file',
            'multi-genotyping-node',
        )

    @property
    def multi_alleles(self):
        for element in self._args_dom.getElementsByTagName('Parameter'):
            name = element.getAttribute('name')
            value = element.getAttribute('currentValue')

            if name != 'process-multi-alleles':
                continue

            if value == 'true':
                return True

        return False

    @property
    def rare_het_adjustment(self):
        if Key.RARE_HET_ADJUSTMENT in self._args:
            return self._args[Key.RARE_HET_ADJUSTMENT]
        else:
            return False

    @rare_het_adjustment.setter
    def rare_het_adjustment(self, value):
        self._args[Key.RARE_HET_ADJUSTMENT] = value

    @property
    def snpolisher(self):
        if Key.SKIP_SNPOLISHER in self._args:
            return self._args[Key.SKIP_SNPOLISHER]
        else:
            return False

    @snpolisher.setter
    def snpolisher(self, value):
        self._args[Key.SKIP_SNPOLISHER] = value

    @property
    def step2_args_file(self):
        return self._args_file

    @property
    def genotype_freq_file(self):
        if Key.GENOTYPE_FREQ_FILE in self._args:
            return self._args[Key.GENOTYPE_FREQ_FILE]
        return utils.find_file(self._lib_dir,
                               '.genotype_frequency.txt',
                               missing_ok=True)

    @genotype_freq_file.setter
    def genotype_freq_file(self, filepath: Path):
        self._args[Key.GENOTYPE_FREQ_FILE] = filepath


class CnvArguments():

    def __init__(
        self,
        lib_dir,
        thresholds,
        cn_controls_file=None,
    ):

        self._lib_dir = lib_dir
        self._cn_controls_file = cn_controls_file
        self._thresholds = thresholds

    @property
    def cn_controls_file(self):
        return self._cn_controls_file

    @property
    def lib_dir(self):
        return self._lib_dir

    @property
    def cnqc_thresholds(self):
        return self._thresholds.cnqc

    @property
    def cnqcc_thresholds(self):
        return self._thresholds.cnqcc

    @property
    def cn_ps1_args_file(self):
        filepath = utils.find_file(
            self._lib_dir,
            '*apt-genotype-axiom.AxiomCN_PS1.apt2.xml',
            missing_ok=True,
        )

        if not filepath:
            filepath = self.cn_gt1_args_file

        return filepath

    @property
    def cn_gt1_args_file(self):
        return utils.find_file(self._lib_dir,
                               '*apt-genotype-axiom.AxiomCN_GT1.apt2.xml')

    @property
    def cnvmix_args_file(self):
        return utils.find_file(
            self._lib_dir,
            '*apt-copynumber-axiom-cnvmix.AxiomCNVmix.apt2.xml')

    @property
    def cnvhmm_args_file(self):
        return utils.find_file(self._lib_dir,
                               '*apt-copynumber-axiom-hmm.AxiomHMM.apt2.xml')

    @property
    def cnref_args_file(self):
        return utils.find_file(
            self._lib_dir, r'*.apt-copynumber-axiom-ref.AxiomCNref.apt2.xml')

    @property
    def cnref_file(self):
        return utils.find_file(self._lib_dir, r'*.cn_models')

    @property
    def annotdb_file(self):
        return utils.find_file(self._lib_dir, r'*.annot.db')

    @property
    def special_snps_file(self):
        return utils.find_file(self._lib_dir, r'*.specialSNPs')

    @property
    def cn_models_template_file(self):
        return utils.find_file(self._lib_dir, r'*.cn_models_template')

    @property
    def mapd_threshold(self):
        return self._thresholds.mapd

    @property
    def waviness_sd_threshold(self):
        return self._thresholds.waviness_sd


class TranslationArguments():

    def __init__(self, lib_dir):
        self._lib_dir = lib_dir
        self._args = dict()

    @property
    def lib_dir(self):
        return self._lib_dir

    @property
    def translation_file(self):
        if Key.TRANSLATION_FILE in self._args:
            return self._args[Key.TRANSLATION_FILE]

        return utils.find_file(self._lib_dir, '*.translation')

    @translation_file.setter
    def translation_file(self, filepath: Path):
        self._args[Key.TRANSLATION_FILE] = filepath

    @property
    def metabolizer_file(self):
        if Key.METABOLIZER_FILE in self._args:
            return self._args[Key.METABOLIZER_FILE]

        return utils.find_file(self._lib_dir, '*.metabolizer')

    @metabolizer_file.setter
    def metabolizer_file(self, filepath: Path):
        self._args[Key.METABOLIZER_FILE] = filepath

    @property
    def dc_annot_file(self):
        return utils.find_file(self._lib_dir, '*.dc_annot.csv')

    @property
    def allele_translation(self):
        if Key.ALLELE_TRANSLATION in self._args:
            return self._args[Key.ALLELE_TRANSLATION]
        return False

    @allele_translation.setter
    def allele_translation(self, value):
        self._args[Key.ALLELE_TRANSLATION] = value


class FormatArguments():

    def __init__(
        self,
        lib_dir: Path,
        export_vcf: bool = False,
        export_plink: bool = False,
        export_igv: bool = False,
        export_axas: bool = False,
    ):
        self._lib_dir = lib_dir
        self._export_vcf = export_vcf
        self._export_plink = export_plink
        self._export_igv = export_igv
        self._export_axas = export_axas

    @property
    def annotdb_file(self):
        return utils.find_file(self._lib_dir, '*annot.db')

    @property
    def ax_thresholds_dir(self):
        return utils.find_file(self._lib_dir, '*annot.db')

    @property
    def export_igv(self):
        return self._export_igv

    @property
    def export_axas(self):
        return self._export_igv

    @property
    def export_vcf(self):
        return self._export_vcf

    @property
    def export_plink(self):
        return self._export_plink


def _get_arg(arg_dom, name, analysis, converter=str, missing_ok=False):
    for element in arg_dom.getElementsByTagName('Parameter'):
        current_name = element.getAttribute('name')
        current_analysis = element.getAttribute('analysis')
        value = element.getAttribute('currentValue')

        if current_name == name and current_analysis == analysis:
            return converter(value)

    if missing_ok:
        return None

    raise Exception(
        f'parameter with name = {name} and analysis = {analysis} not found')
