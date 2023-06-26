import gzip
import json
import logging
import os
import pickle
import re
import shutil
import sys
import xml.dom.minidom
from datetime import datetime
from enum import Enum
from importlib import resources
from pathlib import Path
from subprocess import PIPE, STDOUT, Popen

import pandas as pd

from . import config, utils


class Library():

    def __init__(
        self,
        libpath: Path,
    ):

        self.path = Path(libpath)
        ax_package_file = self.get_filepath('*.ax_package')

        with ax_package_file.open('r') as fh:
            ax_package = json.load(fh)

        self._ax_package = ax_package

        self._lib_set_version = self._ax_package['lib_set_version'].split(
            '.')[0]
        self._array_name = self._ax_package['array_name']

        self._step1_args = Step1Args(self.get_filepath('*Step1*.xml'))
        self._step2_args = Step2Args(self.get_filepath('*Step2*.xml'))

        if self.step2_args.process_multi_alleles:
            self.ps2snp_file = self.get_filepath(
                f'{self.name}.ps2multisnp_map.ps')
        else:
            self.ps2snp_file = self.get_filepath(f'{self.name}.ps2snp_map.ps')

    @property
    def name(self):
        return f'{self._array_name}.{self._lib_set_version}'

    @property
    def psct_file(self):
        return self._psct_file

    @property
    def cnref_file(self):
        return self.get_filepath('*.cn_models')

    @property
    def cn_ps1_args_file(self):
        return self.get_filepath('*apt-genotype-axiom.AxiomCN_PS1.apt2.xml')

    @property
    def cnvmix_args_file(self):
        return self.get_filepath(
            '*apt-copynumber-axiom-cnvmix.AxiomCNVmix.apt2.xml')

    @property
    def step1_args(self):
        return self._step1_args

    @property
    def step2_args(self):
        return self._step2_args

    @property
    def qc_args_file(self):
        return self.get_filepath(f'{self.name}.apt-geno-qc.AxiomQC1.xml')

    @property
    def array_type(self):
        return self._array_name

    @property
    def special_snps_file(self):
        return self.step2_args.special_snps_file

    @property
    def snp_params_file(self):
        return self.step2_args.snp_params_file

    @property
    def snp_priors_file(self):
        return self.step2_args.snp_priors_file

    @property
    def multi_priors_input_file(self):
        if self.step2_args.process_multi_alleles:
            return self.step2_args.multi_priors_input_file
        return None

    @property
    def process_multi_alleles(self):
        return self.step2_args.process_multi_alleles

    @property
    def igender_female_threshold(self):
        return self.step2_args.igender_female_threshold

    @property
    def igender_male_threshold(self):
        return self.step2_args.igender_male_threshold

    @property
    def annotdb_file(self):
        return self.get_filepath('*.annot.db', )

    @property
    def annotdc_file(self):
        return self.get_filepath('*.dc_annot.csv')

    @property
    def translation_file(self):
        return self.get_filepath('*.translation')

    @property
    def metabolizer_file(self):
        return self.get_filepath('*.metabolizer')

    @property
    def use_copynumber_call_codes(self):
        return self.step2_args.use_copynumber_call_codes

    def get_filepath(self, pattern, missing_ok=False):
        root_path = self.path
        filepaths = [x for x in root_path.glob(pattern)]

        if len(filepaths) == 0:
            if missing_ok:
                return None
            utils.error(
                f'{root_path} does not contain file with pattern {pattern}')

        filepath = filepaths[0]

        if len(filepaths) > 1:
            logging.warning(
                f'{root_path} contains more than one file with pattern {pattern}; '
                f'{filepath} is used')

        return filepath


class Step1Args():

    def __init__(self, args_file):
        self.args_file = args_file

        self.ps_file = args_file.parents[0] / self._arg(
            name='probeset-ids',
            analysis='library-file-node',
        )

        self.snp_priors_file = args_file.parents[0] / self._arg(
            name='snp-priors-input-file',
            analysis='genotyping-node',
        )

    def _arg(self, name, analysis):
        args_dom = xml.dom.minidom.parse(str(self.args_file))

        for element in args_dom.getElementsByTagName('Parameter'):
            current_name = element.getAttribute('name')
            current_analysis = element.getAttribute('analysis')
            value = element.getAttribute('currentValue')

            if current_name == name and current_analysis == analysis:
                return value
        utils.error(''
                    f'{self.args_file} does not contain '
                    f'parameter with name = {name} and analysis = {analysis}'
                    '')


class Step2Args():

    def __init__(self, args_file):

        self.args_file = args_file

        self.args_dom = xml.dom.minidom.parse(str(self.args_file))

        dir_path = args_file.parents[0]

        ps_file = self._arg(
            name='probeset-ids',
            analysis='library-file-node',
        )

        if ps_file:
            self.ps_file = dir_path / ps_file
        else:
            self.ps_file = None
        snp_params_file = self._arg(
            name='snp-specific-param-file',
            analysis='library-file-node',
        )

        if snp_params_file:
            snp_params_file = dir_path / snp_params_file

        self.snp_params_file = snp_params_file

        self.special_snps_file = dir_path / self._arg(
            name='special-snps',
            analysis='library-file-node',
        )

        self.snp_priors_file = dir_path / self._arg(
            name='snp-priors-input-file',
            analysis='genotyping-node',
        )

        self.process_multi_alleles = self._process_multi_alleles(
            self.args_dom)

        if self.process_multi_alleles:
            self.multi_priors_input_file = dir_path / self._arg(
                name='multi-priors-input-file',
                analysis='multi-genotyping-node',
            )

        self.igender_female_threshold = float(
            self._arg(
                'igender-female-threshold',
                'raw-gender-node',
            ))

        self.igender_male_threshold = float(
            self._arg(
                'igender-male-threshold',
                'raw-gender-node',
            ))

    @property
    def use_copynumber_call_codes(self):
        return self._arg(
            'use-copynumber-call-codes',
            'probeset-summarize-genotype-node',
        )

    @staticmethod
    def _process_multi_alleles(args_dom):
        for element in args_dom.getElementsByTagName('Parameter'):
            name = element.getAttribute('name')
            value = element.getAttribute('currentValue')

            if name != 'process-multi-alleles':
                continue

            if value == 'true':
                return True

        return False

    def _arg(self, name, analysis):
        args_dom = self.args_dom

        for element in args_dom.getElementsByTagName('Parameter'):
            current_name = element.getAttribute('name')
            current_analysis = element.getAttribute('analysis')
            value = element.getAttribute('currentValue')

            if current_name == name and current_analysis == analysis:
                return value

        return None
