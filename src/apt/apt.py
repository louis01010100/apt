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
from pathlib import Path
from subprocess import PIPE, STDOUT, Popen

import pandas as pd

from . import config, utils


class Apt():

    def __init__(self, apt_bin_dir=None):

        bin_dirs = []

        if apt_bin_dir:
            bin_dirs.append(apt_bin_dir)

        bin_dirs.extend(config.BIN_DIRS)
        bin_dirs.append(os.environ['PATH'])

        apt_cmds = utils.find_apt_cmds(config.APT_PROGRAMS, bin_dirs)

        self.cmd_paths = apt_cmds

        # test to see if enscript works, if not we will use paps.  this primarily an issue
        # under Docker, since enscript refuses to run if the current user is not in /etc/passwd
        if shutil.which('enscript'):
            self._enscript_ok = True
        else:
            self._enscript_ok = False

    @property
    def apt_bin_dir(self):
        return Path(self.cmd_paths['apt-genotype-axiom']).parents[0]

    def export_cnv_axas(
        self,
        cnv_dir: Path,
        output_dir: Path,
        ax_thresholds_dir: Path,
    ):
        executable = self.cmd_paths['apt-package-util']

        log_file = output_dir / 'apt-package-util.log'

        cmd = (''
               f'{executable}\n'
               f'   --copynumber-data-dir       {cnv_dir}\n'
               f'   --batch-folder              {output_dir}\n'
               f'   --ax-thresholds-dir         {ax_thresholds_dir}\n'
               f'   --log-file                  {log_file}\n'
               '')

        return Command(cmd, check_log=True)

    def export_cnv_igv(
        self,
        cndata_file: Path,
        annotdb_file: Path,
        output_dir: Path,
    ):
        executable = self.cmd_paths['apt-format-result']

        log_file = output_dir / 'apt-format-result.log'

        cmd = (''
               f'{executable}\n'
               f'   --run-format-igv             true\n'
               f'   --annotation-file            {annotdb_file}\n'
               f'   --igv-cndata-file            {cndata_file}\n'
               f'   --igv-export-cnv-copynumber  true \n'
               f'   --igv-export-cnv-log2ratio   true \n'
               f'   --igv-export-cnv-baf         true \n'
               f'   --igv-export-cnv-smooth-signal true \n'
               f'   --igv-export-seg-cn          true \n'
               f'   --igv-export-seg-loh         true \n'
               f'   --export-single-samples      true \n'
               f'   --igv-out-dir                {output_dir}\n'
               f'   --log-file                   {log_file}\n'
               '')

        return Command(cmd, check_log=True)

    def unpack_a5_file(self, filepath, output_dir):
        executable = self.cmd_paths['apt2-dset-util']

        log_file = output_dir / 'apt2-dset-util.log'

        cmd = (''
               f"{executable}\n"
               f"    --input-file  {filepath}\n"
               f"    --output-dir  {output_dir}\n"
               f"    --output-type txt\n"
               f"    --log-file    {log_file} \n"
               '')

        return Command(cmd, check_log=True)

    def apt2_dset_util(self, cel_path):
        executable = self.cmd_paths['apt2-dset-util']
        cmd = (''
               f"{executable}\n"
               f"    --dump-gdh    \n"
               f"    --log-file    /dev/null \n"
               f"    --input-file  {cel_path}\n"
               '')

        return Command(cmd, check_log=True)

    def apt_geno_qc(self, lib_dir, args_file, cel_files, output_dir, force):
        executable = self.cmd_paths['apt-geno-qc']
        cmd = (
            ''
            f"{executable}\n"
            f"    --cel-files            {cel_files}\n"
            f"    --analysis-files-path  {lib_dir}\n"
            f"    --xml-file             {args_file}\n"
            f"    --out-file             {output_dir / 'apt-geno-qc.txt'}\n"
            f"    --log-file             {output_dir / 'apt-geno-qc.log'}\n"
            '')
        if force:
            cmd += "    --force\n"

        return Command(cmd)

    def apt_copynumber_axiom_ref(
        self,
        cels_file: Path,
        lib_dir: Path,
        args_file: Path,
        output_dir: Path,
        cnqc_thresholds: str = None,
        force: bool = False,
        wave_correction: bool = None,
    ):

        log_file = output_dir / 'apt-copynumber-axiom-ref.log'
        ref_file = output_dir / config.CNREF_FILENAME

        executable = self.cmd_paths['apt-copynumber-axiom-ref']

        cmd = (''
               f"{executable}\n"
               f"    --cel-files                         {cels_file}\n"
               f"    --analysis-files-path               {lib_dir}\n"
               f"    --arg-file                          {args_file}\n"
               f"    --out-dir                           {output_dir}\n"
               f"    --new-reference-file                {ref_file}\n"
               f"    --log-file                          {log_file}\n"
               '')
        if wave_correction is not None:
            cmd += f"    --use-wave-correction           {wave_correction}\n"
        if force:
            cmd += "    --force\n"
        if cnqc_thresholds:
            cnqc_report_file = output_dir / config.CNQC_REPORT_FILENAME
            cmd += f"    --qc-results-file          {cnqc_report_file}\n"
            cmd += cnqc_thresholds

        return Command(cmd, check_log=True)

    def create_cnref_from_summary(
        self,
        lib_dir: Path,
        args_file: Path,
        summary_file: Path,
        report_file: Path,
        calls_file: Path,
        output_dir: Path,
        annotdb_file: Path,
        special_snps_file: Path,
        cn_models_template_file: Path,
    ):
        log_file = output_dir / 'apt-copynumber-axiom-ref.log'
        ref_file = output_dir / config.CNREF_FILENAME

        executable = self.cmd_paths['apt-copynumber-axiom-ref']

        cmd = (
            ''
            f"{executable}\n"
            f"    --summary-file                         {summary_file}\n"
            f"    --cn-library-file-node:report-file     {report_file}\n"
            f"    --adjusted-intensities-node:calls-file {calls_file}\n"
            f"    --reference-file-template              {cn_models_template_file}\n"
            f"    --skip-genotyping                      True\n"
            f"    --out-dir                              {output_dir}\n"
            f"    --apt-axiom-cnv-node:annotation-file   {annotdb_file}\n"
            f"    --special-snps                         {special_snps_file}\n"
            f"    --new-reference-file                   {ref_file}\n"
            f"    --log-file                             {log_file}\n"
            '')

        return Command(cmd, check_log=True)

        # summary_file = gt_dir / 'AxiomGT1.summary.a5'
        # report_file = gt_dir / 'AxiomGT1.report.txt'
        # calls_file = gt_dir / 'AxiomGT1.calls.txt'
        # log_file = output_dir / 'apt-copynumber-axiom-ref.log'
        #
        # cmd = APT.apt_copynumber_axiom_ref + '\n'
        # cmd += ' --summary-file                         {}\n'.format(
        #     summary_file)
        # cmd += ' --cn-library-file-node:report-file     {}\n'.format(
        #     report_file)
        # cmd += ' --adjusted-intensities-node:calls-file {}\n'.format(
        #     calls_file)
        # cmd += ' --reference-file-template              {}\n'.format(
        #     LIB.cn_model_template)
        # cmd += ' --out-dir                              {}\n'.format(
        #     output_dir)
        # cmd += ' --new-reference-file                   {}\n'.format(
        #     new_ref_file)
        # cmd += ' --log-file                             {}\n'.format(log_file)
        # cmd += ' --skip-genotyping                    true\n'
        # cmd += ' --apt-axiom-cnv-node:annotation-file   {}\n'.format(
        #     LIB.annot_db)
        # cmd += ' --special-snps                         {}\n'.format(
        #     LIB.special_snps)
        #
        # execute(cmd)

    def apt_copynumber_axiom_hmm(
        self,
        summary_file,
        report_file,
        calls_file,
        confidences_file,
        output_dir,
        lib_dir,
        args_file,
        cnref_file,
    ):
        executable = self.cmd_paths['apt-copynumber-axiom-hmm']
        log_file = output_dir / 'apt-copynumber-axiom-hmm.log'

        cmd = (
            ''
            f"{executable}\n"
            f"    --summary-file                  {summary_file}\n"
            f"    --report-file                   {report_file}\n"
            f"    --out-dir                       {output_dir}\n"
            f"    --analysis-files-path           {lib_dir}\n"
            f"    --arg-file                      {args_file}\n"
        # f"    --use-text-format-summary-file  true\n"
            f"    --reference-file                {cnref_file}\n"
            f"    --log-file                      {log_file}\n"
            f"    --loh-calls-file                {calls_file}\n"
            f"    --loh-confidences-file          {confidences_file}\n"
            '')

        return Command(cmd, check_log=True)

    def apt_copynumber_axiom_cnvmix(
        self,
        summary_file,
        report_file,
        output_dir,
        lib_dir,
        args_file,
        cnref_file,
        cn_controls_file,
    ):
        executable = self.cmd_paths['apt-copynumber-axiom-cnvmix']

        cmd = (
            ''
            f"{executable}\n"
            f"    --summary-file           	{summary_file}\n"
            f"    --report-file            	{report_file}\n"
            f"    --out-dir                	{output_dir}\n"
            f"    --analysis-files-path    	{lib_dir}\n"
            f"    --arg-file               	{args_file}\n"
            f"    --reference-file         	{cnref_file}\n"
            f"    --use-text-format-summary-file true\n"
            f"    --log-file               	{output_dir / 'apt-copynumber-axiom-cnvmix.log'}\n"
            '')

        if cn_controls_file:
            cmd += f"    --controls-file        {cn_controls_file}\n"

        return Command(cmd, check_log=True)

    def export_signals(
        self,
        cel_files,
        output_dir,
        lib_dir,
        args_file,
        a5_format=False,
        force=False,
    ):

        executable = self.cmd_paths['apt-genotype-axiom']
        log_file = output_dir / 'apt-genotype-axiom.log'

        cmd = (
            ''
            f"{executable}\n"
            f"    --cel-files                                       {cel_files}\n"
            f"    --out-dir                                         {output_dir}\n"
            f"    --analysis-files-path                             {lib_dir}\n"
            f"    --arg-file                                        {args_file}\n"
            f"    --summaries-only                                  True\n"
            f"    --summary-a5-output                               {a5_format}\n"
            f"    --log-file                                        {log_file}\n"
            f"    --artifact-reduction-trustcheck                   True\n"
            f"    --artifact-reduction-output-trustcheck            True\n"
            f"    --dump-blemishes                                  True\n"
            '')

        if force:
            cmd += "    --force\n"

        return Command(cmd, check_log=True)

    def apt_package_util(
        self,
        genotype_data_dir: Path,
        performance_file: Path,
        batch_folder: Path,
        species: str,
        geno_qc_file: Path = None,
    ):

        executable = self.cmd_paths['apt-package-util']
        log_file = batch_folder / 'apt-package-util.log'
        workflow_type = 'probeset_genotyping'

        cmd = (''
               f"{executable}\n"
               f'    --log-file            {log_file}         \n'
               f'    --workflow-type       {workflow_type}    \n'
               f'    --species             {species}          \n'
               f'    --batch-folder        {batch_folder}     \n'
               f'    --genotype-data-dir   {genotype_data_dir}\n'
               f'    --performance-file    {performance_file} \n'
               f'    --geno-qc-res-file    {geno_qc_file}     \n'
               '')
        if geno_qc_file and geno_qc_file.exists():
            cmd += f'    --geno-qc-res-file    {geno_qc_file}     \n'

        return Command(cmd, check_log=True)

    def apt_dmet_translation(
        self,
        axas_dir: Path,
        cn_region_calls_file: Path,
        marker_list_file: Path,
        lib_dir: Path,
        translation_file: Path,
        metabolizer_file: Path,
        annotation_file: Path,
        output_dir: Path,
    ):

        executable = self.cmd_paths['apt-dmet-translation']

        log_file = output_dir / 'apt-dmet-translation.log'

        cmd = (''
               f"{executable}\n"
               f'    --batch-folder {axas_dir} \n'
               f'    --cn-region-calls-file {cn_region_calls_file} \n'
               f'    --marker-list-file {marker_list_file} \n'
               f'    --analysis-files-path {lib_dir} \n'
               f'    --translate-file {translation_file} \n'
               f'    --metabolizer-file {metabolizer_file} \n'
               f'    --annotation-file {annotation_file}\n'
               f'    --use-first-dup-allele-def true \n'
               f'    --out-dir {output_dir} \n'
               f'    --log-file {log_file}\n'
               '')

        return Command(cmd, check_log=False)

    def apt_genotype_axiom(
        self,
        cel_files,
        output_dir,
        lib_dir,
        args_file,
        snp_priors_file=None,
        snp_params_file=None,
        probeset_ids_file=None,
        rare_het_adjustment=False,
        force=False,
        process_multi_alleles=None,
        use_copynumber_call_codes=None,
        cnpscalls_file=None,
        export_allele_summaries=None,
        dump_blemishes=False,
        export_trustcheck=False,
    ):

        executable = self.cmd_paths['apt-genotype-axiom']

        log_file = output_dir / 'apt-genotype-axiom.log'

        cmd = (
            ''
            f"{executable}\n"
            f"    --cel-files                                       {cel_files}\n"
            f"    --out-dir                                         {output_dir}\n"
            f"    --analysis-files-path                             {lib_dir}\n"
            f"    --arg-file                                        {args_file}\n"
            f"    --log-file                                        {log_file}\n"
            f"    --dual-channel-normalization                      True\n"
            f"    --genotyping-node:snp-posteriors-output           True\n"
            f"    --dump-blemishes                                  {dump_blemishes}\n"
            f"    --artifact-reduction-trustcheck                   {export_trustcheck}\n"
            f"    --artifact-reduction-output-trustcheck            {export_trustcheck}\n"
            '')

        if export_allele_summaries is not None:
            cmd += f"    --allele-summaries                         {export_allele_summaries}\n"

        if process_multi_alleles is not None:
            cmd += f"    --process-multi-alleles                               {process_multi_alleles}\n"
            cmd += f"    --multi-genotyping-node:multi-posteriors-output       {process_multi_alleles}\n"

        if force:
            cmd += "    --force\n"
        if snp_priors_file:
            cmd += f"    --genotyping-node:snp-priors-input-file    {snp_priors_file}\n"
        if snp_params_file:
            cmd += f"    --snp-specific-param-file                  {snp_params_file}\n"
        if probeset_ids_file:
            cmd += f"    --probeset-ids                             {probeset_ids_file}\n"
        if use_copynumber_call_codes is not None:
            cmd += f"    --use-copynumber-call-codes  {use_copynumber_call_codes} \n"
        if rare_het_adjustment:
            cmd += "    --do-rare-het-adjustment                       True\n"
            cmd += "    --genotyping-node:@low-sdDist-minor-cutoff-4   -10\n"
            cmd += "    --genotyping-node:@high-sdDist-minor-cutoff-4  5000\n"
            cmd += "    --genotyping-node:@low-sdDist-major-cutoff-4   -100\n"
            cmd += "    --genotyping-node:@high-sdDist-major-cutoff-4  5000\n"

        if cnpscalls_file and cnpscalls_file.exists():
            cmd += f"    --copynumber-probeset-calls-file {cnpscalls_file} \n"

        return Command(cmd, check_log=True)

    def apt_summary_genotype_axiom(
        self,
        args_file,
        summary_file,
        trustcheck_file,
        priors_file,
        params_file,
        special_snps_file,
        gender_file,
        lib_dir,
        output_dir,
        use_copynumber_call_codes,
        probeset_ids_file=None,
        rare_het_adjustment=False,
    ):
        log_file = output_dir / 'apt-summary-genotype-axiom.log'
        posteriors_file = output_dir / 'AxiomGT1.snp-posteriors.txt'

        executable = self.cmd_paths['apt-summary-genotype-axiom']

        # cmd = (
        #     ''
        #     f'{executable}                                                           \n'
        #     f'    --analysis-files-path                         {lib_dir}            \n'
        #     f"    --arg-file                                    {args_file}          \n"
        #     f'    --out-dir                                     {output_dir}         \n'
        #     f'    --log-file                                    {log_file}           \n'
        #     f'    --snp-posteriors-output                       true                 \n'
        #     f'    --genotyping-node:snp-posteriors-output-file  {posteriors_file}    \n'
        #     f'    --summary-input-file                          {summary_file}       \n'
        #     f'    --artifact-reduction-trustcheck-file          {trustcheck_file}    \n'
        #     f'    --read-genders                                {gender_file}        \n'
        #     f'    --genotyping-node:snp-priors-input-file       {priors_file}        \n'
        #     f'    --snp-specific-param-file                     {params_file}        \n'
        # # f'    --special-snps                                {special_snps_file}  \n'
        #     '')

        cmd = (
            ''
            f'{executable}                                                           \n'
            f'    --analysis-files-path                         {lib_dir}            \n'
            f'    --genotyping-node:brlmmp-clustertype          2                    \n'
            f'    --genotyping-node:brlmmp-MS                   0.15                 \n'
            f'    --genotyping-node:brlmmp-SB                   0.75                 \n'
            f'    --genotyping-node:brlmmp-CSepPen              0.1                  \n'
            f'    --genotyping-node:brlmmp-CSepThr              4                    \n'
            f'    --genotyping-node:brlmmp-copytype             -1                   \n'
            f'    --genotyping-node:brlmmp-ocean                0.00001              \n'
            f'    --out-dir                                     {output_dir}         \n'
            f'    --log-file                                    {log_file}           \n'
            f'    --snp-posteriors-output                       true                 \n'
            f'    --genotyping-node:snp-posteriors-output-file  {posteriors_file}    \n'
            f'    --summary-input-file                          {summary_file}       \n'
            f'    --artifact-reduction-trustcheck-file          {trustcheck_file}    \n'
            f'    --read-genders                                {gender_file}        \n'
            f'    --genotyping-node:snp-priors-input-file       {priors_file}        \n'
            f'    --snp-specific-param-file                     {params_file}        \n'
            f'    --special-snps                                {special_snps_file}  \n'
            '')

        # if rare_het_adjustment:
        #     cmd += "    --do-rare-het-adjustment                       true\n"
        #     cmd += "    --genotyping-node:@low-sdDist-minor-cutoff-4   -10\n"
        #     cmd += "    --genotyping-node:@high-sdDist-minor-cutoff-4  5000\n"
        #     cmd += "    --genotyping-node:@low-sdDist-major-cutoff-4   -100\n"
        #     cmd += "    --genotyping-node:@high-sdDist-major-cutoff-4  5000\n"

        if use_copynumber_call_codes is not None:
            cmd += f"    --use-copynumber-call-codes  {use_copynumber_call_codes} \n"
        if probeset_ids_file:
            cmd += f'    --probeset-ids {probeset_ids_file} \n'
        cmd += '\n'

        return Command(cmd, check_log=True)

    def otv_caller(self, genotype_dir, probesets_file):
        executable = self.cmd_paths['otv-caller']
        cmd = (
            ''
            f"{executable}          \n"
            f"    --pid-file        {probesets_file}\n"
            f"    --posterior-file  {genotype_dir / 'AxiomGT1.snp-posteriors.txt'}\n"
            f"    --call-file       {genotype_dir / 'AxiomGT1.calls.txt'}\n"
            f"    --confidence-file {genotype_dir / 'AxiomGT1.confidences.txt'}\n"
            f"    --summary-file    {genotype_dir / 'AxiomGT1.summary.txt'}\n"
            f"    --log-file        {genotype_dir / 'OTV'/ 'otv-caller.log'}\n"
            f"    --output-dir      {genotype_dir / 'OTV'}\n"
            '')
        #       cmd+=" --output-otv-only=false"

        return Command(cmd, check_log=True)

    def ps_metrics(
        self,
        report_file,
        calls_file,
        summary_file,
        posteriors_file,
        output_dir,
        special_snps_file,
        ps_metrics_thresholds,
        genotype_freq_file=None,
        multi_posteriors_file=None,
    ):

        log = output_dir / 'ps-metrics.log'

        executable = self.cmd_paths['ps-metrics']

        cmd = (''
               f"{executable}           \n"
               f"    --summary-file     {summary_file}\n"
               f"    --report-file      {report_file}\n"
               f"    --call-file        {calls_file}\n"
               f"    --output-dir       {output_dir}\n"
               f"    --log-file         {log}\n"
               f"    --posterior-file   {posteriors_file}\n"
               f"    --special-snps     {special_snps_file}\n"
               '')
        if genotype_freq_file and genotype_freq_file.exists():
            cmd += f"    --genotype-freq-file  {genotype_freq_file}\n"
        if multi_posteriors_file and multi_posteriors_file.exists():
            cmd += f"    --multi-posterior-file  {multi_posteriors_file}\n"

        cmd += ps_metrics_thresholds

        return Command(cmd, check_log=True)

    def ps_classification(
        self,
        metrics_file,
        output_dir,
        psct_file,
        ps2snp_file,
        ps_classification_thresholds,
        multi_metrics_file=None,
    ):

        log = output_dir / 'ps-classification.log'
        executable = self.cmd_paths['ps-classification']

        cmd = (''
               f"{executable}           \n"
               f"    --metrics-file     {metrics_file}\n"
               f"    --output-dir       {output_dir}\n"
               f"    --log-file         {log}\n"
               f"    --ps2snp-file      {ps2snp_file}\n"
               '')

        cmd += ps_classification_thresholds

        if multi_metrics_file and multi_metrics_file.exists():
            cmd += f"    --multi-metrics-file    {multi_metrics_file}\n"
        if psct_file and psct_file.exists():
            cmd += f"    --psct-file             {psct_file}\n"

        return Command(cmd, check_log=True)

    def apt_format_result_cnv_vcf(
        self,
        cnvhmm_a5_file: Path,
        annotation_file: Path,
        output_path,
        export_vcf_file,
    ):
        log_file = output_path / 'apt-format-result.log'
        executable = self.cmd_paths['apt-format-result']
        cmd = (''
               f"{executable}\n"
               f"    --cn-region-calls-file     {cnvhmm_a5_file}\n"
               f"    --annotation-file          {annotation_file}\n"
               f"    --export-chr-shortname     True \n"
               f"    --export-vcf-file          {export_vcf_file}\n"
               f"    --log-file                 {log_file} \n"
               '')

        return Command(cmd, check_log=True)

    def apt_format_result_vcf(
        self,
        output_path,
        calls_file,
        performance_file,
        annotation_file,
        export_vcf_file,
    ):
        executable = self.cmd_paths['apt-format-result']
        cmd = (
            ''
            f"{executable}\n"
            f"    --log-file             {output_path / 'apt-format-result-vcf.log'}\n"
            f"    --calls-file           {calls_file}\n"
            f"    --performance-file     {performance_file}\n"
            f"    --annotation-file      {annotation_file}\n"
            f"    --export-chr-shortname True \n"
            f"    --export-confidence    True \n"
            f"    --export-vcf-file      {export_vcf_file}\n"
            '')

        return Command(cmd, check_log=True)

    def apt_format_result_plink(
        self,
        output_path,
        calls_file,
        annotation_file,
        plink_file,
        pedigree_file,
        snp_list_file=None,
    ):
        executable = self.cmd_paths['apt-format-result']

        cmd = (
            ''
            f"{executable}\n"
            f"    --log-file             {output_path / 'apt-format-result-plink.log'}\n"
            f"    --calls-file           {calls_file}\n"
            f"    --annotation-file      {annotation_file}\n"
            f"    --export-chr-shortname True \n"
            f"    --snp-identifier-column Affy_SNP_ID\n"
            f"    --export-plinkt-file    {plink_file}\n"
            f"    --pedigree-file        {pedigree_file}\n"
            '')
        if snp_list_file:

            cmd += f"    --snp-list-file        {snp_list_file}\n"

        return Command(cmd, check_log=True)


class Command():

    def __init__(self, cmd, check_log=False):
        self._cmd = cmd
        self.check_log = check_log

    def __str__(self):
        return self._cmd

    def execute(self, fatal=True):
        cmd = self._cmd

        cmd = cmd.replace('\n', '')

        msg_buffer = []
        n_errors = None

        with Popen(
                cmd,
                shell=True,
                universal_newlines=True,
                stdout=PIPE,
                stderr=STDOUT,
        ) as proc:
            for line in proc.stdout:
                line = line.strip()
                msg_buffer.append(line)
                if not self.check_log:
                    continue
                match = config.APT_STATUS_LINE.match(line)
                if match:
                    n_errors = int(match.group(1))

            proc.wait()

        msg = '\n'.join(msg_buffer)

        exit_status = proc.returncode

        if exit_status == 0 and not self.check_log:
            return exit_status, msg
        elif exit_status == 0 and n_errors == 0:
            return exit_status, msg
        else:
            print(f'exit_status: {exit_status}, n_errors: {n_errors}')
            assert False

            # raise Exception(f'Error: {cmd}')

    def _apt_crashed(self, apt_log_file):

        with open(apt_log_file, 'rt') as fd:
            for line in fd:
                line = line.strip()
                match = config.APT_STATUS_LINE.match(line)
                if match:
                    n_errors = int(match.group(1))

                    if n_errors == 0:
                        return False
                    else:
                        return True
        return True
