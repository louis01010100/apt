import logging
from pathlib import Path

import pandas as pd

from . import config, utils
from .arguments import SampleQcArguments, SnvArguments
from .commands import Commands
from .library import Library


def dqc(
    samples: pd.DataFrame,
    lib_dir: Path,
    output_dir: Path,
    commands: Commands,
    sqc_args: SampleQcArguments,
    force: bool,
):

    failed = dict()

    output_dir.mkdir(parents=True, exist_ok=True)

    cels_file = output_dir / "dqc_cels.txt"

    utils.export_cels(samples, cels_file)

    cmd = commands.apt_geno_qc(
        lib_dir=lib_dir,
        args_file=sqc_args.geno_qc_args_file,
        cel_files=cels_file,
        output_dir=output_dir,
        force=force,
    )

    cmd.execute()

    outputFile = output_dir / config.GENO_QC_FILENAME

    dqc_report = pd.read_csv(
        outputFile,
        comment='#',
        header=0,
        sep='\t',
        dtype={
            'cel_files': 'str',
            'axiom_dishqc_DQC': 'float64'
        },
        usecols=['cel_files', 'axiom_dishqc_DQC'],
    )

    dqc_report = dqc_report.rename(columns={
        'cel_files': 'cel_name',
        'axiom_dishqc_DQC': 'dqc',
    })

    dqc_report = dqc_report.assign(
        passing_dqc=lambda x: x['dqc'] >= sqc_args.dqc_threshold)

    failed_samples = dqc_report[lambda x: x['passing_dqc'] == False]

    for _, sample in failed_samples.iterrows():
        cel_name = sample['cel_name']
        dqc = sample['dqc']
        logging.info(f'DQC\t{cel_name}\t{dqc:5.2f}\tfailed')

    return dqc_report


def qccr(
    samples: pd.DataFrame,
    lib_dir: Path,
    sqc_args: SampleQcArguments,
    output_dir: Path,
    commands: Commands,
    force: bool,
):

    output_dir.mkdir(parents=True, exist_ok=True)

    cels_file = output_dir / "qccr_cels.txt"

    utils.export_cels(samples, cels_file)

    cmd = commands.apt_genotype_axiom(
        cel_files=cels_file,
        output_dir=output_dir,
        lib_dir=lib_dir,
        args_file=sqc_args.step1_args_file,
        snp_priors_file=sqc_args.snp_priors_file,
        force=force,
    )

    cmd.execute()

    qccr_report = utils.call_rate_from_calls_file(
        output_dir / config.CALLS_FILENAME, )

    qccr_report = qccr_report.rename(columns={'call_rate': 'qccr'})
    qccr_report = qccr_report.round(10)

    qccr_report = qccr_report.assign(
        passing_qccr=lambda x: x['qccr'] >= sqc_args.qccr_threshold)

    failed = qccr_report[lambda x: x['passing_qccr'] == False]

    for _, record in failed.iterrows():
        cel_name = record['cel_name']
        qccr = record['qccr']
        logging.info(f'QCCR\t{cel_name}\t{qccr:5.2f}\tfailed')

    return qccr_report


def plate_qc(
    samples: pd.DataFrame,
    output_dir: Path,
    avg_qccr_threshold: float,
    min_samples_for_plate_qccr: float,
):

    output_dir.mkdir(parents=True, exist_ok=True)

    plate_dqc_report = _report_plate_dqc(samples)
    plate_qccr_report = _report_plate_qccr(
        samples,
        avg_qccr_threshold,
        min_samples_for_plate_qccr,
    )

    return plate_dqc_report, plate_qccr_report


def genotype_summary(
    probeset_ids_file: Path,
    summary_file: Path,
    trustcheck_file: Path,
    gender_file: Path,
    output_dir: Path,
    snv_args: SnvArguments,
    commands: Commands,
    skip_qc: bool = False,
):

    commands.apt_summary_genotype_axiom(
        probeset_ids_file=probeset_ids_file,
        args_file=snv_args.step2_args_file,
        summary_file=summary_file,
        trustcheck_file=trustcheck_file,
        priors_file=snv_args.snp_priors_file,
        params_file=snv_args.snp_params_file,
        special_snps_file=snv_args.special_snps_file,
        gender_file=gender_file,
        lib_dir=snv_args.lib_dir,
        output_dir=output_dir,
        use_copynumber_call_codes=snv_args.copynumber_call_codes,
        rare_het_adjustment=snv_args.rare_het_adjustment,
    ).execute()

    if skip_qc:
        return
    snv_qc(
        summary_file=summary_file,
        report_file=output_dir / config.REPORT_FILENAME,
        calls_file=output_dir / config.CALLS_FILENAME,
        posteriors_file=output_dir / config.POSTERIORS_FILENAME,
        snv_args=snv_args,
        multi_posteriors_file=output_dir / config.MULTI_POSTERIORS_FILENAME,
        output_dir=output_dir,
        commands=commands,
    )


def genotype(
    cels_file: Path,
    snv_args: SnvArguments,
    output_dir: Path,
    cnpscalls_file: Path,
    commands: Commands,
    force: bool,
    skip_qc: bool = False,
):

    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = commands.apt_genotype_axiom(
        cel_files=cels_file,
        output_dir=output_dir,
        lib_dir=snv_args.lib_dir,
        args_file=snv_args.step2_args_file,
        snp_priors_file=snv_args.snp_priors_file,
        snp_params_file=snv_args.snp_params_file,
        probeset_ids_file=snv_args.probeset_ids_file,
        rare_het_adjustment=snv_args.rare_het_adjustment,
        force=force,
        process_multi_alleles=snv_args.multi_alleles,
        use_copynumber_call_codes=snv_args.copynumber_call_codes,
        cnpscalls_file=cnpscalls_file,
        export_allele_summaries=True,
        dump_blemishes=True,
        export_trustcheck=True,
    )

    cmd.execute()

    if skip_qc:
        return

    snv_qc(
        summary_file=output_dir / config.SUMMARY_FILENAME,
        report_file=output_dir / config.REPORT_FILENAME,
        calls_file=output_dir / config.CALLS_FILENAME,
        posteriors_file=output_dir / config.POSTERIORS_FILENAME,
        snv_args=snv_args,
        multi_posteriors_file=output_dir / config.MULTI_POSTERIORS_FILENAME,
        output_dir=output_dir,
        commands=commands,
    )


def snv_qc(
    summary_file: Path,
    report_file: Path,
    calls_file: Path,
    posteriors_file: Path,
    snv_args: SnvArguments,
    output_dir: Path,
    commands: Commands,
    multi_posteriors_file: Path = None,
):

    output_dir.mkdir(parents=True, exist_ok=True)

    ps_metrics_cmd = commands.ps_metrics(
        report_file=report_file,
        calls_file=calls_file,
        summary_file=summary_file,
        posteriors_file=posteriors_file,
        multi_posteriors_file=multi_posteriors_file,
        output_dir=output_dir,
        ps_metrics_thresholds=snv_args.ps_metrics_thresholds,
        genotype_freq_file=snv_args.genotype_freq_file,
        special_snps_file=snv_args.special_snps_file,
    )

    ps_metrics_cmd.execute()

    multi_metrics_file = output_dir / config.MULTI_METRICS_FILENAME

    if not multi_metrics_file.exists():
        multi_metrics_file = None

    ps_classification_cmd = commands.ps_classification(
        metrics_file=output_dir / config.METRICS_FILENAME,
        output_dir=output_dir,
        psct_file=snv_args.psct_file,
        ps2snp_file=snv_args.ps2snp_file,
        multi_metrics_file=multi_metrics_file,
        ps_classification_thresholds=snv_args.ps_classification_thresholds,
    )

    ps_classification_cmd.execute()


def export_signals(
    cels_file: Path,
    args_file: Path,
    lib_dir: Path,
    output_dir: Path,
    commands: Commands,
    force: bool,
):

    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = commands.export_signals(
        cels_file=cels_file,
        output_dir=output_dir,
        lib_dir=lib_dir,
        force=force,
        args_file=args_file,
    )

    cmd.execute()


def _report_plate_dqc(samples: pd.DataFrame):

    bag = []

    for plate, df in samples.groupby('sbarcode', sort=True):

        n_samples = len(df)
        n_samples_dqc_ok = len(df[lambda x: x['passing_dqc'] == True])
        n_samples_dqc_not_ok = n_samples - n_samples_dqc_ok
        dqc_pass_rate = n_samples_dqc_ok / n_samples * 100
        avg_dqc = df['dqc'].mean()

        bag.append({
            'sbarcode': plate,
            'n_samples_dqc_ok': n_samples_dqc_ok,
            'dqc_pass_rate': dqc_pass_rate,
            'avg_dqc': avg_dqc,
        })

    return pd.DataFrame.from_records(bag)


def _report_plate_qccr(
    samples: pd.DataFrame,
    avg_qccr_threshold: float,
    min_samples_for_plate_qccr: int,
):

    bag = []

    for sbarcode, df in samples.groupby('sbarcode', sort=True):

        n_samples = len(samples)

        n_samples_dqc_ok = len(df[lambda x: x['passing_dqc'] == True])
        samples_qccr_ok = df[lambda x: x['passing_qccr'] == True]
        n_samples_qccr_ok = len(samples_qccr_ok)

        avg_qccr = samples_qccr_ok['qccr'].mean()

        n_samples = len(samples)
        plate_pass_rate = n_samples_qccr_ok / n_samples * 100
        qccr_pass_rate = n_samples_qccr_ok / n_samples_dqc_ok * 100

        bag.append({
            'sbarcode': sbarcode,
            'n_samples_qccr_ok': n_samples_qccr_ok,
            'qccr_pass_rate': qccr_pass_rate,
            'plate_pass_rate': plate_pass_rate,
            'avg_qccr': avg_qccr,
        })

        if (n_samples_qccr_ok >= min_samples_for_plate_qccr
                and avg_qccr < avg_qccr_threshold):
            logging.warn(
                f"QCCR_Plate\t{sbarcode}\tsamples {n_samples}\t{avg_qccr:5.2f}\tfailed"
            )

    report = pd.DataFrame.from_records(bag)

    report = report.assign(
        passing_avg_qccr=lambda x: x['avg_qccr'] >= avg_qccr_threshold)

    return report
