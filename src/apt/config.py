import re
from importlib import resources

BIN_DIRS = [
    '/opt/apt/bin',
]

APT_PROGRAMS = [
    'apt-package-util',
    'apt-dmet-translation',
    'apt-geno-qc',
    'apt-genotype-axiom',
    'apt-summary-genotype-axiom',
    'ps-metrics',
    'ps-classification',
    'otv-caller',
    'apt-format-result',
    'apt2-dset-util',
    'apt-copynumber-axiom-cnvmix',
    'apt-copynumber-axiom-ref',
    'apt-copynumber-axiom-hmm',
]

BANNER_WIDTH = 50

RARE_HET_REPORT_FILENAME = 'AxiomGT1.rare_het.report.txt'
FLOATING_POINT_PRECISION = 5
CN_REGION_CALLS_FILENAME = 'AxiomCNVMix.cnregioncalls.txt'
ALLELE_TRANSLATION_DIRNAME = 'allele_translation'
AXAS_DIRNAME = 'axas'
GENO_QC_FILENAME = 'apt-geno-qc.txt'
CALLS_FILENAME = 'AxiomGT1.calls.txt'
REPORT_FILENAME = 'AxiomGT1.report.txt'
PERFORMANCE_FILENAME = 'Ps.performance.txt'
SUMMARY_FILENAME = 'AxiomGT1.summary.txt'
SUMMARY_A5_FILENAME = 'AxiomGT1.summary.a5'
POSTERIORS_FILENAME = 'AxiomGT1.snp-posteriors.txt'
MULTI_POSTERIORS_FILENAME = 'AxiomGT1.snp-posteriors.multi.txt'
METRICS_FILENAME = 'metrics.txt'
MULTI_METRICS_FILENAME = 'multi-metrics.txt'
CONFIDENCES_FILENAME = 'AxiomGT1.confidences.txt'
TRUSTCHECK_FILENAME = 'AxiomGT1.trustcheck.txt'
RECOMMENDED_FILENAME = 'Recommended.ps'
CNPSCALLS_FILENAME = 'AxiomCNVMix.cnpscalls.txt'
HMM_CNV_FILENAME = 'AxiomHMM.cnv.a5'

PLATEMAP_FILENAME = 'platemap.tsv'
GENDER_FILENAME = 'gender.tsv'
SAMPLE_ORDER_FILENAME = 'sample_order.tsv'

PEDIGREE_FILENAME = 'pedigree.tsv'
SUBOPTIMAL_PROBESETS_FILENAME = 'suboptimal.ps'
RECALL_PROBESETS_FILENAME = 'recall.ps'
SIBLING_PROBESETS_FILENAME = 'sibling.ps'
TARGET_PROBESETS_FILENAME = 'target.tsv'
UGF_PROBESETS_FILENAME = 'ugf.ps'
OUTLIER_SUPPORTING_FILENAME = 'outlier.tsv'

UNEXPECTED_GENOTYPE_FREQ = "UnexpectedGenotypeFreq"
GENOTYPE_P_VALUE_CUTOFF_FOR_PNORM = 0.5

CONVERSION_TYPES = [
    "PolyHighResolution",
    "NoMinorHom",
    "MonoHighResolution",
    "CallRateBelowThreshold",
    "OffTargetVariant",
    UNEXPECTED_GENOTYPE_FREQ,
    "Other",
    "OTV",
]

CONVERSION_TYPES_RECOMMENDED = [
    "PolyHighResolution",
    "NoMinorHom",
    "MonoHighResolution",
]

CONVERSION_TYPES_PNORM = [
    "CallRateBelowThreshold",
    "OffTargetVariant",
    "Other",
    UNEXPECTED_GENOTYPE_FREQ,
]

#`#    info     1 816 | 0 error(s) and 670791 warning(s).'
APT_STATUS_LINE = re.compile(
    r'^#.+info.+\| ([0-9]+) error\(s\) and (?:[0-9]+) warning\(s\)\.$')

## plate normalization

RESOURCES_DIR = resources.files('apt') / 'resources'

FIND_OUTLIERS_SCRIPT = RESOURCES_DIR / 'find_outliers.r'
PLATE_NORM_REDUCE_SCRIPT = RESOURCES_DIR / 'platenorm_reduce.r'

PLATEMAP_FILENAME = 'platemap.tsv'
GENDER_FILENAME = 'gender.tsv'
SAMPLE_ORDER_FILENAME = 'sample_order.tsv'

ESSENTIAL_FILES = [
    SUMMARY_FILENAME,
    REPORT_FILENAME,
    CALLS_FILENAME,
    CONFIDENCES_FILENAME,
    POSTERIORS_FILENAME,
    METRICS_FILENAME,
    PERFORMANCE_FILENAME,
]

PNORM_DEFAULT_DIRNAME = 'default'
PNORM_NORM_DIRNAME = 'norm'
IMPROVED_PROBESETS_FILE = "improved.ps"

PLOT2MIN_DQC = 0.80
PLOT4MIN_DQC = 0.82
PLOT4MIN_QCCR = 90
MIN_MEAN_QCCR = 97
plot7MinQCCR = 98.5
plot7MinRecommendedRate = 81

SNPOLISHER_COUNT = 200

SIGNALS_DIRNAME = 'signals'
SAMPLE_QC_REPORT_FILENAME = 'qc.tsv'
CNQC_REPORT_FILENAME = 'cnref_qc.tsv'
CNREF_FILENAME = 'self.cn_models'
CNV_REPORT_FILENAME = 'cnv.tsv'
LOH_REPORT_FILENAME = 'loh.tsv'
TMP_DIRNAME = 'tmp'
CNVHMM_A5_FILENAME = 'AxiomHMM.cnv.a5'
VCF_DIRNAME = 'vcf'
AXAS_DIRNAME = 'axas'
CNV_DIRNAME = 'cnv'
IGV_DIRNAME = 'igv'
SNV_DIRNAME = 'snv'
SAMPLE_QC_DIRNAME = 'qc'
TRANSLATION_DIRNAME = 'tx'
PLINK_DIRNAME = 'plink'
REPORT_DIRNAME = 'report'

EXPECTED_GMMS_FILENAME = 'exp_gmms.tsv'
GMM_SIMILARITIES_FILENAME = 'gmm_similarities.tsv'
BATCH_ID = 'batch_id'
GT_DIR = 'snv_dir'

GMM_COSINE_CUTOFF = 0.955

CONFIDENCE_SCORE_THRESHOLD = 0.15
N_JOBS = 1000

OCEAN_PROBABILITY = 0.00001

CEL_ORDER_PATTERN = re.compile(
    '^#%affymetrix-algorithm-param-apt-opt-cel-([0-9]+)=(.+\.CEL)$')
