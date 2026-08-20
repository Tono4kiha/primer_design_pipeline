# Primer enumeration / single-primer filters
PRIMER_MIN_LEN = 18
PRIMER_MAX_LEN = 27
PRIMER_MIN_TM = 57.0
PRIMER_OPT_TM = 60.0
PRIMER_MAX_TM = 63.0
PRIMER_MIN_GC = 40.0
PRIMER_MAX_GC = 60.0
MAX_HOMOPOLYMER = 4
MIN_SHANNON_ENTROPY = 1.50
MAX_TANDEM_COPIES = 4

# Primer3 thermodynamic cutoffs (deg C)
MAX_HAIRPIN_TM = 47.0
MAX_HOMODIMER_TM = 47.0
MAX_SELF_END_TM = 47.0

# PCR chemistry passed to primer3-py
# Change these to match your actual PCR conditions.
MV_CONC = 50.0      # mM monovalent cation
DV_CONC = 1.5       # mM divalent cation
DNTP_CONC = 0.6     # mM
DNA_CONC = 50.0     # nM oligo concentration used by Primer3 model
THERMO_TEMP_C = 37.0
MAX_LOOP = 30

# Pair filters
PRODUCT_MIN = 50
PRODUCT_OPT = 250
PRODUCT_MAX = 700
MAX_DELTA_TM = 2.5
MAX_HETERODIMER_TM = 47.0
MAX_PAIR_END_TM = 47.0

# BLAST hit acceptance for possible primer-template binding
BLAST_MIN_QUERY_COVERAGE = 0.80
BLAST_MAX_MISMATCH = 3
BLAST_MAX_GAPOPEN = 0
BLAST_3P_WINDOW = 5
BLAST_MAX_3P_MISMATCH = 1

# When reconstructing possible off-target PCR products from primer BLAST hits,
# inspect a wider window than the intended 50-700 bp product range.
OFFTARGET_PRODUCT_MIN = 40
OFFTARGET_PRODUCT_MAX = 2000
