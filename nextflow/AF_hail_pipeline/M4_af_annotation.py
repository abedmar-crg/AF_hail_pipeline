import hail as hl
import logging
import csv
import gzip
import json
import os
import re
from pathlib import Path
from utils import *

config = load_config()

def stats_by_sex(mt):
    """
    Creates stats: number of female, number of male and undefined 
    :params: mt 
    :return: stats about sexes in data 
    """
   
    results_sex_agg = mt.aggregate_cols(hl.agg.counter(mt.imputed_sex))

    num_females = results_sex_agg.get("XX", 0)
    num_males = results_sex_agg.get("XY", 0)
    num_missing = results_sex_agg.get("undefined", 0)

    logging.info(f"Sex Statistics")
    logging.info(f"Number of females: {num_females}")
    logging.info(f"Number of males: {num_males}")
    logging.info(f"Number of samples without sex information: {num_missing}")
    sex_summary = []
    sex_summary.append("")
    sex_summary.append(["Sex Statistics"])
    sex_summary.append(["Number of females", num_females])
    sex_summary.append(["Number of males", num_males])
    sex_summary.append(["Number of samples without sex information:", num_missing])
    csv_writer(sex_summary)

    return mt, results_sex_agg


def stats_by_ancestry(mt):   
    """
    Creates stats: number of samples per population in the data
    :params: mt 
    :return: stats about ancestry groups in the data 
    """
    results_ancestry_agg = mt.aggregate_cols(hl.agg.counter(mt.ancestry))  # agg samples (columns) by ancestry

    if None in results_ancestry_agg:
        logging.error("Some samples without ancestry information")

    logging.info("Ancestry Statistics")
    logging.info(f"Ancestry super-populations: {results_ancestry_agg}") 
    anc_summary = []
    anc_summary.append("")
    anc_summary.append(["Ancestry statistics"])
    for ancestry, count in results_ancestry_agg.items():
        anc_summary.append([ancestry, count])
    csv_writer(anc_summary)
    return mt, results_ancestry_agg

def af_by_sex_ancestry(mt, results_sex_agg, results_ancestry_agg):
    """
    Recalculates AF fields grouped by sex and ancestry 
    :params:
    - mt: Hail MatrixTable containing genotype data
    - results_sex_agg: Dictionary with sex categories present in the dataset
    - results_ancestry_agg: Dictionary with ancestry groups present in the dataset
    
    :return:
    - mt: Updated MatrixTable with calculated statistics
    - AF_total: Overall allele frequency for alternate allele
    - AC_total: Overall allele count for alternate allele  
    - AN_total: Overall allele number (total chromosomes analyzed)
    - nhom_alt_total: Overall count of homozygous alternate individuals
    - AF_sex: Dictionary of sex-stratified frequency statistics
    - AF_ancestry: Dictionary of ancestry-stratified frequency statistics
    """
    logging.info("Re-calculating AFs by sex and ancestry if available")

    has_sex = 'imputed_sex' in mt.col.dtype.keys()

    # PRE-STEP : Ensure haploid calls are correctly written 

    if has_sex:
        mt = mt.annotate_entries(
            GT = hl.case()
            .when(
                # Condition: Male sample AND non-PAR region AND diploid homozygous variant call
                (mt.imputed_sex == "XY") & 
                (~mt.locus.in_autosome_or_par()) & 
                mt.GT.is_diploid() & 
                (mt.GT.is_hom_ref() | mt.GT.is_hom_var()),
                # Convert to haploid: take first allele only
                hl.call(mt.GT[0])
            )
            .default(mt.GT)  # Keep original GT for all other cases
        )


    # STEP 1: Calculate comprehensive genotype statistics
    row_annotations = {
        'gt_stats': hl.agg.call_stats(mt.GT, mt.alleles)
    }
    
    # STEP 2: Conditionally calculate sex-stratified statistics
    if has_sex:
        row_annotations['gt_stats_by_sex'] = hl.agg.group_by(mt.imputed_sex, hl.agg.call_stats(mt.GT, mt.alleles))
    
    # STEP 3: Conditionally calculate ancestry-stratified statistics
    # Only compute ancestry stats if the ancestry field exists in the dataset
    if 'ancestry' in mt.col:
        row_annotations['gt_stats_by_ancestry'] = hl.agg.group_by(mt.ancestry, hl.agg.call_stats(mt.GT, mt.alleles))
        # Only do nested grouping if BOTH exist
        if has_sex:
            row_annotations['gt_stats_by_ancestry_sex'] = hl.agg.group_by(
                hl.struct(ancestry=mt.ancestry, sex=mt.imputed_sex),
                hl.agg.call_stats(mt.GT, mt.alleles)
            )

    mt = mt.annotate_rows(**row_annotations)

    # STEP 4: Extract overall population statistics
    # Index [1] refers to the alternate allele (index [0] would be reference allele)
    AF_total = mt.gt_stats.AF[1] # Alternate allele frequency
    AC_total = mt.gt_stats.AC[1] # Alternate allele count
    AN_total = mt.gt_stats.AN # Total allele number (2 * sample count for diploid)
    nhom_total = mt.gt_stats.homozygote_count[1] # Homozygous alternate count

    # HEMIZYGOTES (Total) - Default to 0 if no sex info
    if has_sex:
        nhemi_total = hl.if_else(~mt.locus.in_autosome_or_par(),
            hl.or_else(mt.gt_stats_by_sex.get("XY").AC[1], 0), 0)
    else:
        nhemi_total = hl.int32(0)

    nhet_total = AC_total - 2*nhom_total - nhemi_total # heterozygous count

    # STEP 5: Extract sex-stratified statistics with hemizygote handling
    AF_sex = {}
    if has_sex:
        sexes_avail = list(results_sex_agg.keys())  # Get available sex categories from input
        logging.warning("Samples where sex couldn't be inferred will be ignored")
        
        for sex in sexes_avail: # create stats by sex
            
            if sex == "undefined":
                continue  # Skip samples with unknown sex
            
            AF_sex[f"AF_{sex}_recalc"] = mt.gt_stats_by_sex[sex].AF[1]
            AF_sex[f"AC_{sex}_recalc"] = mt.gt_stats_by_sex[sex].AC[1]
            AF_sex[f"nhom_{sex}_recalc"] = mt.gt_stats_by_sex[sex].homozygote_count[1]
            AF_sex[f"AN_{sex}_recalc"] = mt.gt_stats_by_sex[sex].AN

            # HEMIZYGOTE HANDLING FOR SEX CHROMOSOMES
            # Add hemizygote allele count if variant is on non-autosomal/non-PAR regions
            # Only XY (male) individuals can be hemizygous for X chromosome variants
            AF_sex[f"nhemi_{sex}_recalc"] = hl.if_else(
                # Check if variant is on sex chromosome outside pseudoautosomal regions
                ~mt.locus.in_autosome_or_par(),
                # If on sex chromosome and this is XY group, use the allele count as hemizygote count
                hl.if_else(
                    sex == "XY", 
                    mt.gt_stats_by_sex[sex].AC[1],  # Males: AC equals hemizygote count
                    hl.int32(0)  # Females: no hemizygotes possible
                ),
                hl.int32(0)  # Autosomal variants: no hemizygotes
            )

            AF_sex[f"nhet_{sex}_recalc"] = AF_sex[f"AC_{sex}_recalc"] - 2*(AF_sex[f"nhom_{sex}_recalc"]) - AF_sex[f"nhemi_{sex}_recalc"] # heterozygous count



    # STEP 6: Extract ancestry-stratified statistics (if available)
    if results_ancestry_agg != "":
        AF_ancestry = {}
        ancestry_avail = list(results_ancestry_agg.keys()) # gets the ancestries in the VCF
        
        for ancestry in ancestry_avail: # create stats by ancestry
            
            if ancestry == "Multi-ancestry":
                continue  # Skip samples with multi-ancestry
            
            AF_ancestry[f"AF_{ancestry}_recalc"] = mt.gt_stats_by_ancestry[ancestry].AF[1]
            AF_ancestry[f"AC_{ancestry}_recalc"] = mt.gt_stats_by_ancestry[ancestry].AC[1]
            AF_ancestry[f"nhom_{ancestry}_recalc"] = mt.gt_stats_by_ancestry[ancestry].homozygote_count[1]
            AF_ancestry[f"AN_{ancestry}_recalc"] = mt.gt_stats_by_ancestry[ancestry].AN

            # HEMIZYGOTE HANDLING FOR ANCESTRY GROUPS
            # For ancestry groups, we need to extract the XY subset within each ancestry
            # This requires a nested grouping: ancestry -> sex -> stats
            # Create the key for this ancestry + XY combination
            if has_sex:
                anc_xy_key = hl.struct(ancestry=ancestry, sex="XY")
                AF_ancestry[f"nhemi_{ancestry}_recalc"] = hl.if_else(
                    ~mt.locus.in_autosome_or_par(),
                    hl.or_else(mt.gt_stats_by_ancestry_sex.get(anc_xy_key).AC[1], 0), 0)
            else:
                AF_ancestry[f"nhemi_{ancestry}_recalc"] = hl.int32(0)

            AF_ancestry[f"nhet_{ancestry}_recalc"] = AF_ancestry[f"AC_{ancestry}_recalc"] - 2*AF_ancestry[f"nhom_{ancestry}_recalc"] - AF_ancestry[f"nhemi_{ancestry}_recalc"]
    else:
        AF_ancestry={}

        # STEP 7: Extract ancestry-by-sex stratified statistics (if both are available)
    # Produces fields like AF_EUR_XY_recalc, AF_EUR_XX_recalc, etc.
    AF_ancestry_sex = {}

    if results_ancestry_agg != "" and has_sex:
        
        for ancestry in ancestry_avail:

            if ancestry == "Multi-ancestry":
                continue

            for sex in sexes_avail:

                if sex == "undefined":
                    continue

                key = hl.struct(ancestry=ancestry, sex=sex)

                AF_ancestry_sex[f"AF_{ancestry}_{sex}_recalc"]   = hl.or_else(mt.gt_stats_by_ancestry_sex.get(key).AF[1],   hl.float64(0))
                AF_ancestry_sex[f"AC_{ancestry}_{sex}_recalc"]   = hl.or_else(mt.gt_stats_by_ancestry_sex.get(key).AC[1],   hl.int32(0))
                AF_ancestry_sex[f"AN_{ancestry}_{sex}_recalc"]   = hl.or_else(mt.gt_stats_by_ancestry_sex.get(key).AN,      hl.int32(0))
                AF_ancestry_sex[f"nhom_{ancestry}_{sex}_recalc"] = hl.or_else(mt.gt_stats_by_ancestry_sex.get(key).homozygote_count[1], hl.int32(0))

                # HEMIZYGOTE HANDLING FOR ANCESTRY x SEX
                # Only XY individuals outside autosomal/PAR regions can be hemizygous
                AF_ancestry_sex[f"nhemi_{ancestry}_{sex}_recalc"] = hl.if_else(
                    ~mt.locus.in_autosome_or_par(),
                    hl.if_else(
                        sex == "XY",
                        hl.or_else(mt.gt_stats_by_ancestry_sex.get(key).AC[1], hl.int32(0)),
                        hl.int32(0)  # XX individuals: no hemizygotes
                    ),
                    hl.int32(0)  # Autosomal variants: no hemizygotes
                )

                AF_ancestry_sex[f"nhet_{ancestry}_{sex}_recalc"] = (
                    AF_ancestry_sex[f"AC_{ancestry}_{sex}_recalc"]
                    - 2 * AF_ancestry_sex[f"nhom_{ancestry}_{sex}_recalc"]
                    - AF_ancestry_sex[f"nhemi_{ancestry}_{sex}_recalc"]
                )

    return mt, AF_total, AC_total, AN_total, nhom_total, nhet_total, nhemi_total, AF_sex, AF_ancestry, AF_ancestry_sex

def to_allele_counts(fields: dict) -> dict:
    """
    Converts nhom/nhet/nhemi genotype count fields to allele count equivalents.
    AC_hom = 2*nhom, AC_het = nhet, AC_hemi = nhemi
    """
    converted = {}
    for key, val in fields.items():
        if key.startswith("nhom_"):
            converted[key.replace("nhom_", "AC_hom_")] = 2 * val
        elif key.startswith("nhet_"):
            converted[key.replace("nhet_", "AC_het_")] = val
        elif key.startswith("nhemi_"):
            converted[key.replace("nhemi_", "AC_hemi_")] = val
        else:
            converted[key] = val
    return converted

def annotate_new_vcf(mt, AF_total, AC_total, AN_total, homozygote_count_total, 
                     heterozygous_count, hemizygotes_total, AF_sex, AF_ancestry, AF_ancestry_sex):

    if config['use_allele_counts']:
        nhom_out   = 2 * homozygote_count_total   # AC_hom
        nhet_out   = heterozygous_count            # AC_het (nhet == AC_het)
        nhemi_out  = hemizygotes_total             # AC_hemi
        AF_sex           = to_allele_counts(AF_sex)
        AF_ancestry      = to_allele_counts(AF_ancestry)
        AF_ancestry_sex  = to_allele_counts(AF_ancestry_sex)
        hom_key, het_key, hemi_key = "AC_hom_total_recalc", "AC_het_total_recalc", "AC_hemi_total_recalc"
    else:
        nhom_out, nhet_out, nhemi_out = homozygote_count_total, heterozygous_count, hemizygotes_total
        hom_key, het_key, hemi_key = "nhom_alt_total_recalc", "nhet_total_recalc", "nhemi_total_recalc"

    mt_af = mt.annotate_rows(
        info=mt.info.annotate(
            AF_total_recalc=AF_total,
            AC_total_recalc=AC_total,
            AN_total_recalc=AN_total,
            **{hom_key: nhom_out},
            **{het_key: nhet_out},
            **{hemi_key: nhemi_out},
            **AF_sex,
            **AF_ancestry,
            **AF_ancestry_sex
        )
    )
    return mt_af


def export_new_vcf(mt_af, ac0_applied=False):
    """
    Export mt with ancestry and sex AF fields to VCF 
    :params: mt with AF fields by sex and ancestry
    :return: VCF with AF fields by sex and ancestry
    """
    # Get all VCF file paths in the directory
    vcfs = sorted([
        os.path.join(config['vcf_dir'], f)
        for f in os.listdir(config['vcf_dir'])
        if f.endswith('.vcf') or f.endswith('.vcf.bgz') or f.endswith('.vcf.gz')
    ])
    
    header_metadata = hl.get_vcf_metadata(vcfs[0])  # get original header
    if ac0_applied:
        # Hail requires explicit FILTER metadata to document a new VCF flag.
        header_metadata.setdefault("filter", {})["AC0"] = {
            "Description": (
                "Adjusted ALT allele count is zero after genotype QC "
                "(GQ>=20; DP>=10 for diploid or >=5 for haploid; "
                "AB>0.2 for heterozygous calls)"
            )
        }
    
    if config['summary_VCF'] == True:
        logging.info("Exporting final summary VCF with recalculated AF fields and no sample information")
        # Remove all sample columns
        variants_table = mt_af.rows()
        hl.export_vcf(
        variants_table,
        config['final_vcf_AF'],
        metadata=header_metadata,
        tabix=True
    )

    else:
        logging.info("Exporting final VCF with recalculated AF fields and sample information")
        hl.export_vcf(
        mt_af,
        config['final_vcf_AF'],
        metadata=header_metadata,
        tabix=True
        )
    



# ============================================================
# Optional AC0 FILTER flagging (integrated into M4)
# ============================================================

REQUIRED_FORMAT_FIELDS = frozenset(("GT", "GQ", "DP", "AD"))
# gnomAD adjusted-genotype defaults (as of gnomad_methods get_adj_expr).
GQ_MIN = 20
DIPLOID_DP_MIN = 10
HAPLOID_DP_MIN = 5
HET_AB_MIN = 0.2

# Not Hail-specific; unit tests can run without a Spark/Hail installation.
FORMAT_DECLARATION = re.compile(r"(?:^|[,<])ID=([^,>]+)")


def list_input_vcfs(vcf_dir):
    """Use the same accepted extensions as M1 and M4, excluding index files."""
    return sorted(
        str(p) for p in Path(vcf_dir).iterdir()
        if p.is_file() and p.name.endswith((".vcf", ".vcf.bgz", ".vcf.gz"))
    )


def inspect_vcf_header(vcf_path):
    """Read just the header: detects absent FORMAT declarations per input VCF."""
    compressed = str(vcf_path).endswith((".bgz", ".gz"))
    opener = gzip.open if compressed else open
    declared = set()
    samples_present = False
    header_found = False
    with opener(vcf_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            fmt = FORMAT_DECLARATION.search(line) if line.startswith("##FORMAT=<") else None
            if fmt:
                declared.add(fmt.group(1))
            elif line.startswith("#CHROM\t"):
                header_found = True
                # The FORMAT column and at least one actual sample are required.
                cols = line.rstrip("\r\n").split("\t")
                samples_present = len(cols) >= 10 and cols[8] == "FORMAT"
                break
    if not header_found:
        raise ValueError("VCF has no #CHROM header line")
    return {"fields": sorted(declared), "has_sample_columns": samples_present}


def inspect_input_formats(vcf_dir):
    """Return per-file FORMAT evidence and missing-field lists.

    If one input lacks a field, AC0 is skipped for the entire merged dataset,
    rather than risking false AC0 calls on the incomplete input.
    """
    audits = []
    for fn in list_input_vcfs(vcf_dir):
        info = {"file": fn, "missing_fields": [], "error": None}
        try:
            hdr = inspect_vcf_header(fn)
            info["present_fields"] = hdr["fields"]
            info["has_sample_columns"] = hdr["has_sample_columns"]
            info["missing_fields"] = sorted(REQUIRED_FORMAT_FIELDS - set(hdr["fields"]))
            if not hdr["has_sample_columns"]:
                info["missing_fields"] = sorted(set(info["missing_fields"]) | {"SAMPLE_COLUMNS"})
        except (OSError, EOFError, UnicodeError, ValueError) as exc:
            info["error"] = str(exc)
            info["present_fields"] = []
            info["missing_fields"] = sorted(REQUIRED_FORMAT_FIELDS)
        audits.append(info)
    return audits


def new_report(enabled):
    return {
        "enabled": bool(enabled),
        "status": "DISABLED" if not enabled else "PENDING",
        "reason": "AC0 flagging disabled in config" if not enabled else "",
        "required_format_fields": sorted(REQUIRED_FORMAT_FIELDS),
        "thresholds": {
            "min_GQ": GQ_MIN,
            "min_diploid_DP": DIPLOID_DP_MIN,
            "min_haploid_DP": HAPLOID_DP_MIN,
            "min_het_AB": HET_AB_MIN,
            "AB_denominator": "FORMAT/DP",
        },
        "input_files": [],
        "missing_fields": [],
        "missing_values": {},
        "variants_total": None,
        "variants_evaluated": None,
        "variants_skipped_missing_values": None,
        "variants_flagged_ac0": None,
        "multiallelic_sites": None,
        "additional_notes": [
            "AC0 only adds a FILTER tag; no variant rows are removed.",
            "Original VCF FILTER tags are retained, including pre-existing AC0 tags.",
            "AC0 is computed on the remaining samples and GT calls at the start of M4.",
            "Missing FORMAT declarations skip AC0 for the dataset; missing values skip only affected variants.",
            "No-call genotypes are allowed, but a variant with no called GT is not evaluated.",
            "For unsplit multiallelic sites AC0 means zero high-quality copies of ANY ALT allele.",
        ],
    }


def _format_schema_problems(mt, hl):
    """Protect against schemas that differ from the VCF declarations."""
    fields = mt.entry.dtype
    missing = sorted(REQUIRED_FORMAT_FIELDS - set(fields.keys()))
    invalid_types = []
    if "GT" in fields and fields["GT"] != hl.tcall:
        invalid_types.append("GT (expected Hail call)")
    for name in ("GQ", "DP"):
        if name in fields and fields[name] not in (
            hl.tint32, hl.tint64, hl.tfloat32, hl.tfloat64
        ):
            invalid_types.append(f"{name} (expected numeric)")
    if "AD" in fields and fields["AD"] not in (
        hl.tarray(hl.tint32), hl.tarray(hl.tint64),
        hl.tarray(hl.tfloat32), hl.tarray(hl.tfloat64),
    ):
        invalid_types.append("AD (expected array of numeric depths)")
    return missing, invalid_types


def _count_incomplete_genotypes(mt, hl):
    """Count missing data needed to evaluate called genotypes.

    Missing GT is normal for no-calls: those entries contribute no observed ALT
    copies, and are reported but do not invalidate AC0. GQ/DP are required for
    called genotypes; AD is required only for heterozygotes.
    """
    called = hl.is_defined(mt.GT)
    heterozygote = called & mt.GT.is_het()
    ad_valid = hl.or_else(
        hl.is_defined(mt.AD)
        & (hl.len(mt.AD) == hl.len(mt.alleles))
        & hl.all(lambda x: hl.is_defined(x) & (x >= 0), mt.AD),
        False,
    )
    return mt.aggregate_entries(
        hl.struct(
            called_gt=hl.agg.count_where(called),
            missing_gt=hl.agg.count_where(~called),
            missing_gq=hl.agg.count_where(called & hl.is_missing(mt.GQ)),
            missing_dp=hl.agg.count_where(called & hl.is_missing(mt.DP)),
            missing_or_invalid_het_ad=hl.agg.count_where(
                heterozygote & ~ad_valid
            ),
        )
    )


def _flag_when_ac_zero(mt, hl):
    """Add AC0 only to sites with complete, evaluable genotype-quality values.

    FORMAT declarations are checked separately by maybe_flag_ac0(). A missing
    quality value on one *called* genotype makes its row ineligible for AC0,
    not every variant in the merged dataset. No-call GT entries are allowed,
    but a row containing only no-calls is not evaluable. No entries, alleles,
    or variant rows are filtered or modified by this function.

    NOTE: each annotate_rows() creates a new MatrixTable. Expressions must be
    defined using the same current MatrixTable; never carry GT/GQ/AD/DP
    expression handles from before an annotate_rows() into later expressions.
    """
    # Phase 1: make a per-variant eligibility mask. All entry expressions
    # below derive from the same source MT.
    called = hl.is_defined(mt.GT)
    heterozygous = hl.or_else(mt.GT.is_het(), False)
    ad_valid = hl.or_else(
        hl.is_defined(mt.AD)
        & (hl.len(mt.AD) == hl.len(mt.alleles))
        & hl.all(lambda x: hl.is_defined(x) & (x >= 0), mt.AD),
        False,
    )
    entry_complete = (
        ~called
        | (
            hl.is_defined(mt.GQ)
            & hl.is_defined(mt.DP)
            & (~heterozygous | ad_valid)
        )
    )
    mt = mt.annotate_rows(
        ac0_evaluable_tmp=(
            hl.agg.all(entry_complete)
            & hl.agg.any(called)
        )
    )

    # Phase 2: recompute all entry expressions from the updated MT. The
    # earlier `called`, `heterozygous`, etc. must NOT be reused here.
    gt = mt.GT
    dp = mt.DP
    ad = mt.AD
    ad_valid = hl.or_else(
        hl.is_defined(ad)
        & (hl.len(ad) == hl.len(mt.alleles))
        & hl.all(lambda x: hl.is_defined(x) & (x >= 0), ad),
        False,
    )
    het = hl.or_else(gt.is_het(), False)

    # For 0/1 and 1|0 use AD of the non-reference allele. For 1/2
    # both ALT depths must pass the allele-balance threshold.
    ab_pass = hl.if_else(
        het,
        hl.if_else(
            hl.or_else(ad_valid & hl.is_defined(dp) & (dp > 0), False),
            hl.if_else(
                hl.or_else(gt.is_het_ref(), False),
                ad[hl.max(gt[0], gt[1])] / dp > HET_AB_MIN,
                (ad[gt[0]] / dp > HET_AB_MIN)
                & (ad[gt[1]] / dp > HET_AB_MIN),
            ),
            False,
        ),
        True,
    )

    adequate_depth = hl.if_else(
        gt.is_haploid(), dp >= HAPLOID_DP_MIN, dp >= DIPLOID_DP_MIN
    )
    high_quality = hl.or_else(
        hl.is_defined(gt)
        & (mt.GQ >= GQ_MIN)
        & adequate_depth
        & (dp > 0)
        & ab_pass,
        False,
    )
    mt = mt.annotate_rows(
        ac0_adj_stats_tmp=hl.agg.call_stats(
            hl.or_missing(high_quality, gt), mt.alleles
        )
    )

    # Phase 3: each expression below references the MT created in the
    # immediately preceding stage, avoiding Hail cross-source expressions.
    mt = mt.annotate_rows(
        ac0_flag_tmp=(
            mt.ac0_evaluable_tmp
            & (hl.sum(mt.ac0_adj_stats_tmp.AC[1:]) == 0)
        )
    )
    mt = mt.annotate_rows(
        filters=hl.if_else(
            mt.ac0_flag_tmp,
            hl.or_else(mt.filters, hl.empty_set(hl.tstr)).add("AC0"),
            mt.filters,
        )
    )

    # Aggregation forces evaluation before optional QC's nonfatal error
    # handler allows the MatrixTable to advance to normal M4 AF calculation.
    totals = mt.aggregate_rows(
        hl.struct(
            variants=hl.agg.count(),
            evaluated=hl.agg.count_where(mt.ac0_evaluable_tmp),
            skipped=hl.agg.count_where(~mt.ac0_evaluable_tmp),
            flagged=hl.agg.count_where(mt.ac0_flag_tmp),
            multiallelic_sites=hl.agg.count_where(hl.len(mt.alleles) > 2),
            actual_ac0_filter=hl.agg.count_where(mt.filters.contains("AC0")),
        )
    )
    return mt.drop("ac0_adj_stats_tmp", "ac0_flag_tmp", "ac0_evaluable_tmp"), totals


def maybe_flag_ac0(mt, enabled, vcf_dir):
    """Best-effort AC0. On any AC0 problem return the original MT unchanged."""
    report = new_report(enabled)
    if not enabled:
        logging.info("AC0 status=DISABLED (ac0_flagging=false)")
        return mt, report

    try:
        audits = inspect_input_formats(vcf_dir)
        report["input_files"] = audits
        if not audits:
            report["status"] = "SKIPPED_INPUT_AUDIT"
            report["reason"] = "No input VCFs found for AC0 FORMAT inspection"
            logging.warning("AC0 skipped: %s", report["reason"])
            return mt, report

        missing = sorted({
            field for item in audits for field in item["missing_fields"]
        })
        report["missing_fields"] = missing
        problems = [item for item in audits if item["error"]]
        if problems:
            report["status"] = "SKIPPED_INPUT_AUDIT"
            report["reason"] = "Unable to inspect every input VCF header"
            logging.warning("AC0 skipped: %s; files=%s", report["reason"], problems)
            return mt, report
        if missing:
            report["status"] = "SKIPPED_MISSING_FORMAT_FIELDS"
            report["reason"] = "Required FORMAT declarations absent from one or more VCFs"
            logging.warning("AC0 skipped: %s; fields=%s", report["reason"], ", ".join(missing))
            for item in audits:
                if item["missing_fields"]:
                    logging.warning("  VCF %s missing: %s", item["file"], ", ".join(item["missing_fields"]))
            return mt, report

        import hail as hl
        missing, invalid_types = _format_schema_problems(mt, hl)
        report["missing_fields"] = missing
        if missing or invalid_types:
            report["status"] = "SKIPPED_MATRIX_SCHEMA"
            report["reason"] = "Hail entry schema lacks usable required fields"
            report["invalid_types"] = invalid_types
            logging.warning("AC0 skipped: schema missing=%s; invalid=%s", missing, invalid_types)
            return mt, report

        counts = _count_incomplete_genotypes(mt, hl)
        report["missing_values"] = {
            "called_gt": int(counts.called_gt),
            "missing_gt_no_calls": int(counts.missing_gt),
            "missing_gq": int(counts.missing_gq),
            "missing_dp": int(counts.missing_dp),
            "missing_or_invalid_het_ad": int(counts.missing_or_invalid_het_ad),
        }
        # Individual missing FORMAT values are NOT dataset-wide blockers.
        # A site with incomplete called genotypes is excluded from the AC0
        # decision, while other sites are still evaluated and can be flagged.
        annotated_mt, totals = _flag_when_ac_zero(mt, hl)
        report["variants_total"] = int(totals.variants)
        report["variants_evaluated"] = int(totals.evaluated)
        report["variants_skipped_missing_values"] = int(totals.skipped)
        report["variants_flagged_ac0"] = int(totals.flagged)
        report["multiallelic_sites"] = int(totals.multiallelic_sites)

        if totals.evaluated == 0:
            report["status"] = "SKIPPED_MISSING_FORMAT_VALUES"
            report["reason"] = (
                "No variants have sufficient called genotype-quality "
                "information for AC0 evaluation"
            )
            # Preserve the unmodified MT when no sites were evaluated.
            return mt, report
        elif totals.skipped > 0:
            report["status"] = "PARTIALLY_APPLIED_MISSING_FORMAT_VALUES"
            report["reason"] = (
                f"AC0 evaluated for {totals.evaluated} variants; "
                f"{totals.skipped} skipped due to missing genotype-quality values"
            )
        else:
            report["status"] = "APPLIED"
            report["reason"] = (
                "AC0 evaluated for all variants and added to FILTER where "
                "the quality-adjusted ALT allele count is zero"
            )
        logging.info(
            "AC0 status=%s; variants=%s; evaluated=%s; skipped=%s; flagged=%s",
            report["status"], totals.variants, totals.evaluated,
            totals.skipped, totals.flagged,
        )
        return annotated_mt, report

    except Exception as exc:
        # AC0 is best-effort by contract.  Errors related to this optional
        # calculation must not interrupt other pipeline stages.
        report["status"] = "SKIPPED_AC0_ERROR"
        report["reason"] = f"AC0 evaluation failed ({type(exc).__name__}): {exc}"
        logging.exception("AC0 skipped due to optional-QC error; continuing unchanged")
        return mt, report


def write_ac0_reports(report, output_prefix="ac0_status"):
    """Always write a human-readable TSV and detailed machine-readable JSON."""
    Path(f"{output_prefix}.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    with open(f"{output_prefix}.tsv", "w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(("key", "value"))
        for key in (
            "enabled", "status", "reason", "missing_fields", "variants_total",
            "variants_evaluated", "variants_skipped_missing_values",
            "variants_flagged_ac0", "multiallelic_sites",
        ):
            value = report.get(key, "")
            if isinstance(value, list):
                value = ", ".join(str(v) for v in value)
            writer.writerow((key, "" if value is None else value))
        for key, value in report.get("missing_values", {}).items():
            writer.writerow((key, value))
        for item in report.get("input_files", []):
            missing = ", ".join(item["missing_fields"]) or "none"
            writer.writerow((f"VCF: {item['file']}", f"missing: {missing}; error: {item['error'] or 'none'}"))


def append_to_workflow_csv(report, csv_writer):
    """Integrate into existing WorkflowStats_*.csv diagnostics."""
    rows = [
        [],
        ["AC0 flagging status"],
        ["Configured", str(report["enabled"])],
        ["Status", report["status"]],
        ["Reason", report["reason"]],
        ["Required FORMAT", ", ".join(report["required_format_fields"])],
        ["Missing FORMAT", ", ".join(report["missing_fields"]) or "none"],
    ]
    for label, key in (
        ("Variants total", "variants_total"),
        ("Variants evaluated", "variants_evaluated"),
        ("Variants skipped missing values", "variants_skipped_missing_values"),
        ("Variants flagged AC0", "variants_flagged_ac0"),
    ):
        if report.get(key) is not None:
            rows.append([label, report[key]])
    for key, count in report.get("missing_values", {}).items():
        rows.append([key, count])
    for item in report.get("input_files", []):
        if item.get("missing_fields") or item.get("error"):
            rows.append(["Input VCF", item["file"], "Missing: " + ", ".join(item["missing_fields"]), item.get("error") or ""])
    csv_writer(rows)
