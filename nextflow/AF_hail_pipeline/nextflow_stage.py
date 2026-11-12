#!/usr/bin/env python3

import argparse
import logging
import os
from pathlib import Path

import hail as hl

import utils
import M1_preprocessing as m1
import M2_unrelated_samples as m2
import M3_ancestry as m3
import M4_af_annotation as m4


# All Nextflow tasks receive config.yaml in their working directory.
config = utils.load_config("config.yaml")


def init_hail():
    """
    Initialise Hail using the same logic as the original
    vcf-af-pipeline.py.
    """

    if config.get("cluster", False):
        hl.init(
            spark_conf={
                "spark.driver.memory":
                    config["spark_driver_memory"],

                "spark.executor.memory":
                    config["spark_executor_memory"],

                "spark.executor.instances":
                    config["spark_executor_instances"],

                "spark.executor.cores":
                    config["spark_executor_cores"],

                "spark.rpc.askTimeout":
                    config["spark_rpc_askTimeout"],

                "spark.sql.shuffle.partitions":
                    config["spark_sql_shuffle_partitions"],

                "spark.local.dir":
                    config["spark_local_dir"],

                "spark.network.timeout":
                    config["spark_network_timeout"],
            },
            tmp_dir=config["tmp_dir"],
            local_tmpdir=config["local_tmpdir"],
        )

    else:
        hl.init()


def initialise_stage(stage_name):
    """
    Create logging/statistics files for each independent process.
    """

    utils.load_logging()
    utils.csv_creator()

    logging.info("=" * 60)
    logging.info(f"NEXTFLOW STAGE: {stage_name}")
    logging.info("=" * 60)


def read_mt(path):
    logging.info(f"Reading MatrixTable: {path}")
    return hl.read_matrix_table(path)


def write_mt(mt, path):
    logging.info(f"Writing MatrixTable: {path}")
    mt.write(path, overwrite=True)


# ============================================================
# M1
# ============================================================

def run_convert(output_path):

    initialise_stage("M1_CONVERT_AND_MERGE_VCFS")

    # Import VCF(s), but let this wrapper perform the final write.
    mt = m1.convert_and_merge_vcfs(
        config["vcf_dir"],
        output_path,
        config["ref_gen"],
        write_output=False,
    )

    # This happened between conversion and split_multiallelic
    # in the original pipeline.
    mt = utils.rename_chr(
        mt,
        config["ref_gen"]
    )

    write_mt(mt, output_path)


def run_split(input_path, output_path):

    initialise_stage("M1_SPLIT_MULTIALLELIC")

    mt = read_mt(input_path)

    if config.get("split_multiallelic", False):

        original_size = mt.count()

        mt = m1.split_multiallelic(
            mt,
            original_size
        )

    else:
        logging.warning(
            "split_multiallelic=false; passing MatrixTable unchanged"
        )

    write_mt(mt, output_path)


def run_genotype_filtering(input_path, output_path):

    initialise_stage("M1_GENOTYPE_FILTERING")

    mt = read_mt(input_path)

    if config.get("genotype_filtering", False):

        mt = m1.genotype_filtering(mt)

    else:
        logging.warning(
            "genotype_filtering=false; passing MatrixTable unchanged"
        )

    write_mt(mt, output_path)


def run_variant_filtering(input_path, output_path):

    initialise_stage("M1_VARIANT_FILTERING")

    mt = read_mt(input_path)

    if config.get("variant_filtering", False):

        mt = m1.variant_filtering(mt)

    else:
        logging.warning(
            "variant_filtering=false; passing MatrixTable unchanged"
        )

    write_mt(mt, output_path)


def run_sample_hard(input_path, output_path):

    initialise_stage("M1_SAMPLE_FILTERING_HARD_THRESHOLDS")

    mt = read_mt(input_path)

    if config.get("sample_filtering", False):

        basename = os.path.basename(
            config["mt_from_vcf"].rstrip("/")
        )

        mt = m1.sample_filtering_hard_thresholds(
            mt,
            config["seq_type"],
            basename,
        )

    else:
        logging.warning(
            "sample_filtering=false; hard filtering skipped"
        )

    write_mt(mt, output_path)


def run_sample_mad(input_path, output_path):

    initialise_stage("M1_SAMPLE_FILTERING_MAD_THRESHOLDS")

    mt = read_mt(input_path)

    if config.get("sample_filtering", False):

        basename = os.path.basename(
            config["mt_from_vcf"].rstrip("/")
        )

        mt = m1.sample_filtering_mad_thresholds(
            mt,
            config["seq_type"],
            basename,
        )

    else:
        logging.warning(
            "sample_filtering=false; MAD filtering skipped"
        )

    # This reproduces the final QC sanity check in
    # the original pipeline.
    n_rows, n_cols = mt.count()

    logging.info(
        f"Final M1 dataset: {n_rows} variants, {n_cols} samples"
    )

    if n_rows == 0:
        raise ValueError(
            "All variants were removed during M1 preprocessing"
        )

    if n_cols == 0:
        raise ValueError(
            "All samples were removed during M1 preprocessing"
        )

    write_mt(mt, output_path)


# ============================================================
# M2
# ============================================================

def run_m2(input_path, output_path):

    initialise_stage("M2_RELATEDNESS")

    mt = read_mt(input_path)

    if config.get("delete_related", False):

        _, n_samples = mt.count()

        if n_samples == 1:

            logging.info(
                "Only one sample present; relatedness filtering skipped"
            )

        else:

            mt = m2.delete_related_samples(mt)

    else:

        logging.warning(
            "delete_related=false; M2 passes MatrixTable unchanged"
        )

    write_mt(mt, output_path)


# ============================================================
# M3
# ============================================================

def ensure_m3_resources():
    """
    M3 currently accesses GrafAnc resources through relative paths.

    Each Nextflow task has its own working directory, therefore
    create task-local links to resources contained in the image.
    """

    grafanc_snps = Path("GrafAnc_SNPs")

    if not grafanc_snps.exists():
        grafanc_snps.symlink_to(
            "/app/GrafAnc_SNPs",
            target_is_directory=True
        )

    grafanc_binary = Path("grafanc")

    if not grafanc_binary.exists():
        grafanc_binary.symlink_to(
            "/usr/local/bin/grafanc"
        )


def run_m3(input_path, output_path):

    initialise_stage("M3_ANCESTRY")

    mt = read_mt(input_path)

    infer = config.get("infer_ancestry", False)
    submitted = config.get("submit_ancestry", False)

    if infer and submitted:
        raise ValueError(
            "infer_ancestry and submit_ancestry cannot both be true"
        )

    if infer:

        ensure_m3_resources()

        # subset_matrix() uses the configured MT path to determine
        # where its temporary ancestry VCF should be written.
        #
        # Point it at this Nextflow task instead of the old shared
        # /data/work MatrixTable.
        current_mt = os.path.abspath(input_path)

        m3.config["mt_afterQC"] = current_mt
        m3.config["mt_from_vcf"] = current_mt
        m3.config["preprocessing"] = True

        ancestry_vcf, ancestry_count = m3.subset_matrix(mt)

        logging.info(
            f"Ancestry SNPs found: {ancestry_count}"
        )

        if ancestry_count > 0:

            ancestry_results = m3.call_grafanc(
                ancestry_vcf,
                current_mt
            )

            if ancestry_results is not None:

                mt = m3.annotate_ancestry(
                    ancestry_results,
                    mt
                )

        else:

            logging.warning(
                "No ancestry SNPs found; ancestry annotation skipped"
            )

    elif submitted:

        logging.info(
            f"Using submitted ancestry: "
            f"{config['ancestry_information']}"
        )

        mt = m3.annotate_ancestry(
            config["ancestry_information"],
            mt
        )

    else:

        logging.warning(
            "No ancestry processing requested; "
            "M3 passes MatrixTable unchanged"
        )

    write_mt(mt, output_path)


# ============================================================
# M4
# ============================================================

def run_m4(input_path, output_path):

    initialise_stage("M4_AF_ANNOTATION")

    if not config.get("af_annotation", False):

        raise ValueError(
            "af_annotation=false. "
            "M4_AF_ANNOTATION requires af_annotation=true."
        )

    mt = read_mt(input_path)

    # Optional AC0 flagging is independent of M1 genotype filtering and
    # NEVER removes variants. Missing input FORMAT information is non-fatal.
    mt, ac0_report = m4.maybe_flag_ac0(
        mt, config.get("ac0_flagging", False), config["vcf_dir"]
    )
    m4.write_ac0_reports(ac0_report)
    m4.append_to_workflow_csv(ac0_report, utils.csv_writer)

    # --------------------------------------------------------
    # Sex
    # --------------------------------------------------------

    if config.get("infer_sex", False):

        mt = utils.impute_sex(mt)

        mt, results_sex_agg = m4.stats_by_sex(mt)

    else:

        logging.warning(
            "Sex inference disabled; sex-stratified AFs skipped"
        )

        results_sex_agg = {}


    # --------------------------------------------------------
    # Ancestry
    # --------------------------------------------------------

    if (
        config.get("infer_ancestry", False)
        or config.get("submit_ancestry", False)
    ):

        try:

            mt, results_ancestry_agg = m4.stats_by_ancestry(mt)

        except Exception as exc:

            logging.warning(
                f"Ancestry-stratified AFs skipped: {exc}"
            )

            results_ancestry_agg = ""

    else:

        results_ancestry_agg = ""


    # --------------------------------------------------------
    # AF calculation
    # --------------------------------------------------------

    (
        mt,
        AF_total,
        AC_total,
        AN_total,
        nhom_total,
        nhet_total,
        nhemi_total,
        AF_sex,
        AF_ancestry,
        AF_ancestry_sex,
    ) = m4.af_by_sex_ancestry(
        mt,
        results_sex_agg,
        results_ancestry_agg,
    )

    mt = m4.annotate_new_vcf(
        mt,
        AF_total,
        AC_total,
        AN_total,
        nhom_total,
        nhet_total,
        nhemi_total,
        AF_sex,
        AF_ancestry,
        AF_ancestry_sex,
    )

    # Override the old config output path.
    # The final VCF should now belong to the Nextflow task.
    m4.config["final_vcf_AF"] = os.path.abspath(
        output_path
    )

    m4.export_new_vcf(mt, ac0_applied=(ac0_report["status"] in (
        "APPLIED", "PARTIALLY_APPLIED_MISSING_FORMAT_VALUES"
    )))


# ============================================================
# CLI
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "stage",
        choices=[
            "convert",
            "split",
            "genotype",
            "variant",
            "sample_hard",
            "sample_mad",
            "m2",
            "m3",
            "m4",
        ],
    )

    parser.add_argument(
        "--input",
        required=False,
    )

    parser.add_argument(
        "--output",
        required=True,
    )

    args = parser.parse_args()

    if args.stage != "convert" and not args.input:

        parser.error(
            f"--input is required for stage {args.stage}"
        )

    init_hail()

    try:

        if args.stage == "convert":

            run_convert(
                args.output
            )

        elif args.stage == "split":

            run_split(
                args.input,
                args.output
            )

        elif args.stage == "genotype":

            run_genotype_filtering(
                args.input,
                args.output
            )

        elif args.stage == "variant":

            run_variant_filtering(
                args.input,
                args.output
            )

        elif args.stage == "sample_hard":

            run_sample_hard(
                args.input,
                args.output
            )

        elif args.stage == "sample_mad":

            run_sample_mad(
                args.input,
                args.output
            )

        elif args.stage == "m2":

            run_m2(
                args.input,
                args.output
            )

        elif args.stage == "m3":

            run_m3(
                args.input,
                args.output
            )

        elif args.stage == "m4":

            run_m4(
                args.input,
                args.output
            )

    finally:

        hl.stop()


if __name__ == "__main__":
    main()
