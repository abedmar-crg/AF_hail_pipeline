nextflow.enable.dsl=2


/*
 * ============================================================
 * M1.1 - CONVERT AND MERGE VCFs
 * ============================================================
 */

process M1_CONVERT_AND_MERGE_VCFS {

    tag "M1 - convert and merge VCFs"

    input:
    path input_root
    path config_file, name: 'config.yaml'

    output:
    path 'm1_convert.mt', emit: mt

    path '*.log',
        optional: true,
        emit: logs

    path '*.csv',
        optional: true,
        emit: csvs

    path '*.png',
        optional: true,
        emit: pngs
    
    path '*.html',
        optional: true,
        emit: htmls    

    script:
    """
    echo "M1.1 - Convert and merge VCFs"

    test -d "${input_root}"

    python /app/nextflow_stage.py \
        convert \
        --output m1_convert.mt
    """
}


/*
 * ============================================================
 * M1.2 - SPLIT MULTIALLELIC
 * ============================================================
 */

process M1_SPLIT_MULTIALLELIC {

    tag "M1 - split multiallelic"

    input:
    path input_mt
    path config_file, name: 'config.yaml'

    output:
    path 'm1_split.mt', emit: mt

    path '*.log',
        optional: true,
        emit: logs

    path '*.csv',
        optional: true,
        emit: csvs

    path '*.png',
        optional: true,
        emit: pngs
    
    path '*.html',
        optional: true,
        emit: htmls    

    script:
    """
    echo "M1.2 - Split multiallelic variants"

    python /app/nextflow_stage.py \
        split \
        --input "${input_mt}" \
        --output m1_split.mt
    """
}


/*
 * ============================================================
 * M1.3 - GENOTYPE FILTERING
 * ============================================================
 */

process M1_GENOTYPE_FILTERING {

    tag "M1 - genotype filtering"

    input:
    path input_mt
    path config_file, name: 'config.yaml'

    output:
    path 'm1_genotype.mt', emit: mt

    path '*.log',
        optional: true,
        emit: logs

    path '*.csv',
        optional: true,
        emit: csvs

    path '*.png',
        optional: true,
        emit: pngs
    
    path '*.html',
        optional: true,
        emit: htmls    

    script:
    """
    echo "M1.3 - Genotype filtering"

    python /app/nextflow_stage.py \
        genotype \
        --input "${input_mt}" \
        --output m1_genotype.mt
    """
}


/*
 * ============================================================
 * M1.4 - VARIANT FILTERING
 * ============================================================
 */

process M1_VARIANT_FILTERING {

    tag "M1 - variant filtering"

    input:
    path input_mt
    path config_file, name: 'config.yaml'

    output:
    path 'm1_variant.mt', emit: mt

    path '*.log',
        optional: true,
        emit: logs

    path '*.csv',
        optional: true,
        emit: csvs

    path '*.png',
        optional: true,
        emit: pngs
    
    path '*.html',
        optional: true,
        emit: htmls    

    script:
    """
    echo "M1.4 - Variant filtering"

    python /app/nextflow_stage.py \
        variant \
        --input "${input_mt}" \
        --output m1_variant.mt
    """
}


/*
 * ============================================================
 * M1.5 - SAMPLE FILTERING: HARD THRESHOLDS
 * ============================================================
 */

process M1_SAMPLE_FILTERING_HARD {

    tag "M1 - sample hard thresholds"

    input:
    path input_mt
    path config_file, name: 'config.yaml'

    output:
    path 'm1_sample_hard.mt', emit: mt

    path '*.log',
        optional: true,
        emit: logs

    path '*.csv',
        optional: true,
        emit: csvs

    path '*.png',
        optional: true,
        emit: pngs
    
    path '*.html',
        optional: true,
        emit: htmls    

    script:
    """
    echo "M1.5 - Sample filtering: hard thresholds"

    python /app/nextflow_stage.py \
        sample_hard \
        --input "${input_mt}" \
        --output m1_sample_hard.mt
    """
}


/*
 * ============================================================
 * M1.6 - SAMPLE FILTERING: MAD THRESHOLDS
 * ============================================================
 */

process M1_SAMPLE_FILTERING_MAD {

    tag "M1 - sample MAD thresholds"

    input:
    path input_mt
    path config_file, name: 'config.yaml'

    output:
    path 'm1_sample_mad.mt', emit: mt

    path '*.log',
        optional: true,
        emit: logs

    path '*.csv',
        optional: true,
        emit: csvs

    path '*.png',
        optional: true,
        emit: pngs
    
    path '*.html',
        optional: true,
        emit: htmls    

    script:
    """
    echo "M1.6 - Sample filtering: MAD thresholds"

    python /app/nextflow_stage.py \
        sample_mad \
        --input "${input_mt}" \
        --output m1_sample_mad.mt
    """
}


/*
 * ============================================================
 * M2 - RELATEDNESS
 * ============================================================
 */

process M2_RELATEDNESS {

    tag "M2 - relatedness"

    input:
    path input_mt
    path config_file, name: 'config.yaml'

    output:
    path 'm2_unrelated.mt', emit: mt

    path '*.log',
        optional: true,
        emit: logs

    path '*.csv',
        optional: true,
        emit: csvs

    path '*.png',
        optional: true,
        emit: pngs
    
    path '*.html',
        optional: true,
        emit: htmls    

    script:
    """
    echo "M2 - Relatedness"

    python /app/nextflow_stage.py \
        m2 \
        --input "${input_mt}" \
        --output m2_unrelated.mt
    """
}


/*
 * ============================================================
 * M3 - ANCESTRY
 * ============================================================
 */

process M3_ANCESTRY {

    tag "M3 - ancestry"

    input:
    path input_mt
    path config_file, name: 'config.yaml'

    output:
    path 'm3_ancestry.mt', emit: mt

    /*
     * These files exist only when ancestry inference
     * actually produces them.
     */
    path 'ancestrySNPs.vcf',
        optional: true,
        emit: ancestry_vcf

    path '*-GrafAnc_results',
        optional: true,
        emit: ancestry_results

    path '*.log',
        optional: true,
        emit: logs

    path '*.csv',
        optional: true,
        emit: csvs

    path '*.png',
        optional: true,
        emit: pngs
    
    path '*.html',
        optional: true,
        emit: htmls    

    script:
    """
    echo "M3 - Ancestry"

    python /app/nextflow_stage.py \
        m3 \
        --input "${input_mt}" \
        --output m3_ancestry.mt
    """
}


/*
 * ============================================================
 * M4 - AF ANNOTATION
 * ============================================================
 */

process M4_AF_ANNOTATION {

    tag "M4 - AF annotation"

    input:
    path input_mt
    path config_file, name: 'config.yaml'

    output:
    path "${params.final_vcf_name}", emit: vcf

    path "${params.final_vcf_name}.tbi", emit: index

    path '*.log',
        optional: true,
        emit: logs

    path '*.csv',
        optional: true,
        emit: csvs

    // Dedicated AC0 status reports are always retained, including when
    // AC0 was disabled or skipped because FORMAT information was missing.
    path 'ac0_status.tsv', emit: ac0_status
    path 'ac0_status.json', emit: ac0_details

    path '*.png',
        optional: true,
        emit: pngs
    
    path '*.html',
        optional: true,
        emit: htmls    

    script:
    """
    echo "M4 - AF annotation"

    python /app/nextflow_stage.py \
        m4 \
        --input "${input_mt}" \
        --output "${params.final_vcf_name}"
    """
}


/*
 * ============================================================
 * WORKFLOW
 * ============================================================
 */

workflow {

    main:

    /*
     * --------------------------------------------------------
     * INPUTS
     * --------------------------------------------------------
     */

    input_root_ch = channel.value(
        file(params.input_dir)
    )

    config_ch = channel.value(
        file(params.hail_config)
    )


    /*
     * --------------------------------------------------------
     * M1
     * --------------------------------------------------------
     *
     * These are sequential because every stage operates on
     * the MatrixTable produced by the preceding stage.
     */

    M1_CONVERT_AND_MERGE_VCFS(
        input_root_ch,
        config_ch
    )

    M1_SPLIT_MULTIALLELIC(
        M1_CONVERT_AND_MERGE_VCFS.out.mt,
        config_ch
    )

    M1_GENOTYPE_FILTERING(
        M1_SPLIT_MULTIALLELIC.out.mt,
        config_ch
    )

    M1_VARIANT_FILTERING(
        M1_GENOTYPE_FILTERING.out.mt,
        config_ch
    )

    M1_SAMPLE_FILTERING_HARD(
        M1_VARIANT_FILTERING.out.mt,
        config_ch
    )

    M1_SAMPLE_FILTERING_MAD(
        M1_SAMPLE_FILTERING_HARD.out.mt,
        config_ch
    )


    /*
     * --------------------------------------------------------
     * M2
     * --------------------------------------------------------
     */

    M2_RELATEDNESS(
        M1_SAMPLE_FILTERING_MAD.out.mt,
        config_ch
    )


    /*
     * --------------------------------------------------------
     * M3
     * --------------------------------------------------------
     */

    M3_ANCESTRY(
        M2_RELATEDNESS.out.mt,
        config_ch
    )


    /*
     * --------------------------------------------------------
     * M4
     * --------------------------------------------------------
     */

    M4_AF_ANNOTATION(
        M3_ANCESTRY.out.mt,
        config_ch
    )

    /*
     * ========================================================
     * GRAPH CHANNELS
     * ========================================================
     *
     * Combine PNG and HTML plots from each process.
     */
    
    m1_convert_plots_ch = M1_CONVERT_AND_MERGE_VCFS.out.pngs.mix(M1_CONVERT_AND_MERGE_VCFS.out.htmls)
    
    m1_split_plots_ch = M1_SPLIT_MULTIALLELIC.out.pngs.mix(M1_SPLIT_MULTIALLELIC.out.htmls)
    
    m1_genotype_plots_ch = M1_GENOTYPE_FILTERING.out.pngs.mix(M1_GENOTYPE_FILTERING.out.htmls)
    
    m1_variant_plots_ch = M1_VARIANT_FILTERING.out.pngs.mix(M1_VARIANT_FILTERING.out.htmls)
    
    m1_hard_plots_ch = M1_SAMPLE_FILTERING_HARD.out.pngs.mix(M1_SAMPLE_FILTERING_HARD.out.htmls)
    
    m1_mad_plots_ch = M1_SAMPLE_FILTERING_MAD.out.pngs.mix(M1_SAMPLE_FILTERING_MAD.out.htmls)
    
    m2_plots_ch = M2_RELATEDNESS.out.pngs.mix(M2_RELATEDNESS.out.htmls)
    
    m3_plots_ch = M3_ANCESTRY.out.pngs.mix(M3_ANCESTRY.out.htmls)
    
    m4_plots_ch = M4_AF_ANNOTATION.out.pngs.mix(M4_AF_ANNOTATION.out.htmls)
    
    /*
     * ========================================================
     * PUBLISH
     * ========================================================
     *
     * The processes ALWAYS create the MatrixTables needed by
     * downstream stages.
     *
     * These parameters only decide whether a copy should also
     * be retained under results/.
     *
     * channel.empty() means:
     *    keep it in work/
     *    but do not publish it to results/
     *
     * Keep every assignment on ONE line.
     */

    publish:

    original_matrix = params.save_original_matrix ? M1_CONVERT_AND_MERGE_VCFS.out.mt : channel.empty()

    qced_matrix = params.save_qced_matrix ? M1_SAMPLE_FILTERING_MAD.out.mt : channel.empty()

    relatedness_matrix = params.save_relatedness_matrix ? M2_RELATEDNESS.out.mt : channel.empty()


    ancestry_vcf = params.save_ancestry_results ? M3_ANCESTRY.out.ancestry_vcf : channel.empty()

    ancestry_results = params.save_ancestry_results ? M3_ANCESTRY.out.ancestry_results : channel.empty()


    m1_convert_logs = params.save_logs ? M1_CONVERT_AND_MERGE_VCFS.out.logs : channel.empty()

    m1_convert_csvs = params.save_logs ? M1_CONVERT_AND_MERGE_VCFS.out.csvs : channel.empty()


    m1_split_logs = params.save_logs ? M1_SPLIT_MULTIALLELIC.out.logs : channel.empty()

    m1_split_csvs = params.save_logs ? M1_SPLIT_MULTIALLELIC.out.csvs : channel.empty()


    m1_genotype_logs = params.save_logs ? M1_GENOTYPE_FILTERING.out.logs : channel.empty()

    m1_genotype_csvs = params.save_logs ? M1_GENOTYPE_FILTERING.out.csvs : channel.empty()


    m1_variant_logs = params.save_logs ? M1_VARIANT_FILTERING.out.logs : channel.empty()

    m1_variant_csvs = params.save_logs ? M1_VARIANT_FILTERING.out.csvs : channel.empty()


    m1_hard_logs = params.save_logs ? M1_SAMPLE_FILTERING_HARD.out.logs : channel.empty()

    m1_hard_csvs = params.save_logs ? M1_SAMPLE_FILTERING_HARD.out.csvs : channel.empty()


    m1_mad_logs = params.save_logs ? M1_SAMPLE_FILTERING_MAD.out.logs : channel.empty()

    m1_mad_csvs = params.save_logs ? M1_SAMPLE_FILTERING_MAD.out.csvs : channel.empty()


    m2_logs = params.save_logs ? M2_RELATEDNESS.out.logs : channel.empty()

    m2_csvs = params.save_logs ? M2_RELATEDNESS.out.csvs : channel.empty()


    m3_logs = params.save_logs ? M3_ANCESTRY.out.logs : channel.empty()

    m3_csvs = params.save_logs ? M3_ANCESTRY.out.csvs : channel.empty()


    m4_logs = params.save_logs ? M4_AF_ANNOTATION.out.logs : channel.empty()

    m4_csvs = params.save_logs ? M4_AF_ANNOTATION.out.csvs : channel.empty()


    // AC0 status reports are published even if save_logs=false.
    ac0_status = M4_AF_ANNOTATION.out.ac0_status
    ac0_details = M4_AF_ANNOTATION.out.ac0_details

    final_vcf = params.save_final_vcf ? M4_AF_ANNOTATION.out.vcf : channel.empty()

    final_vcf_index = params.save_final_vcf ? M4_AF_ANNOTATION.out.index : channel.empty()

    m1_convert_plots = params.save_plots ? m1_convert_plots_ch : channel.empty()
    m1_split_plots = params.save_plots ? m1_split_plots_ch : channel.empty()
    m1_genotype_plots = params.save_plots ? m1_genotype_plots_ch : channel.empty()
    m1_variant_plots = params.save_plots ? m1_variant_plots_ch : channel.empty()
    m1_hard_plots = params.save_plots ? m1_hard_plots_ch : channel.empty()
    m1_mad_plots = params.save_plots ? m1_mad_plots_ch : channel.empty()
    m2_plots = params.save_plots ? m2_plots_ch : channel.empty()
    m3_plots = params.save_plots ? m3_plots_ch : channel.empty()
    m4_plots = params.save_plots ? m4_plots_ch : channel.empty()
}


/*
 * ============================================================
 * RESULTS DIRECTORY STRUCTURE
 * ============================================================
 */

output {

    /*
     * --------------------------------------------------------
     * MATRICES
     * --------------------------------------------------------
     */

    original_matrix {
        path 'matrices/original'
        mode 'copy'
    }

    qced_matrix {
        path 'matrices/qced'
        mode 'copy'
    }

    relatedness_matrix {
        path 'matrices/relatedness'
        mode 'copy'
    }


    /*
     * --------------------------------------------------------
     * ANCESTRY
     * --------------------------------------------------------
     */

    ancestry_vcf {
        path 'ancestry'
        mode 'copy'
    }

    ancestry_results {
        path 'ancestry'
        mode 'copy'
    }


    /*
     * --------------------------------------------------------
     * M1 CONVERT / MERGE LOGS
     * --------------------------------------------------------
     */

    m1_convert_logs {
        path 'logs/M1_CONVERT'
        mode 'copy'
    }

    m1_convert_csvs {
        path 'logs/M1_CONVERT'
        mode 'copy'
    }


    /*
     * --------------------------------------------------------
     * M1 SPLIT LOGS
     * --------------------------------------------------------
     */

    m1_split_logs {
        path 'logs/M1_SPLIT'
        mode 'copy'
    }

    m1_split_csvs {
        path 'logs/M1_SPLIT'
        mode 'copy'
    }


    /*
     * --------------------------------------------------------
     * M1 GENOTYPE LOGS
     * --------------------------------------------------------
     */

    m1_genotype_logs {
        path 'logs/M1_GENOTYPE'
        mode 'copy'
    }

    m1_genotype_csvs {
        path 'logs/M1_GENOTYPE'
        mode 'copy'
    }


    /*
     * --------------------------------------------------------
     * M1 VARIANT LOGS
     * --------------------------------------------------------
     */

    m1_variant_logs {
        path 'logs/M1_VARIANT'
        mode 'copy'
    }

    m1_variant_csvs {
        path 'logs/M1_VARIANT'
        mode 'copy'
    }


    /*
     * --------------------------------------------------------
     * M1 SAMPLE HARD LOGS
     * --------------------------------------------------------
     */

    m1_hard_logs {
        path 'logs/M1_SAMPLE_HARD'
        mode 'copy'
    }

    m1_hard_csvs {
        path 'logs/M1_SAMPLE_HARD'
        mode 'copy'
    }


    /*
     * --------------------------------------------------------
     * M1 SAMPLE MAD LOGS
     * --------------------------------------------------------
     */

    m1_mad_logs {
        path 'logs/M1_SAMPLE_MAD'
        mode 'copy'
    }

    m1_mad_csvs {
        path 'logs/M1_SAMPLE_MAD'
        mode 'copy'
    }


    /*
     * --------------------------------------------------------
     * M2 LOGS
     * --------------------------------------------------------
     */

    m2_logs {
        path 'logs/M2_RELATEDNESS'
        mode 'copy'
    }

    m2_csvs {
        path 'logs/M2_RELATEDNESS'
        mode 'copy'
    }


    /*
     * --------------------------------------------------------
     * M3 LOGS
     * --------------------------------------------------------
     */

    m3_logs {
        path 'logs/M3_ANCESTRY'
        mode 'copy'
    }

    m3_csvs {
        path 'logs/M3_ANCESTRY'
        mode 'copy'
    }


    /*
     * --------------------------------------------------------
     * M4 LOGS
     * --------------------------------------------------------
     */

    m4_logs {
        path 'logs/M4_AF'
        mode 'copy'
    }

    m4_csvs {
        path 'logs/M4_AF'
        mode 'copy'
    }


    /*
     * --------------------------------------------------------
     * FINAL RECALCULATED VCF
     * --------------------------------------------------------
     */

    final_vcf {
        path '.'
        mode 'copy'
    }

    final_vcf_index {
        path '.'
        mode 'copy'
    }

        /*
     * ========================================================
     * plots
     * ========================================================
     */
    
    m1_convert_plots {
        path 'plots/M1_CONVERT'
        mode 'copy'
    }
    
    m1_split_plots {
        path 'plots/M1_SPLIT'
        mode 'copy'
    }
    
    m1_genotype_plots {
        path 'plots/M1_GENOTYPE'
        mode 'copy'
    }
    
    m1_variant_plots {
        path 'plots/M1_VARIANT'
        mode 'copy'
    }
    
    m1_hard_plots {
        path 'plots/M1_SAMPLE_HARD'
        mode 'copy'
    }
    
    m1_mad_plots {
        path 'plots/M1_SAMPLE_MAD'
        mode 'copy'
    }
    
    m2_plots {
        path 'plots/M2_RELATEDNESS'
        mode 'copy'
    }
    
    m3_plots {
        path 'plots/M3_ANCESTRY'
        mode 'copy'
    }
    
    m4_plots {
        path 'plots/M4_AF'
        mode 'copy'
    }
    /* AC0 audit is independent of optional log retention. */
    ac0_status {
        path 'qc/ac0'
        mode 'copy'
    }

    ac0_details {
        path 'qc/ac0'
        mode 'copy'
    }

}
