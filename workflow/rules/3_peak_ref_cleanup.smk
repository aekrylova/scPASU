configfile: "../config/config.yaml"

rule create_turef:
    input:
        gtf=f"{config['cellrangerrefpath']}genes/genes.gtf",
        work_dir=config["work_dir"]
    params:
        flank=5000,
        final_annotation="gene_symbol",
        goldmine_path=config["goldminepath"]
    output: 
        directory("outputs/hg38_turef/")
    conda:
        "../envs/module_3.yaml"
    log:
        "logs/create_turef.log"
    shell:
        """
        echo 'Assigning TU to peaks...' &> {log}
        Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/create_turef.R -d {input.work_dir}/scripts/3_peak_ref_cleanup/ -g {input.gtf} -o {output}/ -m {params.goldmine_path} -f {params.flank} -a {params.final_annotation} &>> {log}
        """

rule assign_tu_all_filtered_reads:
    input: 
        turef_dir="outputs/hg38_turef/",
        peak_ref_dir="outputs/2d_intersect_run/all_filtered_reads/",
        work_dir=config["work_dir"]
    params:
        chrs=config["chrs"],
        fprefix=f"{config['compartment']}_all_filtered_reads_1sttu_assigned",
        goldmine_path=config["goldminepath"]
    output: 
        directory("outputs/3a_1st_assign_tu/all_filtered_reads/")
    conda:
        "../envs/module_3.yaml"
    log:
        "logs/3a_assign_tu/all_filtered_reads.log"
    shell: 
        """
        echo 'Assigning TU to all filtered reads...' &> {log}
        Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3a_assign_tu.R -d {input.work_dir}/scripts/3_peak_ref_cleanup/ -p {input.peak_ref_dir}/ -m {params.goldmine_path} -f {params.fprefix} -o {output}/ -t {input.turef_dir}/turef.rds -c {params.chrs} &>> {log}

        echo 'Creating BED file...' &>> {log}
        ref={output}/{params.fprefix}_peak_universe_updated.txt
        Rscript {input.work_dir}/scripts/create_bed_file.R -r ${{ref}} -o {output}/ -f {params.fprefix} &>> {log}
        """

rule assign_tu_polyA_reads:
    input: 
        turef_dir="outputs/hg38_turef/",
        peak_ref_dir="outputs/2d_intersect_run/polyA_reads/",
        work_dir=config["work_dir"]
    params:
        chrs=config["chrs"],
        fprefix=f"{config['compartment']}_polyA_reads_1sttu_assigned",
        goldmine_path=config["goldminepath"]
    output: 
        directory("outputs/3a_1st_assign_tu/polyA_reads/")
    conda:
        "../envs/module_3.yaml"
    log:
        "logs/3a_assign_tu/polyA_reads.log"
    shell: 
        """
        echo 'Assigning TU to polyA reads...' &> {log}
        Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3a_assign_tu.R -d {input.work_dir}/scripts/3_peak_ref_cleanup/ -p {input.peak_ref_dir}/ -m {params.goldmine_path} -f {params.fprefix} -o {output}/ -t {input.turef_dir}/turef.rds -c {params.chrs} &>> {log}

        echo 'Creating BED file...' &>> {log}
        ref={output}/{params.fprefix}_peak_universe_updated.txt
        Rscript {input.work_dir}/scripts/create_bed_file.R -r ${{ref}} -o {output}/ -f {params.fprefix} &>> {log}
        """

rule update_peak_ref_all_filtered_reads:
    input: 
        fasta=f"{config['cellrangerrefpath']}fasta/genome.fa",
        bam="outputs/1d_merged_bam/dedup_uniq_genomicAfiltered_merged.bam",
        ref_dir="outputs/3a_1st_assign_tu/all_filtered_reads/",
        work_dir=config["work_dir"]
    params:
        chrs=config["chrs"],
        maxwidth=1000,
        PAS_filtering=config["PAS_filtering"],
        pr_extn="5,50",
        tes_prox_distance=100,
        min_cov=10,
        cores=16,
        fprefix=f"{config['compartment']}_all_filtered_reads",
        goldmine_path=config["goldminepath"]
    output: 
        pr_filter=directory("outputs/3b1_filtered_by_PR_width/all_filtered_reads/"),
        peak_cov=directory("outputs/3b3_peak_coverage_stats/all_filtered_reads/"),
        peak_cov_filt=directory("outputs/3b4_filtered_by_peak_coverage/all_filtered_reads/")
    conda:
        "../envs/module_3.yaml"
    log:
        "logs/3b_update_peak_ref/all_filtered_reads.log"
    shell: 
        """
        set +o pipefail;

        mkdir -p {output.pr_filter}

        echo Update peak reference... &> {log}
        Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3b1_filter_by_pr_width.R -r {input.ref_dir}/{params.fprefix}_1sttu_assigned_peak_universe_updated.txt -o {output.pr_filter}/ -f {params.fprefix} -m {params.maxwidth} &>> {log}

        ref1={input.work_dir}/{output.pr_filter}/{params.fprefix}_pr_filtered_peak_ref.txt
        turef={input.work_dir}/{input.ref_dir}/{params.fprefix}_1sttu_assigned_turef_flankupdated.rds
        
        # Filter by polyA hexamer signals (PAS.fa) and proximity to annotated transcript ends

        if [[ "{params.PAS_filtering}" == "true" ]]; then
            echo Perform filtering by PAS... &>> {log}
            echo Prepare bed input for hexamer checking... &>> {log}
            out2={input.work_dir}/outputs/3b2_filtered_by_PAS_and_TES/all_filtered_reads/
            mkdir -p ${{out2}}
            Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3b2_hexamer_check_bed_input.R -r ${{ref1}} -o ${{out2}} -f {params.fprefix} -e {params.pr_extn} &>> {log}

            bed_input=${{out2}}{params.fprefix}_hexamer_check_input.bed
            fa_output=${{out2}}{params.fprefix}_hexamer_check_input.fa
            echo Extract genomic ranges of peaks from reference genome... &>> {log}
            bedtools getfasta -fi {input.fasta} -bed $bed_input -s -nameOnly -fo $fa_output &>> {log}

            peaks_with_hexamers=${{out2}}{params.fprefix}_peaks_with_hexamer.txt
            peaks_wo_hexamers=${{out2}}{params.fprefix}_peaks_wo_hexamer.txt
            PAS_file={input.work_dir}/scripts/PAS.fa

            echo Check PAS presence in peaks... &>> {log}
            grep -B 1 -f $PAS_file $fa_output | grep '^>' > $peaks_with_hexamers 2>> {log}
            diff -u <(grep '^>' $fa_output| sort) <(sort $peaks_with_hexamers) | grep '^->' > $peaks_wo_hexamers 2>> {log}

            echo Check proximity to TES in peaks without PAS... &>> {log}
            Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3b2_TES_proximity_check.R -b {params.tes_prox_distance} -t $turef -f {params.fprefix} -o $out2 -r $ref1 &>> {log}
            ref2a=${{out2}}{params.fprefix}_peaks_PAS_and_TES_filtered.txt
        else
            echo Skip filtering by PAS... &>> {log}
            out2={input.work_dir}/{output.pr_filter}/
            ref2a=${{out2}}{params.fprefix}_pr_filtered_peak_ref.txt
            fprefix1={params.fprefix}_noPASfiltering_all_filtered_reads
        fi

        # Assign TU - to update peak ref

        echo 2nd TU assignment and peak ref update... &>> {log}
        fprefix2={params.fprefix}_2ndtu_assigned
        Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3a_assign_tu.R -d {input.work_dir}/scripts/3_peak_ref_cleanup/ -m {params.goldmine_path} -f ${{fprefix2}} -o ${{out2}} -t ${{turef}} -b ${{ref2a}} -c {params.chrs} &>> {log}

        echo Create bed file for UCSC uploads... &>> {log}
        ref2b=${{out2}}${{fprefix2}}_peak_universe_updated.txt
        Rscript {input.work_dir}/scripts/create_bed_file.R -r ${{ref2b}} -o ${{out2}} -f ${{fprefix2}} &>> {log}

        # Calculate peak counts

        echo Get peak counts using bam... &>> {log}
        ref_saf=${{out2}}${{fprefix2}}_peak_universe_updated.saf
        mkdir -p {output.peak_cov}
        Rscript {input.work_dir}/scripts/feature_counts.R -b {input.bam} -r ${{ref_saf}} -o {output.peak_cov}/ -f {params.fprefix} -c {params.cores} -i no &>> {log}

        # Calculate per TU peak coverage

        echo Get peak coverage... &>> {log}
        counts_file={input.work_dir}/{output.peak_cov}/{params.fprefix}_peak_count.rds
        Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3b3_peaks_coverage.R -r ${{ref2b}} -c ${{counts_file}} -o {output.peak_cov}/ &>> {log}

        # Plot peak coverage percent

        echo Peak coverage plot... &>> {log}
        pcov_file={input.work_dir}/{output.peak_cov}/{params.fprefix}_peak_count_updated.rds
        w=11
        h=8
        x=100
        y=5000
        Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3b3_plot_peak_cov.R -p ${{pcov_file}} -o {output.peak_cov}/ -m {params.goldmine_path} -f {params.fprefix} -a ${{w}},${{h}},${{x}},${{y}} -s {input.work_dir}/scripts/3_peak_ref_cleanup/ &>> {log}
        
        # Filter peak ref
        echo Filter peak ref by coverage... &>> {log}
        mkdir -p {output.peak_cov_filt}
        Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3b4_filter_peak_ref.R -r ${{ref2b}} -p ${{pcov_file}} -o {output.peak_cov_filt}/ -f {params.fprefix} -m {params.min_cov} &>> {log}

        # Assign TU
        echo 3rd TU assignment... &>> {log}
        ref3={input.work_dir}/{output.peak_cov_filt}/{params.fprefix}_cov_filtered_peak_ref.txt
        fprefix3={params.fprefix}_3rdtu_assigned
        turef=$(ls ${{out2}}*_2ndtu_assigned_turef_flankupdated.rds)
        Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3a_assign_tu.R -d {input.work_dir}/scripts/3_peak_ref_cleanup/ -m {params.goldmine_path} -f ${{fprefix3}} -o {output.peak_cov_filt}/ -t ${{turef}} -b ${{ref3}} -c {params.chrs} &>> {log}

        # Create bed file for UCSC uploads
        echo Create bed file for UCSC uploads... &>> {log}
        ref4={input.work_dir}/{output.peak_cov_filt}/${{fprefix3}}_peak_universe_updated.txt
        Rscript {input.work_dir}/scripts/create_bed_file.R -r ${{ref4}} -o {output.peak_cov_filt}/ -f ${{fprefix3}} &>> {log}
        """

rule update_peak_ref_polyA_reads:
    input: 
        fasta=f"{config['cellrangerrefpath']}fasta/genome.fa",
        bam="outputs/1d_merged_bam/dedup_uniq_genomicAfiltered_merged.bam",
        ref_dir="outputs/3a_1st_assign_tu/polyA_reads/",
        work_dir=config["work_dir"]
    params:
        chrs=config["chrs"],
        maxwidth=1000,
        PAS_filtering=config["PAS_filtering"],
        pr_extn="5,50",
        tes_prox_distance=100,
        min_cov=10,
        cores=16,
        fprefix=f"{config['compartment']}_polyA_reads",
        goldmine_path=config["goldminepath"]
    output: 
        pr_filter=directory("outputs/3b1_filtered_by_PR_width/polyA_reads/"),
        peak_cov=directory("outputs/3b3_peak_coverage_stats/polyA_reads/"),
        peak_cov_filt=directory("outputs/3b4_filtered_by_peak_coverage/polyA_reads/")
    conda:
        "../envs/module_3.yaml"
    log:
        "logs/3b_update_peak_ref/polyA_reads.log"
    shell: 
        """
        set +o pipefail;

        mkdir -p {output.pr_filter}

        echo Update peak reference... &> {log}
        Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3b1_filter_by_pr_width.R -r {input.ref_dir}/{params.fprefix}_1sttu_assigned_peak_universe_updated.txt -o {output.pr_filter}/ -f {params.fprefix} -m {params.maxwidth} &>> {log}

        ref1={input.work_dir}/{output.pr_filter}/{params.fprefix}_pr_filtered_peak_ref.txt
        turef={input.work_dir}/{input.ref_dir}/{params.fprefix}_1sttu_assigned_turef_flankupdated.rds

        # Filter by polyA hexamer signals (PAS.fa) and proximity to annotated transcript ends

        if [[ "{params.PAS_filtering}" == "true" ]]; then
            echo Perform filtering by PAS... &>> {log}
            echo Prepare bed input for hexamer checking... &>> {log}
            out2={input.work_dir}/outputs/3b2_filtered_by_PAS_and_TES/polyA_reads/
            mkdir -p ${{out2}}
            Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3b2_hexamer_check_bed_input.R -r ${{ref1}} -o ${{out2}} -f {params.fprefix} -e {params.pr_extn} &>> {log}

            bed_input=${{out2}}/{params.fprefix}_hexamer_check_input.bed
            fa_output=${{out2}}/{params.fprefix}_hexamer_check_input.fa
            echo Extract genomic ranges of peaks from reference genome... &>> {log}
            bedtools getfasta -fi {input.fasta} -bed $bed_input -s -nameOnly -fo $fa_output &>> {log}

            peaks_with_hexamers=${{out2}}/{params.fprefix}_peaks_with_hexamer.txt
            peaks_wo_hexamers=${{out2}}{params.fprefix}_peaks_wo_hexamer.txt
            PAS_file={input.work_dir}/scripts/PAS.fa

            echo Check PAS presence in peaks... &>> {log}
            grep -B 1 -f $PAS_file $fa_output | grep '^>' > $peaks_with_hexamers 2>> {log}
            diff -u <(grep '^>' $fa_output| sort) <(sort $peaks_with_hexamers) | grep '^->' > $peaks_wo_hexamers 2>> {log}

            echo Check proximity to TES in peaks without PAS... &>> {log}
            Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3b2_TES_proximity_check.R -b {params.tes_prox_distance} -t $turef -f {params.fprefix} -o $out2 -r $ref1 &>> {log}
            ref2a=${{out2}}{params.fprefix}_peaks_PAS_and_TES_filtered.txt
        else
            echo Skip filtering by PAS... &>> {log}
            out2={input.work_dir}/{output.pr_filter}/
            ref2a=${{out2}}{params.fprefix}_pr_filtered_peak_ref.txt
            fprefix1={params.fprefix}_noPASfiltering_polyA_reads
        fi

        # Assign TU - to update peak ref

        echo 2nd TU assignment... &>> {log}
        fprefix2={params.fprefix}_2ndtu_assigned
        Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3a_assign_tu.R -d {input.work_dir}/scripts/3_peak_ref_cleanup/ -m {params.goldmine_path} -f ${{fprefix2}} -o ${{out2}} -t ${{turef}} -b ${{ref2a}} -c {params.chrs} &>> {log}

        echo Create bed file for UCSC uploads... &>> {log}
        ref2b=${{out2}}${{fprefix2}}_peak_universe_updated.txt
        Rscript {input.work_dir}/scripts/create_bed_file.R -r ${{ref2b}} -o ${{out2}} -f ${{fprefix2}} &>> {log}

        # Calculate peak counts

        echo Get peak counts using bam... &>> {log}
        ref_saf=${{out2}}${{fprefix2}}_peak_universe_updated.saf
        mkdir -p {output.peak_cov}
        Rscript {input.work_dir}/scripts/feature_counts.R -b {input.bam} -r ${{ref_saf}} -o {output.peak_cov}/ -f {params.fprefix} -c {params.cores} -i no &>> {log}

        # Calculate per TU peak coverage

        echo Get peak coverage... &>> {log}
        counts_file={input.work_dir}/{output.peak_cov}/{params.fprefix}_peak_count.rds
        Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3b3_peaks_coverage.R -r ${{ref2b}} -c ${{counts_file}} -o {output.peak_cov}/ &>> {log}

        # Plot peak coverage percent

        echo Peak coverage plot... &>> {log}
        pcov_file={input.work_dir}/{output.peak_cov}/{params.fprefix}_peak_count_updated.rds
        w=11
        h=8
        x=100
        y=5000
        Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3b3_plot_peak_cov.R -p ${{pcov_file}} -o {output.peak_cov}/ -m {params.goldmine_path} -f {params.fprefix} -a ${{w}},${{h}},${{x}},${{y}} -s {input.work_dir}/scripts/3_peak_ref_cleanup/ &>> {log}
        
        # Filter peak ref
        echo Filter peak ref by coverage... &>> {log}
        mkdir -p {output.peak_cov_filt}
        Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3b4_filter_peak_ref.R -r ${{ref2b}} -p ${{pcov_file}} -o {output.peak_cov_filt}/ -f {params.fprefix} -m {params.min_cov} &>> {log}

        # Assign TU
        echo 3rd TU assignment... &>> {log}
        ref3={input.work_dir}/{output.peak_cov_filt}/{params.fprefix}_cov_filtered_peak_ref.txt
        fprefix3={params.fprefix}_3rdtu_assigned
        turef=$(ls ${{out2}}*_2ndtu_assigned_turef_flankupdated.rds)
        Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3a_assign_tu.R -d {input.work_dir}/scripts/3_peak_ref_cleanup/ -m {params.goldmine_path} -f ${{fprefix3}} -o {output.peak_cov_filt}/ -t ${{turef}} -b ${{ref3}} -c {params.chrs} &>> {log}

        # Create bed file for UCSC uploads
        echo Create bed file for UCSC uploads... &>> {log}
        ref4={input.work_dir}/{output.peak_cov_filt}/${{fprefix3}}_peak_universe_updated.txt
        Rscript {input.work_dir}/scripts/create_bed_file.R -r ${{ref4}} -o {output.peak_cov_filt}/ -f ${{fprefix3}} &>> {log}
        """

rule merge_two_prongs:
    input: 
        merge_dir_polyA="outputs/3b4_filtered_by_peak_coverage/polyA_reads/",
        merge_dir_all_filtered="outputs/3b4_filtered_by_peak_coverage/all_filtered_reads/",
        turef_dir="outputs/hg38_turef/",
        work_dir=config["work_dir"]
    params:
        chrs=config["chrs"],
        PAS_filtering=config["PAS_filtering"],
        ncores=config["ncores"],
        fprefix=config["compartment"],
        goldmine_path=config["goldminepath"]
    output: 
        directory("outputs/3c_merged_prongs/")
    conda:
        "../envs/module_3.yaml"
    log:
        "logs/3c_merge_two_prongs.log"
    shell: 
        """
        if [ {params.PAS_filtering} == 'true' ]
            then file_prefix={params.fprefix}
        else
            file_prefix={params.fprefix}_noPASfiltering
        fi
   
        echo Merge peaks from two prongs... &> {log}
        Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3c_merge_two_prongs.R -n {params.ncores} -i outputs/3b4_filtered_by_peak_coverage/ -m {params.goldmine_path} -f ${{file_prefix}} -o {output}/ -c {params.chrs} &>> {log}

        echo Assigning TU to peaks... &>> {log}
        ref={output}/${{file_prefix}}_both_prongs_merged_peak_universe_updated.txt
        fprefix=${{file_prefix}}_4thtu_assigned
        Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3a_assign_tu.R -d {input.work_dir}/scripts/3_peak_ref_cleanup/ -m {params.goldmine_path} -f ${{fprefix}} -o {output}/ -t {input.turef_dir}/turef.rds -b ${{ref}} -c {params.chrs} &>> {log}

        echo Create bed file for UCSC uploads.. &>> {log}
        ref={output}/${{fprefix}}_peak_universe_updated.txt
        Rscript {input.work_dir}/scripts/create_bed_file.R -r ${{ref}} -o {output}/ -f ${{fprefix}} &>> {log}
        """

rule peak_classification:
    input: 
        indir="outputs/3c_merged_prongs/",
        work_dir=config["work_dir"]
    params:
        PAS_filtering=config["PAS_filtering"],
        dist_pct=0.1,
        frag_length=500,
        ncores=config["ncores"],
        fprefix=config["compartment"],
        goldmine_path=config["goldminepath"]
    output: 
        directory("outputs/3d_classified_peaks/")
    conda:
        "../envs/module_3.yaml"
    log:
        "logs/3d_peak_classification.log"
    shell: 
        """
        if [ {params.PAS_filtering} == 'true' ]
            then file_prefix={params.fprefix}
        else
            file_prefix={params.fprefix}_noPASfiltering
        fi

        turef={input.indir}/${{file_prefix}}_4thtu_assigned_turef_flankupdated.rds
        ref={input.indir}/${{file_prefix}}_4thtu_assigned_peak_universe_updated.txt

        echo 'Classify peaks and identify split peaks to merge...' &> {log}
        Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3d_peak_classification.R -b $ref -o {output}/ -m {params.goldmine_path} -f ${{file_prefix}} -t $turef -p {params.dist_pct} -l {params.frag_length} -n {params.ncores} -z &>> {log}
        """

rule identify_fragmented_peaks:
    input: 
        fasta=f"{config['cellrangerrefpath']}fasta/genome.fa",
        indir="outputs/3d_classified_peaks/",
        refdir="outputs/3c_merged_prongs/",
        work_dir=config["work_dir"]
    params:
        PAS_filtering=config["PAS_filtering"],
        binsize=1,
        gap_threshold=0,
        min_coverage=10,
        spliced_read_pct_thres=40,
        realign_peak_read_pct_thres=40,
        ncores=config["ncores"],
        fprefix=config["compartment"],
        goldmine_path=config["goldminepath"]
    output: 
        directory("outputs/3e_fragmented_peaks_to_merge/")
    conda:
        "../envs/module_3.yaml"
    log:
        "logs/3e_identify_fragmented_peaks.log"
    shell: 
        r"""
        if [ {params.PAS_filtering} == 'true' ]
            then fprefix={params.fprefix}
        else
            fprefix={params.fprefix}_noPASfiltering
        fi

        peak_saf={input.work_dir}/{input.refdir}/${{fprefix}}_4thtu_assigned_peak_universe_updated.saf
        
        # Functions
        check_empty_files() {{
        local files=("$@")
        local is_empty=1

        for file in "${{files[@]}}"; do
            if [ ! -s "$file" ]; then
                is_empty=0
                break
            fi
        done

        return $is_empty
        }}

        mkdir -p {output} &> {log}

        transcripts={input.work_dir}/{input.indir}/${{fprefix}}_transcripts_causing_peak_fragmentation.txt
        if [ $(tail -n+2 $transcripts | wc -l) -eq 0 ]; then
            echo "There are no transcripts likely causing fragmentation. Exiting script." &>> {log}
            exit 1
        fi
        echo Extract transcripts likely causing peak fragmentation... &>> {log}

        for strand in plus minus; do
            bed_input={input.work_dir}/{output}/${{fprefix}}_transcripts_causing_peak_fragmentation_exons_${{strand}}.bed
            fa_output={input.work_dir}/{output}/${{fprefix}}_transcripts_causing_peak_fragmentation_exons_${{strand}}.fa

            if [ $strand == 'plus' ]; then 
                tail -n+2 $transcripts | awk -v OFS='\t' '{{print $1,$2,$3,$10,0,$5}}' | awk '$6 =="+" {{print}}' > $bed_input 2>> {log}
            elif [ $strand == 'minus' ]; then 
                tail -n+2 $transcripts | awk -v OFS='\t' '{{print $1,$2,$3,$10,0,$5}}' | awk '$6 =="-" {{print}}' > $bed_input 2>> {log}
            fi

            echo Extract exons of these transcripts from reference genome &>> {log}

            bedtools getfasta -fi {input.fasta} -bed $bed_input -s -nameOnly -fo $fa_output &>> {log}

            echo Stitch the exons to obtain these transcripts &>> {log}

            Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3e_stitch_exons.R -i $fa_output -o {input.work_dir}/{output}/${{fprefix}}_transcripts_causing_peak_fragmentation_${{strand}}.fa &>> {log}
        done

        echo Filter bam file for spliced reads overlapping with merge candidates... &>> {log}

        cd {input.work_dir}/{output}/
        
        for strand in plus minus; do
            bamfile={input.work_dir}/outputs/1d_merged_bam/dedup_uniq_genomicAfiltered_merged_${{strand}}.bam
            peak_bed={input.work_dir}/{input.refdir}/${{fprefix}}_4thtu_assigned_${{strand}}.bed
            subset_peak_bed=${{fprefix}}_4thtu_assigned_${{strand}}_subset.bed
            subset_peak_saf=${{fprefix}}_4thtu_assigned_${{strand}}_subset.saf
            subset_bam_file=${{fprefix}}_dedup_uniq_genomicAfiltered_${{strand}}_subset.bam

            if [ $strand == 'plus' ]; then 
                tail -n+2 $transcripts | awk '$5 == "+"{{print}}' | cut -f11 | grep -v '^[[:space:]]*$' | tr ',' '\n' | sort | uniq > ${{subset_peak_bed}}_temp 2>> {input.work_dir}/{log}
            elif [ $strand == 'minus' ]; then 
                tail -n+2 $transcripts | awk '$5 == "-"{{print}}' | cut -f11 | grep -v '^[[:space:]]*$' | tr ',' '\n' | sort | uniq > ${{subset_peak_bed}}_temp 2>> {input.work_dir}/{log}
            fi

            grep -w -f ${{subset_peak_bed}}_temp $peak_bed > $subset_peak_bed 2>> {input.work_dir}/{log}
            echo GeneID$'\t'Chr$'\t'Start$'\t'End$'\t'Strand > $subset_peak_saf 2>> {input.work_dir}/{log}
            grep -w -f ${{subset_peak_bed}}_temp $peak_saf >> $subset_peak_saf 2>> {input.work_dir}/{log}

            echo Count reads for merge candidates... &>> {input.work_dir}/{log}

            samtools view -b -h -L $subset_peak_bed $bamfile > $subset_bam_file 2>> {input.work_dir}/{log}

            Rscript {input.work_dir}/scripts/feature_counts.R -r $subset_peak_saf -b $subset_bam_file -o {input.work_dir}/{output}/ -f ${{fprefix}}_${{strand}}_subset -c {params.ncores} -i no &>> {input.work_dir}/{log}
            echo Count spliced reads for merge candidates... &>> {input.work_dir}/{log}
            spliced_only_sam_temp=${{subset_bam_file%.bam}}_spliced_only.sam

            samtools view -H $subset_bam_file > $spliced_only_sam_temp 2>> {input.work_dir}/{log}
            samtools view $subset_bam_file | awk '($6 ~ /N/)' >> $spliced_only_sam_temp 2>> {input.work_dir}/{log}
            samtools view -bS $spliced_only_sam_temp > ${{spliced_only_sam_temp%.sam}}.bam 2>> {input.work_dir}/{log}
            samtools view $subset_bam_file | awk '($6 ~ /N/)' | awk '{{print "@"$1"\n"$10"\n+\n"$11}}' > ${{subset_bam_file%.bam}}_spliced_only.fastq 2>> {input.work_dir}/{log}

            Rscript {input.work_dir}/scripts/feature_counts.R -r $subset_peak_saf -b ${{spliced_only_sam_temp%.sam}}.bam -o {input.work_dir}/{output}/ -f ${{fprefix}}_${{strand}}_subset_spliced_reads -c {params.ncores} -i no &>> {input.work_dir}/{log}
   
            rm ${{subset_peak_bed}}_temp $spliced_only_sam_temp $subset_bam_file
        done

        echo Call peaks on the realignment of spliced reads from merge candidates against their likely exomes... &>> {input.work_dir}/{log}
        for strand in plus minus; do
            ### Realign spliced reads from merge candidates against transcripts likely causing peak fragmentation
            bowtie2_prefix=${{fprefix}}_transcripts_causing_peak_fragmentation_${{strand}}
            realign_prefix=${{fprefix}}_realign_merge_candidates_${{strand}}
   
            if ! check_empty_files ${{bowtie2_prefix}}.fa; then
                bowtie2-build ${{bowtie2_prefix}}.fa ${{bowtie2_prefix}} &>> {input.work_dir}/{log}
                fastq_file=${{fprefix}}_dedup_uniq_genomicAfiltered_${{strand}}_subset_spliced_only.fastq
                bowtie2 --local --threads {params.ncores} -x ${{bowtie2_prefix}} -X 2000 -U ${{fastq_file}} -S ${{realign_prefix}}.sam &>> {input.work_dir}/{log}
   
                ### Preparing files for IGV 

                samtools faidx ${{bowtie2_prefix}}.fa &>> {input.work_dir}/{log}
                samtools view -bS ${{realign_prefix}}.sam > ${{realign_prefix}}.bam 2>> {input.work_dir}/{log}
                samtools sort ${{realign_prefix}}.bam -@ {params.ncores} -o ${{realign_prefix}}.sorted.bam > ${{realign_prefix}}.sorted.bam 2>> {input.work_dir}/{log}
                samtools index ${{realign_prefix}}.sorted.bam &>> {input.work_dir}/{log}
                rm ${{realign_prefix}}.bam ${{realign_prefix}}.sam

                ### Call peaks
                bamCoverage -b ${{realign_prefix}}.sorted.bam -o ${{realign_prefix}}.bw -p {params.ncores} -v --binSize {params.binsize} &>> {input.work_dir}/{log}
                bigWigToBedGraph ${{realign_prefix}}.bw ${{realign_prefix}}.bedGraph &>> {input.work_dir}/{log}
            else
                echo chr$'\t'start$'\t'end$'\t'score > ${{realign_prefix}}.bedGraph 2>> {input.work_dir}/{log}
            fi 
        done

        Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3e_peak_calling_from_bedGraph.R -f $fprefix -b {input.work_dir}/{output}/ -o {input.work_dir}/{output}/ -m {params.goldmine_path} -g {params.gap_threshold} -p {params.min_coverage} &>> {input.work_dir}/{log}

        realign_peak=${{fprefix}}_spliced_reads_realign_peaks.txt
        for strand in plus minus; do   
            ### Count reads for realign peaks from spliced reads
            saf=${{fprefix}}_spliced_reads_realign_peaks_${{strand}}.saf
            realign_prefix=${{fprefix}}_realign_merge_candidates_${{strand}}
            bamfile=${{realign_prefix}}.sorted.bam
            echo GeneID$'\t'Chr$'\t'Start$'\t'End$'\t'Strand > $saf 2>> {input.work_dir}/{log}
            if [ $strand == 'plus' ]; then
                awk -v strand="+" -v OFS='\t' '$6 == strand {{print $4,$1,$2,$3,$6}}' $realign_peak >> $saf 2>> {input.work_dir}/{log}
            elif [ $strand == 'minus' ]; then 
                awk -v strand="-" -v OFS='\t' '$6 == strand {{print $4,$1,$2,$3,$6}}' $realign_peak >> $saf 2>> {input.work_dir}/{log}
            fi
            Rscript {input.work_dir}/scripts/feature_counts.R -r $saf -b $bamfile -o {input.work_dir}/{output}/ -f ${{fprefix}}_spliced_reads_realign_peaks_${{strand}} -c {params.ncores} -i no &>> {input.work_dir}/{log}
        done

        echo Identify fragmented peaks to merge... &>> {input.work_dir}/{log}
        ref={input.work_dir}/{input.indir}/${{fprefix}}_peak_universe_classified.txt 
        Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/3e_identify_fragmented_peaks.R -t $transcripts -p $realign_peak -d {input.work_dir}/{output}/ -m {params.goldmine_path} -f $fprefix -o {input.work_dir}/{output}/ -s {params.spliced_read_pct_thres} -r {params.realign_peak_read_pct_thres} -b $ref -n {params.ncores} &>> {input.work_dir}/{log}
        """ 

rule make_filtered_tracks:
    input: 
        plus="outputs/1d_merged_bam/dedup_uniq_genomicAfiltered_merged_plus.bam",
        minus="outputs/1d_merged_bam/dedup_uniq_genomicAfiltered_merged_minus.bam",
        ref_dir="outputs/3d_classified_peaks/",
    output: 
        directory("outputs/filtered_tracks/")
    params:
        fprefix=config["compartment"],
        PAS_filtering=config["PAS_filtering"],
    conda:
        "../envs/module_3.yaml"
    log:
        "logs/make_filtered_tracks.log"
    shell: 
        """
        mkdir -p {output} 

        if [ {params.PAS_filtering} == 'true' ]
            then file_prefix={params.fprefix}
        else
            file_prefix={params.fprefix}_noPASfiltering
        fi

        txt_ref={input.ref_dir}/${{file_prefix}}_peak_universe_classified.txt
        bed={input.ref_dir}/${{file_prefix}}_peak_universe_classified.bed
        ref={input.ref_dir}/${{file_prefix}}_peak_universe_classified_sorted.bed

        cut -f1,2,3,21 $txt_ref > $bed 2> {log}
        sed -i '1d' $bed &>> {log}
        sort -k1,1 -k2,2n $bed > $ref 2>> {log}

        # Minus
        
        bam_output={output}/${{file_prefix}}_peakref_minus.bam
        samtools view -b -h -L $ref {input.minus} > $bam_output  2>> {log}
        samtools index $bam_output &>> {log}
        bamCoverage -b $bam_output -o ${{bam_output%.bam}}_polyAsubset.bw -p 36 -v --binSize 1 &>> {log}

        # Plus

        bam_output={output}/${{file_prefix}}_peakref_plus.bam
        samtools view -b -h -L $ref {input.plus} > $bam_output 2>> {log}
        samtools index $bam_output &>> {log}
        bamCoverage -b $bam_output -o ${{bam_output%.bam}}_polyAsubset.bw -p 36 -v --binSize 1 &>> {log}
        """

rule plot_nuc_freq:
    input: 
        ref_dir="outputs/3c_merged_prongs/",
        work_dir=config["work_dir"]
    output: 
        directory(f"outputs/nuc_freq_plot_{config['kmer_size']}-mer/")
    params:
        genome=config['genomename'], 
        fprefix=config["compartment"],
        PAS_filtering=config["PAS_filtering"],
        kmer_to_plot="first",
        kmer_size=config["kmer_size"],
        ncores=config["ncores"],
        goldmine_path=config["goldminepath"]
    conda:
        "../envs/module_3.yaml"
    log:
        "logs/plot_nuc_freq.log"
    shell: 
        """
        if [ {params.PAS_filtering} == 'true' ]
            then file_prefix={params.fprefix}
        else
            file_prefix={params.fprefix}_noPASfiltering
        fi

        ref_file={input.work_dir}/{input.ref_dir}/${{file_prefix}}_4thtu_assigned_peak_universe_updated.txt
        echo Nuc frequency plotting... $> {log}
        Rscript {input.work_dir}/scripts/3_peak_ref_cleanup/nuc_freq_plot.R -c {params.ncores} -f ${{file_prefix}} -o {output}/ -m {params.goldmine_path} -r $ref_file -k {params.kmer_size} -b -n {params.kmer_to_plot} -g {params.genome} &>> {log}
        """

rule pa_ref_stats:
    input:
        final_folder="outputs/3e_fragmented_peaks_to_merge/",
        work_dir=config["work_dir"]
    output: 
        f"outputs/{config['compartment']}_PA_ref_stats.txt"
    params:
        fprefix=config["compartment"],
        PAS_filtering=config["PAS_filtering"],
        chrs=config["chrs"]
    conda:
        "../envs/module_3.yaml"
    log:
        "logs/pa_ref_stats.log"
    shell:
        """
        set +u

        if [ {params.PAS_filtering} == 'true' ]
            then PAS_status=PASfiltering
        else
            PAS_status=noPASfiltering
        fi

        chrs=$(echo {params.chrs} | sed 's/,/\\\\|/g')
            
        refstats={input.work_dir}/{output} &> {log}

        if test -f $refstats; then echo Delete the existing $(basename "$refstats") &>> {log} ; rm $refstats; fi

        echo 'Strand'$'\t''Plus'$'\t''Minus'$'\t''Both' > $refstats 2>> {log}

        for strand in plus minus; do
            export bam_file=$(ls {input.work_dir}/outputs/1d_merged_bam/*${{strand}}.bam | grep -v 'genomicAfiltered')
            export bam_stats=${{bam_file%.bam}}.stats.txt
            samtools stats $bam_file > $bam_stats 2>> {log}
            if [ $strand == 'plus' ]; then
                export num_plus_dedup_uniq_reads=$(grep 'raw total sequences:' ${{bam_stats}} | cut -f3)
            elif [ $strand == 'minus' ]; then
                export num_minus_dedup_uniq_reads=$(grep 'raw total sequences:' ${{bam_stats}} | cut -f3)
            fi
        done

        echo 'Number of deduplicated, uniquely mapped reads:'$'\t'${{num_plus_dedup_uniq_reads}}$'\t'${{num_minus_dedup_uniq_reads}}$'\t'$((num_plus_dedup_uniq_reads+num_minus_dedup_uniq_reads)) >> $refstats 2>> {log}

        for prong in all_filtered_reads polyA_reads; do
            echo $prong >> $refstats 2>> {log}
            for strand in plus minus; do
                if [ $prong == 'all_filtered_reads' ]; then
                    export bam_dir=outputs/1d_merged_bam/
                    export bam_file=$(ls ${{bam_dir}}/dedup_uniq_genomicAfiltered_merged_${{strand}}.bam)
                    export polyA_dir={input.work_dir}/outputs/2b_polya/after_genomicAfiltering/
                elif [ $prong == 'polyA_reads' ]; then
                    export bam_dir={input.work_dir}/outputs/2b_polya/before_genomicAfiltering/bam_files/
                    export bam_file=$(ls ${{bam_dir}}*${{strand}}*nowrongstrand.bam)
                    export polyA_dir={input.work_dir}/outputs/2b_polya/before_genomicAfiltering/
                fi

                bam_stats=${{bam_file%.bam}}.stats.txt
                samtools stats $bam_file > $bam_stats 2>> {log}
            
                if [ $strand == 'plus' ]; then
                    export num_plus_reads=$(grep 'raw total sequences:' ${{bam_stats}} | cut -f3)
                    export num_polyA_plus_reads=$(cat ${{polyA_dir}}*${{strand}}_nowrongstrand.bed | wc -l)
                elif [ $strand == 'minus' ]; then
                    export num_minus_reads=$(grep 'raw total sequences:' ${{bam_stats}} | cut -f3)
                    export num_polyA_minus_reads=$(cat ${{polyA_dir}}*${{strand}}_nowrongstrand.bed | wc -l)
                fi
        
            done
        
            echo 'Number of reads used for peak calling:'$'\t'${{num_plus_reads}}$'\t'${{num_minus_reads}}$'\t'$((num_plus_reads+num_minus_reads)) >> $refstats 2>> {log}
            echo 'Number of polyA junction reads:'$'\t'${{num_polyA_plus_reads}}$'\t'${{num_polyA_minus_reads}}$'\t'$((num_polyA_plus_reads+num_polyA_minus_reads)) >> $refstats 2>> {log}
        
            for strand in plus minus; do
                new_peaks=() 
                peaks=$(grep -w $chrs {input.work_dir}/outputs/2a_peaks/${{prong}}/*${{strand}}*peaks_summits.txt | cut -f10)
                for peak in "${{peaks[@]}}"; do
                    new_peak=$(echo "$peak" | awk '{{ if (substr($0,length($0),1) ~ /[[:alpha:]]/) {{ sub(/[[:alpha:]]$/, "", $0) }} print }}') &>> {log}
                    new_peaks+=("$new_peak")
                done
            
                if [ $strand == 'plus' ]; then
                    export num_plus_peaks=$(echo $new_peaks | tr ' ' '\n' | sort | uniq | wc -l)
                elif [ $strand == 'minus' ]; then
                    export num_minus_peaks=$(echo $new_peaks | tr ' ' '\n' | sort | uniq | wc -l)
                fi
            done      
        
            echo 'Number of MACS2 peaks:'$'\t'${{num_plus_peaks}}$'\t'${{num_minus_peaks}}$'\t'$((num_plus_peaks+num_minus_peaks)) >> $refstats
        
            num_plus_peaks_aftersplitting=$(cat {input.work_dir}/outputs/2c2_split_peaks/${{prong}}/*plus_peaks.narrowPeak | wc -l)
            num_minus_peaks_aftersplitting=$(cat {input.work_dir}/outputs/2c2_split_peaks/${{prong}}/*minus_peaks.narrowPeak | wc -l) 
            echo 'Number of MACS2 peaks after splitting:'$'\t'${{num_plus_peaks_aftersplitting}}$'\t'${{num_minus_peaks_aftersplitting}}$'\t'$((num_plus_peaks_aftersplitting+num_minus_peaks_aftersplitting)) >> $refstats
            
            num_plus_polyA_peaks=$(cut -f5 {input.work_dir}/outputs/3a_1st_assign_tu/${{prong}}/*_peak_universe.txt | tail -n+2 | grep -w '+' | wc -l)
            num_minus_polyA_peaks=$(cut -f5 {input.work_dir}/outputs/3a_1st_assign_tu/${{prong}}/*_peak_universe.txt | tail -n+2 | grep -w '-' | wc -l) 
            echo 'Number of polyA-supported peaks:'$'\t'${{num_plus_polyA_peaks}}$'\t'${{num_minus_polyA_peaks}}$'\t'$((num_plus_polyA_peaks+num_minus_polyA_peaks)) >> $refstats
            
            num_plus_TUassigned_peaks=$(cat {input.work_dir}/outputs/3a_1st_assign_tu/${{prong}}/*_plus.bed | wc -l)
            num_minus_TUassigned_peaks=$(cat {input.work_dir}/outputs/3a_1st_assign_tu/${{prong}}/*_minus.bed | wc -l)
            echo 'Number of TU-assigned peaks:'$'\t'${{num_plus_TUassigned_peaks}}$'\t'${{num_minus_TUassigned_peaks}}$'\t'$((num_plus_TUassigned_peaks+num_minus_TUassigned_peaks)) >> $refstats
            
            num_plus_PRfiltered_peaks=$(cat {input.work_dir}/outputs/3b1_filtered_by_PR_width/${{prong}}/*_pr_filtered_peak_ref.txt | tail -n+2 | grep -w '+' | wc -l)
            num_minus_PRfiltered_peaks=$(cat {input.work_dir}/outputs/3b1_filtered_by_PR_width/${{prong}}/*_pr_filtered_peak_ref.txt | tail -n+2 | grep -w '-' | wc -l)
            echo 'Number of PR-filtered peaks:'$'\t'${{num_plus_PRfiltered_peaks}}$'\t'${{num_minus_PRfiltered_peaks}}$'\t'$((num_plus_PRfiltered_peaks+num_minus_PRfiltered_peaks)) >> $refstats
            
            if [ {params.PAS_filtering} == 'true' ]; then    
                num_plus_PAS_peaks=$(tail -n+2 {input.work_dir}/outputs/3b2_filtered_by_PAS_and_TES/${{prong}}/*peaks_with_hexamer.txt | cut -f4 | grep -w '+' | wc -l)
                num_minus_PAS_peaks=$(tail -n+2 {input.work_dir}/outputs/3b2_filtered_by_PAS_and_TES/${{prong}}/*peaks_with_hexamer.txt | cut -f4 | grep -w '-' | wc -l)
                num_plus_noPAS_nearTES_peaks=$(tail -n+2 {input.work_dir}/outputs/3b2_filtered_by_PAS_and_TES/${{prong}}/*peaks_wo_hexamer_near_TES.txt | cut -f4 | grep -w '+' | wc -l)
                num_minus_noPAS_nearTES_peaks=$(tail -n+2 {input.work_dir}/outputs/3b2_filtered_by_PAS_and_TES/${{prong}}/*peaks_wo_hexamer_near_TES.txt | cut -f4 | grep -w '-' | wc -l)
                num_plus_PASfiltered_peaks=$(cat {input.work_dir}/outputs/3b2_filtered_by_PAS_and_TES/${{prong}}/*_plus.bed | wc -l)
                num_minus_PASfiltered_peaks=$(cat {input.work_dir}/outputs/3b2_filtered_by_PAS_and_TES/${{prong}}/*_minus.bed | wc -l)
                echo 'Number of peaks with PAS:'$'\t'${{num_plus_PAS_peaks}}$'\t'${{num_minus_PAS_peaks}}$'\t'$((num_plus_PAS_peaks+num_minus_PAS_peaks)) >> $refstats
                echo 'Number of peaks without PAS but near a TES:'$'\t'${{num_plus_noPAS_nearTES_peaks}}$'\t'${{num_minus_noPAS_nearTES_peaks}}$'\t'$((num_plus_noPAS_nearTES_peaks+num_minus_noPAS_nearTES_peaks)) >> $refstats
                echo 'Number of PAS-filtered peaks:'$'\t'${{num_plus_PASfiltered_peaks}}$'\t'${{num_minus_PASfiltered_peaks}}$'\t'$((num_plus_PASfiltered_peaks+num_minus_PASfiltered_peaks)) >> $refstats
                fprefix={params.fprefix}
            else
                fprefix={params.fprefix}_noPASfiltering
            fi

            num_plus_CovFiltered_peaks=$(cat {input.work_dir}/outputs/3b4_filtered_by_peak_coverage/${{prong}}/${{fprefix}}_${{prong}}_3rdtu_assigned_plus.bed | wc -l)
            num_minus_CovFiltered_peaks=$(cat {input.work_dir}/outputs/3b4_filtered_by_peak_coverage/${{prong}}/${{fprefix}}_${{prong}}_3rdtu_assigned_minus.bed | wc -l)
            echo 'Number of coverage-filtered peaks:'$'\t'${{num_plus_CovFiltered_peaks}}$'\t'${{num_minus_CovFiltered_peaks}}$'\t'$((num_plus_CovFiltered_peaks+num_minus_CovFiltered_peaks)) >> $refstats
        done

        echo Final peaks >> $refstats
        num_plus_final_peaks=$(cat {input.work_dir}/outputs/3c_merged_prongs/${{fprefix}}_4thtu_assigned_plus.bed | wc -l)
        num_minus_final_peaks=$(cat {input.work_dir}/outputs/3c_merged_prongs/${{fprefix}}_4thtu_assigned_minus.bed | wc -l)
        echo 'Number of final peaks:'$'\t'${{num_plus_final_peaks}}$'\t'${{num_minus_final_peaks}}$'\t'$((num_plus_final_peaks+num_minus_final_peaks)) >> $refstats

        echo Peak Classification: >> $refstats
        tail -n+2 {input.work_dir}/outputs/3d_classified_peaks/${{fprefix}}_peak_universe_classified.txt | cut -f31 | sort | uniq -c >> $refstats

        num_plus_fragmented_peaks=$(tail -n+2 {input.work_dir}/{input.final_folder}/${{fprefix}}_final_peak_universe_updated.txt | awk -v col=4 '{{ if ($col == "+") {{print $0}} }}' | cut -f32 | grep -v 'NA' | wc -l)
        num_minus_fragmented_peaks=$(tail -n+2 {input.work_dir}/{input.final_folder}/${{fprefix}}_final_peak_universe_updated.txt | awk -v col=4 '{{ if ($col == "-") {{print $0}} }}' | cut -f32 | grep -v 'NA' | wc -l)
        echo 'Number of fragmented peaks:'$'\t'${{num_plus_fragmented_peaks}}$'\t'${{num_minus_fragmented_peaks}}$'\t'$((num_plus_fragmented_peaks+num_minus_fragmented_peaks)) >> $refstats

        """