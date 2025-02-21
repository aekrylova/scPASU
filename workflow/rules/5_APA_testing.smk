configfile: "../config/config.yaml"

rule merge_per_sample_peakcounts:
    input: 
        counts_dir="outputs/4b_cellranger_peakcount/",
        work_dir=config["work_dir"]
    output: 
        directory("outputs/5a_merged_cellranger_peakcount/")
    params:
        seurat_obj=config["seurat_obj_file"],
        file_prefix=config["compartment"],
        samples=config["samples"],
    conda:
        "../envs/module_5.yaml"
    log:
        "logs/5a_merge_peakcounts.log"
    script: 
        "../scripts/5_APA_testing/5a_merge_per_sample_counts.R"

rule merge_per_sample_genecounts:
    input: 
        counts_dir="outputs/4b_cellranger_genecount/",
        work_dir=config["work_dir"]
    output: 
        directory("outputs/5a_merged_cellranger_genecount/")
    params:
        seurat_obj=config["seurat_obj_file"],
        file_prefix=config["compartment"],
        samples=config["samples"],
    conda:
        "../envs/module_5.yaml"
    log:
        "logs/5a_merge_genecounts.log"
    script: 
        "../scripts/5_APA_testing/5a_merge_per_sample_counts.R"

rule APA_testing:
    input: 
        counts_dir="outputs/5a_merged_cellranger_peakcount/",
        meta_file=config["metadata_file"],
        peak_ref_dir="outputs/3e_fragmented_peaks_to_merge/",
        work_dir=config["work_dir"]
    output: 
        directory("outputs/5b_APA_testing/")
    params:
        group_one=config["group_one"],
        group_two=config["group_two"],
        meta_testing_column=config["meta_testing_column"],
        file_prefix=config["compartment"],
        min_cell_expr_pct=10
    conda:
        "../envs/module_5.yaml"
    log:
        "logs/5b_APA_testing.log"
    script: 
        "../scripts/5_APA_testing/5b1_APA_testing.R"

rule DEG_testing:
    input: 
        counts_dir="outputs/5a_merged_cellranger_genecount/",
        meta_file=config["metadata_file"],
        work_dir=config["work_dir"]
    output: 
        directory("outputs/5b_DEG_testing/")
    params:
        group_one=config["group_one"],
        group_two=config["group_two"],
        meta_testing_column=config["meta_testing_column"],
        file_prefix=config["compartment"],
    conda:
        "../envs/module_5.yaml"
    log:
        "logs/5b_DEG_testing.log"
    script: 
        "../scripts/5_APA_testing/5b2_DEG_testing.R"

rule make_fractional_tracks:
    input: 
        bed_dir="outputs/5b_APA_testing/",
        work_dir=config["work_dir"]
    output: 
        directory("outputs/fractional_usage_tracks/")
    params:
        group_one=config["group_one"],
        group_two=config["group_two"],
        fasta=f"{config['cellrangerrefpath']}fasta/genome.fa",
        refpath=config["cellrangerrefpath"]
    conda:
        "../envs/module_5.yaml"
    log:
        "logs/make_fractional_tracks.log"
    shell: 
        """
        mkdir -p {output} &> {log}
        comparison={params.group_one}_v_{params.group_two}
        genome_fai={params.fasta}.fai
        chrome_sizes={params.refpath}/fasta/GRCh38_v2024_Chrom.sizes
        cut -f1,2 $genome_fai > $chrome_sizes 2>> {log}

        # Group one
        strand=plus
        bed4={input.bed_dir}/${{comparison}}-{params.group_one}_${{strand}}.bedGraph
        bigwig={output}/${{comparison}}-{params.group_one}_${{strand}}.bw

        sort -k1,1 -k2,2n $bed4 > ${{bed4%.bedGraph}}_sorted.bedGraph 2>> {log}
        bedGraphToBigWig ${{bed4%.bedGraph}}_sorted.bedGraph $chrome_sizes ${{bigwig}} &>> {log}

        strand=minus
        bed4={input.bed_dir}/${{comparison}}-{params.group_one}_${{strand}}.bedGraph
        bigwig={output}/${{comparison}}-{params.group_one}_${{strand}}.bw

        sort -k1,1 -k2,2n $bed4 > ${{bed4%.bedGraph}}_sorted.bedGraph 2>> {log}
        bedGraphToBigWig ${{bed4%.bedGraph}}_sorted.bedGraph $chrome_sizes ${{bigwig}} &>> {log}

        # Group two
        strand=plus
        bed4={input.bed_dir}/${{comparison}}-{params.group_two}_${{strand}}.bedGraph
        bigwig={output}/${{comparison}}-{params.group_two}_${{strand}}.bw

        sort -k1,1 -k2,2n $bed4 > ${{bed4%.bedGraph}}_sorted.bedGraph 2>> {log}
        bedGraphToBigWig ${{bed4%.bedGraph}}_sorted.bedGraph $chrome_sizes ${{bigwig}} &>> {log}

        strand=minus
        bed4={input.bed_dir}/${{comparison}}-{params.group_two}_${{strand}}.bedGraph
        bigwig={output}/${{comparison}}-{params.group_two}_${{strand}}.bw

        sort -k1,1 -k2,2n $bed4 > ${{bed4%.bedGraph}}_sorted.bedGraph 2>> {log}
        bedGraphToBigWig ${{bed4%.bedGraph}}_sorted.bedGraph $chrome_sizes ${{bigwig}} &>> {log}
        """
rule make_track_annotations:
    input: 
        bed_dir="outputs/3c_merged_prongs/",
        work_dir=config["work_dir"]
    output: 
        directory("outputs/track_annotations/")
    params:
        fprefix=config["compartment"],
        PAS_filtering=config["PAS_filtering"],
        genomerefpath=config["cellrangerrefpath"]
    conda:
        "../envs/module_5.yaml"
    log:
        "logs/make_track_annotations.log"
    shell: 
        """
        mkdir -p {output} &> {log}

        if [ {params.PAS_filtering} == 'true' ]
            then file_prefix={params.fprefix}
        else
            file_prefix={params.fprefix}_noPASfiltering
        fi

        genome_fai={params.genomerefpath}/fasta/genome.fa.fai
        chrome_sizes={params.genomerefpath}/fasta/GRCh38_v2024_Chrom.sizes
        cut -f1,2 $genome_fai > $chrome_sizes 2>> {log}

        bed={input.bed_dir}/${{file_prefix}}_4thtu_assigned_plus.bed
        bigBed={output}/${{file_prefix}}_plus.bigBed
        sort -k1,1 -k2,2n $bed > ${{bed%.bed}}_sorted.bed 2>> {log}
        bedToBigBed -tab ${{bed%.bed}}_sorted.bed $chrome_sizes ${{bigBed}} &>> {log}

        bed={input.bed_dir}/${{file_prefix}}_4thtu_assigned_minus.bed
        bigBed={output}/${{file_prefix}}_minus.bigBed
        sort -k1,1 -k2,2n $bed > ${{bed%.bed}}_sorted.bed 2>> {log}
        bedToBigBed -tab ${{bed%.bed}}_sorted.bed $chrome_sizes ${{bigBed}} &>> {log}
        """

