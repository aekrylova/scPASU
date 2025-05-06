configfile: "../config/config.yaml"

rule create_gtf_db:
    input:
        work_dir=config["work_dir"],
        gtf=f"{config['cellrangerrefpath']}genes/genes.gtf", 
        polyA_dir=config["polyAfilterpath"]
    output: 
        db="outputs/1b_filtered_bam/gtfdb.db"
    conda:
        "../envs/module_1.yaml"
    log:
        "logs/create_gtf_db.log"
    shell: 
        """
        python {input.polyA_dir}/polyAfilter.py createDB -v {input.gtf} {input.work_dir}/{output.db} &>> {log}
        """

rule preprocess_bam_dedup:
    input:
        f"{config['cellrangeroutputpath']}/{{sample}}/outs/possorted_genome_bam.bam",
    output:
        "outputs/1a1_dedup_bam/{sample}_dedup.bam"
    conda:
        "../envs/module_1.yaml"
    log:
        "logs/1a1_preprocess_bam_dedup/{sample}.log"
    shell: 
        """
        echo 'Deduplicate reads...' &> {log}
        umi_tools dedup -I {input} -S {output} --method=unique --extract-umi-method=tag --umi-tag=UB --cell-tag=CB &>> {log}
        echo 'Indexing the last bam file...' &>> {log}
        samtools index {output} &>> {log}
        """

rule preprocess_bam_cleanup:
    input:
        "outputs/1a1_dedup_bam/{sample}_dedup.bam"
    output:
        "outputs/1a2_clean_bam/{sample}_dedup_uniq.bam"
    conda:
        "../envs/module_1.yaml" 
    log:
        "logs/1a2_preprocess_bam_cleanup/{sample}.log"
    shell: 
        """
        echo 'Clean up bam...' &> {log}
        samtools view -b -q 1 {input} > {output} 2>> {log}
        echo 'Indexing the last bam file...' &>> {log}
        samtools index {output} &>> {log}
        """

rule genomicA_filtering_create_trans_file:
    input: 
        bams="outputs/1a2_clean_bam/{sample}_dedup_uniq.bam",
        db="outputs/1b_filtered_bam/gtfdb.db",
        polyA_dir=config["polyAfilterpath"]
    output: 
        intron="outputs/1b_filtered_bam/{sample}_trans_intron.pkl",
        nointron="outputs/1b_filtered_bam/{sample}_trans_nointron.pkl"
    conda:
        "../envs/module_1.yaml"
    log:
        "logs/1b_genomicA_filtering_trans_file/{sample}.log"
    shell: 
        """
        echo 'Create the TRANS file including intronic coverage...' &> {log}
        python {input.polyA_dir}/polyAfilter.py createTRANS -i {input.db} {input.bams} {output.intron} &>> {log}
        
        echo 'Create the TRANS file without including intronic coverage...' &>> {log}
        python {input.polyA_dir}/polyAfilter.py createTRANS {input.db} {input.bams} {output.nointron} &>> {log}

        """

rule genomicA_filtering_first_round:
    input: 
        bams="outputs/1a2_clean_bam/{sample}_dedup_uniq.bam",
        fasta=f"{config['cellrangerrefpath']}fasta/genome.fa",
        intron="outputs/1b_filtered_bam/{sample}_trans_intron.pkl",
        polyA_dir=config["polyAfilterpath"]
    output: 
        "outputs/1b_filtered_bam/{sample}_dedup_uniq_filtered_1stround.bam"
    params:
        MISM=config["mism"],
        NTHREADS=config["nthreads"],
        COVLEN=config["cov_len"],
        MINSNRLEN=config["min_snr_len"]
    conda:
        "../envs/module_1.yaml"
    log:
        "logs/1b_genomicA_filtering_first_round/{sample}.log"
    shell: 
        """
        echo 'Filter out the genomicA reads including intronic coverage...' &> {log}
        python {input.polyA_dir}/polyAfilter.py BAMfilter -m {params.MISM} -i -v -p {params.NTHREADS} {params.COVLEN} {params.MINSNRLEN} {input.bams} {input.fasta} {input.intron} -o {output} &>> {log}
        
        echo 'Indexing the last bam file...' &>> {log}
        samtools index {output} &>> {log}
        """

rule genomicA_filtering_second_round:
    input: 
        firstrd="outputs/1b_filtered_bam/{sample}_dedup_uniq_filtered_1stround.bam",
        fasta=f"{config['cellrangerrefpath']}fasta/genome.fa",
        nointron="outputs/1b_filtered_bam/{sample}_trans_nointron.pkl",
        polyA_dir=config["polyAfilterpath"]
    output: 
        "outputs/1b_filtered_bam/{sample}_dedup_uniq_filtered_2ndround.bam"
    conda:
        "../envs/module_1.yaml"
    params:
        MISM=config["mism"],
        NTHREADS=config["nthreads"],
        COVLEN=config["cov_len"],
        MINSNRLEN=config["min_snr_len"]
    log:
        "logs/1b_genomicA_filtering_second_round/{sample}.log"
    shell: 
        """
        echo 'Filter out the genomicA reads without including intronic coverage...' &> {log}
        python {input.polyA_dir}/polyAfilter.py BAMfilter -m {params.MISM} -v -p {params.NTHREADS} {params.COVLEN} {params.MINSNRLEN} {input.firstrd} {input.fasta} {input.nointron} -o {output} &>> {log}
        
        echo 'Indexing the last bam file...' &>> {log}
        samtools index {output} &>> {log}
        """

rule subset_bam:
    input:
        before_filt="outputs/1a2_clean_bam/{sample}_dedup_uniq.bam",
        after_filt="outputs/1b_filtered_bam/{sample}_dedup_uniq_filtered_2ndround.bam",
        work_dir=config["work_dir"]
    output:
        subset_before_filt="outputs/1c_subset_bam/{sample}_dedup_uniq_fibroblasts.bam",
        subset_after_filt="outputs/1c_subset_bam/{sample}_dedup_uniq_filtered_fibroblasts.bam"
    params:
        ncores=config["ncores"],
        barcodes=f"{config['barcodepath']}/{{sample}}_subset_barcodes.tsv",
        subset=config["subset"]
    conda:
        "../envs/module_1.yaml"
    log:
        "logs/1c_subset_bam/{sample}.log"
    shell:
        """
        if [ {params.subset} == 'true' ]; then
            echo 'Subset-bam for the fibroblast compartment before genomicA filtering...' &> {log}
            subset-bam -b {input.before_filt} -c {params.barcodes} -o {output.subset_before_filt} --log-level debug --cores {params.ncores} &>> {log}
            echo 'Indexing the last bam file...' &>> {log}
            samtools index {output.subset_before_filt} &>> {log}

            echo 'Subset-bam for the fibroblast compartment after genomicA filtering...' &>> {log}
            subset-bam -b {input.after_filt} -c {params.barcodes} -o {output.subset_after_filt} ---log-level debug --cores {params.ncores} &>> {log}
            echo 'Indexing the last bam file...' &>> {log}
            samtools index {output.subset_after_filt} &>> {log}
        else
            echo 'No subsetting performed, bams from folders 1a2 and 1b will be copied directly to 1c' &> {log}
            cp {input.work_dir}/{input.before_filt} {input.work_dir}/{output.subset_before_filt} &>> {log}
            cp {input.work_dir}/{input.after_filt} {input.work_dir}/{output.subset_after_filt} &>> {log}
            
            echo 'Indexing the copied files...'
            samtools index {input.work_dir}/{output.subset_before_filt} &>> {log}
            samtools index {input.work_dir}/{output.subset_after_filt} &>> {log}
        fi

        """

rule merge_bam_before_filt:
    input: 
        expand("outputs/1c_subset_bam/{sample}_dedup_uniq_fibroblasts.bam", sample=config["samples"])
    output: 
        total="outputs/1d_merged_bam/dedup_uniq_merged.bam",
        plus="outputs/1d_merged_bam/dedup_uniq_merged_plus.bam",
        minus="outputs/1d_merged_bam/dedup_uniq_merged_minus.bam"
    conda:
        "../envs/module_1.yaml"
    log:
        "logs/1d_merged_bam/merged_before_filt.log"
    shell:
        """
        echo 'Merging the bam files' &> {log}
        samtools merge {output.total} {input} &>> {log}
        echo 'Indexing the merged bam file...' &>> {log}
        samtools index {output.total} &>> {log}

        echo 'Splitting the plus bam file...' &>> {log}
        samtools view -F 16 -b {output.total} > {output.plus} 2>> {log}
        echo 'Splitting the minus bam file...' &>> {log}
        samtools view -f 16 -b {output.total} > {output.minus} 2>> {log}

        echo 'Indexing the plus bam file...' &>> {log}
        samtools index {output.plus} &>> {log}
        echo 'Indexing the minus bam file...' &>> {log}
        samtools index {output.minus} &>> {log}
        """

rule merge_bam_after_filt:
    input: 
        expand("outputs/1c_subset_bam/{sample}_dedup_uniq_filtered_fibroblasts.bam", sample=config["samples"])
    output: 
        total="outputs/1d_merged_bam/dedup_uniq_genomicAfiltered_merged.bam",
        plus="outputs/1d_merged_bam/dedup_uniq_genomicAfiltered_merged_plus.bam",
        minus="outputs/1d_merged_bam/dedup_uniq_genomicAfiltered_merged_minus.bam"
    conda:
        "../envs/module_1.yaml"
    log:
        "logs/1d_merged_bam/merged_after_filt.log"
    shell:
        """
        echo 'Merging the bam files' &> {log}
        samtools merge {output.total} {input} &>> {log}
        echo 'Indexing the merged bam file...' &>> {log}
        samtools index {output.total} &>> {log}

        echo 'Splitting the plus bam file...' &>> {log}
        samtools view -F 16 -b {output.total} > {output.plus} 2>> {log}
        echo 'Splitting the minus bam file...' &>> {log}
        samtools view -f 16 -b {output.total} > {output.minus} 2>> {log}

        echo 'Indexing the plus bam file...' &>> {log}
        samtools index {output.plus} &>> {log}
        echo 'Indexing the minus bam file...' &>> {log}
        samtools index {output.minus} &>> {log}
        """
rule bw_track_after_filt:
    input: 
        plus="outputs/1d_merged_bam/dedup_uniq_genomicAfiltered_merged_plus.bam",
        minus="outputs/1d_merged_bam/dedup_uniq_genomicAfiltered_merged_minus.bam",
    output: 
        plus="outputs/qc_ucsc_tracks/dedup_uniq_genomicAfiltered_merged_plus_raw.bw",
        minus="outputs/qc_ucsc_tracks/dedup_uniq_genomicAfiltered_merged_minus_raw.bw",
    conda:
        "../envs/module_1.yaml"
    log:
        "logs/qc_bw_tracks/qc_bw_tracks_after_filt.log"
    shell:
        """
        mkdir -p outputs/qc_ucsc_tracks/
        bamCoverage -b {input.plus} -o {output.plus} -p 36 -v --binSize 1 &> {log}
        bamCoverage -b {input.minus} -o {output.minus} -p 36 -v --binSize 1 &>> {log}
        """
rule bw_track_before_filt:
    input: 
        plus="outputs/1d_merged_bam/dedup_uniq_merged_plus.bam",
        minus="outputs/1d_merged_bam/dedup_uniq_merged_minus.bam",
    output: 
        plus="outputs/qc_ucsc_tracks/dedup_uniq_merged_plus_raw.bw",
        minus="outputs/qc_ucsc_tracks/dedup_uniq_merged_minus_raw.bw",
    conda:
        "../envs/module_1.yaml"
    log:
        "logs/qc_bw_tracks/qc_bw_tracks_before_filt.log"
    shell:
        """
        mkdir -p outputs/qc_ucsc_tracks/
        bamCoverage -b {input.plus} -o {output.plus} -p 36 -v --binSize 1 &> {log}
        bamCoverage -b {input.minus} -o {output.minus} -p 36 -v --binSize 1 &>> {log}
        """
