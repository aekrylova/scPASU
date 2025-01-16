configfile: "../config/config.yaml"

rule cellranger_make_scPASU_ref:
    input: 
        indir="outputs/3e_fragmented_peaks_to_merge/",
        fasta=f"{config['cellrangerrefpath']}fasta/genome.fa",
        work_dir=config["work_dir"]
    params:
        version=2024,
        PAS_filtering=config["PAS_filtering"],
        ncores=config["ncores"],
        memgb=160,
        compartment=config['compartment'],
        genome=config['genomename'], 
        cellranger_path=config["cellrangerpath"],
    output: 
        directory(f"outputs/{config['genomename']}_{config['compartment']}/")
    conda:
        "../envs/module_4.yaml"
    log:
        "logs/cellranger_make_scPASU_ref.log"
    shell:  
        """
        export PATH={params.cellranger_path}/:$PATH

        ref={input.work_dir}/{input.indir}/{params.compartment}_$(if [ {params.PAS_filtering} == 'true' ]; then continue; else echo noPASfiltering_; fi)final_peak_universe_updated.txt
        gtf={input.work_dir}/{input.indir}/{params.compartment}_$(if [ {params.PAS_filtering} == 'true' ]; then continue; else echo noPASfiltering_; fi)final_peak_universe_updated.gtf

        echo Convert final polyA site reference to gtf file format... &> {log}
        Rscript {input.work_dir}/scripts/ref2gtf.R -b $ref -f {params.compartment} -o {input.indir}/ &>> {log}

        echo Create reference package... &>> {log}
        cd {input.work_dir}/outputs/
        cellranger mkref --ref-version={params.version} --genome={params.genome}_{params.compartment} --fasta={input.work_dir}/{input.fasta} --genes=$gtf --nthreads={params.ncores} --memgb={params.memgb} &>> {input.work_dir}/{log}
        """

rule bam_to_fastq:
    input: 
        wait_dir="outputs/3e_fragmented_peaks_to_merge/",
        work_dir=config["work_dir"]
    output:
        directory("outputs/4a_bam_to_fastq/") 
    params:
        ncores=config["ncores"],
        samples=config["samples"],
    conda:
        "../envs/module_4.yaml"
    log:
        "logs/4a_bam_to_fastq.log"
    shell: 
        """
        set +o pipefail;
        mkdir -p {output}

        echo Convert bam to fastq... &> {log}

        for sample in {params.samples}; do
            bamtofastq --nthreads={params.ncores} {input.work_dir}/outputs/1c_subset_bam/${{sample}}_dedup_uniq_filtered_fibroblasts.bam {input.work_dir}/{output}/${{sample}} &>> {input.work_dir}/{log}
            fastq_dir={input.work_dir}/{output}/${{sample}}/
            sample_dir=$(ls $fastq_dir | grep ${{sample}})
            fastq_dir=${{fastq_dir}}${{sample_dir}}/
            cd $fastq_dir

            if ! ls | grep -q "^bamtofastq"; then 
                echo "Filenames do not need to be modified" &>> {input.work_dir}/{log}
            else 
                for f in $(ls bamtofastq*)
                    do mv $f ${{sample}}${{f#"bamtofastq"}} &>> {input.work_dir}/{log}
                done
            fi
        done
        """

rule cellranger_genecount:
    input: 
        cellranger_ref_path=config["cellrangerrefpath"],
        fastqdir="outputs/4a_bam_to_fastq/",
        work_dir=config["work_dir"]
    output: 
        directory("outputs/4b_cellranger_genecount/")    
    params:
        localcores=16,
        localmem=200,
        samples=config["samples"],
        cellranger_path=config["cellrangerpath"],
    conda:
        "../envs/module_4.yaml"
    log:
        "logs/4b_cellranger_genecount.log"
    shell:  
        """
        export PATH={params.cellranger_path}:$PATH
        mkdir -p {output} &> {log}
        cd {output}/
        for sample in {params.samples}; do
            fastq_dir={input.work_dir}/{input.fastqdir}/${{sample}}/
            sample_dir=$(ls $fastq_dir | grep $sample)
            fastq_dir=${{fastq_dir}}${{sample_dir}} 
            cellranger count --id $sample --fastqs=${{fastq_dir}} --sample=${{sample}} --localcores={params.localcores} --localmem={params.localmem} --transcriptome={input.work_dir}/{input.cellranger_ref_path}/ &>> {input.work_dir}/{log}
        done
        """

rule cellranger_peakcount:
    input: 
        peak_refdir=f"outputs/{config['genomename']}_{config['compartment']}/",
        fastqdir="outputs/4a_bam_to_fastq/",
        work_dir=config["work_dir"]
    output: 
        directory("outputs/4b_cellranger_peakcount/")    
    params:
        PAS_filtering=config["PAS_filtering"],
        localcores=16,
        localmem=200,
        samples=config["samples"],
        compartment=config["compartment"],
        genome=config['genomename'], 
        cellranger_path=config["cellrangerpath"],
    conda:
        "../envs/module_4.yaml"
    log:
        "logs/4b_cellranger_peakcount.log"
    shell:  
        """
        if [ {params.PAS_filtering} == 'true' ]
            then fprefix={params.compartment}
        else
            fprefix={params.compartment}_noPASfiltering
        fi

        genome={params.genome}_${{fprefix}}
        peak_refdir={input.work_dir}/outputs/${{genome}}

        export PATH={params.cellranger_path}:$PATH
        mkdir -p {output} &> {log}
        cd {output}/
        for sample in {params.samples}; do
            fastq_dir={input.work_dir}/{input.fastqdir}/${{sample}}/
            sample_dir=$(ls $fastq_dir | grep $sample)
            fastq_dir=${{fastq_dir}}${{sample_dir}} 
            cellranger count --id $sample --fastqs=${{fastq_dir}} --sample=${{sample}} --localcores={params.localcores} --localmem={params.localmem} --transcriptome={input.work_dir}/{input.peak_refdir}/ &>> {input.work_dir}/{log}
        done
        """