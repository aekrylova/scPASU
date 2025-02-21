configfile: "../config/config.yaml"

rule all_filtered_reads_peak_call_summits:
    input:
        plus="outputs/1d_merged_bam/dedup_uniq_genomicAfiltered_merged_plus.bam",
        minus="outputs/1d_merged_bam/dedup_uniq_genomicAfiltered_merged_minus.bam"
    output:
        directory("outputs/2a_peaks/all_filtered_reads/")
    params:
        gsize="hs",
        extsize=200
    conda:
        "../envs/module_2.yaml"
    log:
        "logs/2a_all_filtered_reads_peak_call_summits.log"
    shell:
        """
        bam_prefix_plus=dedup_uniq_genomicAfiltered_merged_plus
        macs2 callpeak \\
            --seed 123 \\
            -t {input.plus} \\
            --keep-dup all \\
            --gsize {params.gsize} \\
            --format BAM \\
            --nomodel \\
            --extsize {params.extsize} \\
            --outdir {output} \\
            --name $bam_prefix_plus \\
            -B \\
            --verbose 2 \\
            --call-summits &> {log}
        
        bam_prefix_minus=dedup_uniq_genomicAfiltered_merged_minus
        macs2 callpeak \\
            --seed 123 \\
            -t {input.minus} \\
            --keep-dup all \\
            --gsize {params.gsize} \\
            --format BAM \\
            --nomodel \\
            --extsize {params.extsize} \\
            --outdir {output} \\
            --name $bam_prefix_minus \\
            -B \\
            --verbose 2 \\
            --call-summits &>> {log}
        """

rule polyAreads_sam_to_polyA:
    input: 
        plus="outputs/1d_merged_bam/dedup_uniq_merged_plus.bam",
        minus="outputs/1d_merged_bam/dedup_uniq_merged_minus.bam",
        work_dir=config["work_dir"]
    output: 
        plus="outputs/2b_polya/before_genomicAfiltering/polyAreads_beforegenomicAfiltering_plus.bed",
        minus="outputs/2b_polya/before_genomicAfiltering/polyAreads_beforegenomicAfiltering_minus.bed",
        plus_nws="outputs/2b_polya/before_genomicAfiltering/polyAreads_beforegenomicAfiltering_plus_nowrongstrand.bed",
        minus_nws="outputs/2b_polya/before_genomicAfiltering/polyAreads_beforegenomicAfiltering_minus_nowrongstrand.bed"
    conda:
        "../envs/module_2.yaml"
    params:
        min_clipped=9,
        min_a_content=0.85
    log:
        "logs/2b_polyAreads_sam_to_polyA.log"
    shell: 
        """
        echo Finding polyA sites in plus.bam and writing output to plus.bed &> {log}
        samtools view {input.plus} | {input.work_dir}/scripts/2_peak_ref/2a_samToPolyA.pl --minClipped={params.min_clipped} --minAcontent={params.min_a_content} - > {output.plus} 2>> {log}

        echo Remove wrong-strand polyA reads &>> {log}
        grep -w '+' {output.plus} > {output.plus_nws} 2>> {log}
        
        echo Finding polyA sites in minus.bam and writing output to minus.bed &>> {log}
        samtools view {input.minus} | {input.work_dir}/scripts/2_peak_ref/2a_samToPolyA.pl --minClipped={params.min_clipped} --minAcontent={params.min_a_content} - > {output.minus} 2>> {log}

        echo Remove wrong-strand polyA reads &>> {log}
        grep -w '-' {output.minus} > {output.minus_nws} 2>> {log}
        """ 

rule all_filtered_reads_sam_to_polyA:
    input: 
        plus="outputs/1d_merged_bam/dedup_uniq_genomicAfiltered_merged_plus.bam",
        minus="outputs/1d_merged_bam/dedup_uniq_genomicAfiltered_merged_minus.bam",
        work_dir=config["work_dir"]
    output: 
        total=directory("outputs/2b_polya/after_genomicAfiltering/"),
        plus="outputs/2b_polya/after_genomicAfiltering/polyAreads_aftergenomicAfiltering_plus.bed",
        minus="outputs/2b_polya/after_genomicAfiltering/polyAreads_aftergenomicAfiltering_minus.bed",
    conda:
        "../envs/module_2.yaml"
    params:
        min_clipped=9,
        min_a_content=0.85
    log:
        "logs/2b_all_filtered_reads_sam_to_polyA.log"
    shell: 
        """
        echo Finding polyA sites in plus.bam and writing output to plus.bed &> {log}
        samtools view {input.plus} | {input.work_dir}/scripts/2_peak_ref/2a_samToPolyA.pl --minClipped={params.min_clipped} --minAcontent={params.min_a_content} - > {output.plus} 2>> {log}

        echo Remove wrong-strand polyA reads &>> {log}
        polya_file_raw={output.plus}
        polya_file=${{polya_file_raw%.bed}}_nowrongstrand.bed
        grep -w '+' $polya_file_raw > $polya_file 2>> {log}
        
        echo Finding polyA sites in minus.bam and writing output to minus.bed  &>> {log}
        samtools view {input.minus} | {input.work_dir}/scripts/2_peak_ref/2a_samToPolyA.pl --minClipped={params.min_clipped} --minAcontent={params.min_a_content} - > {output.minus} 2>> {log}

        echo Remove wrong-strand polyA reads &>> {log}
        polya_file_raw={output.minus}
        polya_file=${{polya_file_raw%.bed}}_nowrongstrand.bed
        grep -w '-' $polya_file_raw > $polya_file 2>> {log}
        """  

rule polyAreads_filterbam_byreadnames:
    input: 
        plus="outputs/1d_merged_bam/dedup_uniq_merged_plus.bam",
        minus="outputs/1d_merged_bam/dedup_uniq_merged_minus.bam",
        polya_plus="outputs/2b_polya/before_genomicAfiltering/polyAreads_beforegenomicAfiltering_plus_nowrongstrand.bed",
        polya_minus="outputs/2b_polya/before_genomicAfiltering/polyAreads_beforegenomicAfiltering_minus_nowrongstrand.bed"
    output:
        readname_minus="outputs/2b_polya/before_genomicAfiltering/read_names/polyAreads_beforegenomicAfiltering_minus_nowrongstrand.csv",
        bam_minus="outputs/2b_polya/before_genomicAfiltering/bam_files/polyAreads_beforegenomicAfiltering_minus_nowrongstrand.bam",
        readname_plus="outputs/2b_polya/before_genomicAfiltering/read_names/polyAreads_beforegenomicAfiltering_plus_nowrongstrand.csv",
        bam_plus="outputs/2b_polya/before_genomicAfiltering/bam_files/polyAreads_beforegenomicAfiltering_plus_nowrongstrand.bam"
    conda:
        "../envs/module_2.yaml"
    log:
        "logs/2b_polyAreads_filterbam_byreadnames.log"
    shell:  
        """
        mkdir -p outputs/2b_polya/before_genomicAfiltering/read_names/
        mkdir -p outputs/2b_polya/before_genomicAfiltering/bam_files/

        cut -f4 {input.polya_plus} > {output.readname_plus} 2> {log}
        
        echo Filter sam reads by read names... &>> {log}
        picard FilterSamReads I={input.plus} O={output.bam_plus} READ_LIST_FILE={output.readname_plus} FILTER=includeReadList &>> {log}
   
        echo Indexing the last bam file... &>> {log}
        samtools index {output.bam_plus} &>> {log}

        cut -f4 {input.polya_minus} > {output.readname_minus} 2>> {log}

        echo Filter sam reads by read names... &>> {log}
        picard FilterSamReads I={input.minus} O={output.bam_minus} READ_LIST_FILE={output.readname_minus} FILTER=includeReadList &>> {log}
   
        echo Indexing the last bam file... &>> {log}
        samtools index {output.bam_minus} &>> {log}
        """

rule polyAreads_peak_call_summits:
    input: 
        minus="outputs/2b_polya/before_genomicAfiltering/bam_files/polyAreads_beforegenomicAfiltering_minus_nowrongstrand.bam",
        plus="outputs/2b_polya/before_genomicAfiltering/bam_files/polyAreads_beforegenomicAfiltering_plus_nowrongstrand.bam"
    output: 
        directory("outputs/2a_peaks/polyA_reads/")
    conda:
        "../envs/module_2.yaml"
    log:
        "logs/2a_polyAreads_peak_call_summits.log"
    params:
        gsize="hs",
        extsize=200
    shell:
        """
        bam_prefix_plus=polyAreads_beforegenomicAfiltering_plus_nowrongstrand
        macs2 callpeak \\
            --seed 123 \\
            -t {input.plus} \\
            --keep-dup all \\
            --gsize {params.gsize} \\
            --format BAM \\
            --nomodel \\
            --extsize {params.extsize} \\
            --outdir {output} \\
            --name $bam_prefix_plus \\
            -B \\
            --verbose 2 \\
            --call-summits &> {log}
        
        bam_prefix_minus=polyAreads_beforegenomicAfiltering_minus_nowrongstrand
        macs2 callpeak \\
            --seed 123 \\
            -t {input.minus} \\
            --keep-dup all \\
            --gsize {params.gsize} \\
            --format BAM \\
            --nomodel \\
            --extsize {params.extsize} \\
            --outdir {output} \\
            --name $bam_prefix_minus \\
            -B \\
            --verbose 2 \\
            --call-summits &>> {log}
        """ 

rule bam_coverage_all_filtered_reads:
    input: 
        summit_dir="outputs/2a_peaks/all_filtered_reads/",
        bam_plus="outputs/1d_merged_bam/dedup_uniq_genomicAfiltered_merged_plus.bam",
        bam_minus="outputs/1d_merged_bam/dedup_uniq_genomicAfiltered_merged_minus.bam",
        work_dir=config["work_dir"]
    params:
        binsize=1,
        fragmentsize=200,
        ncores=config["ncores"],
    output:
        bw_plus="outputs/ucsc_tracks/bigwig_tracks/1bp/all_filtered_reads/dedup_uniq_genomicAfiltered_merged_plus_ext.bw",
        bw_minus="outputs/ucsc_tracks/bigwig_tracks/1bp/all_filtered_reads/dedup_uniq_genomicAfiltered_merged_minus_ext.bw",
        bw_plus_raw="outputs/ucsc_tracks/bigwig_tracks/1bp/all_filtered_reads/dedup_uniq_genomicAfiltered_merged_plus_raw.bw",
        bw_minus_raw="outputs/ucsc_tracks/bigwig_tracks/1bp/all_filtered_reads/dedup_uniq_genomicAfiltered_merged_minus_raw.bw",
        bg_plus="outputs/2c1_alignment_signals/all_filtered_reads/dedup_uniq_genomicAfiltered_merged_plus.bedGraph",
        bg_minus="outputs/2c1_alignment_signals/all_filtered_reads/dedup_uniq_genomicAfiltered_merged_minus.bedGraph",
        bg_dir=directory("outputs/2c1_alignment_signals/all_filtered_reads/")
    conda:
        "../envs/module_2.yaml"
    log:
        "logs/2c1_bam_coverage/all_filtered_reads.log" 
    shell: 
        """
        bamCoverage -b {input.bam_plus} -o {output.bw_plus} -p {params.ncores} -v --binSize {params.binsize} --extendReads {params.fragmentsize} &> {log}
        bigWigToBedGraph {output.bw_plus} {output.bg_plus} &>> {log}
        echo Also creating raw big Wig track for track visualization... &>> {log}
        bamCoverage -b {input.bam_plus} -o {output.bw_plus_raw} -p {params.ncores} -v --binSize {params.binsize} &>> {log}

        bamCoverage -b {input.bam_minus} -o {output.bw_minus} -p {params.ncores} -v --binSize {params.binsize} --extendReads {params.fragmentsize} &>> {log}
        bigWigToBedGraph {output.bw_minus} {output.bg_minus} &>> {log}
        echo Also creating raw big Wig track for track visualization... &>> {log}
        bamCoverage -b {input.bam_minus} -o {output.bw_minus_raw} -p {params.ncores} -v --binSize {params.binsize} &>> {log}

        cd {input.summit_dir}
        plus_summit=$(ls *plus*peaks.xls)
        grep -v '^#\|^[[:space:]]*$' ${{plus_summit}} > ${{plus_summit%_peaks.xls}}_peaks_summits.txt 2>> {input.work_dir}/{log}
        minus_summit=$(ls *minus*peaks.xls) 
        grep -v '^#\|^[[:space:]]*$' ${{minus_summit}} > ${{minus_summit%_peaks.xls}}_peaks_summits.txt 2>> {input.work_dir}/{log}
        """

rule bam_coverage_polyA_reads:
    input: 
        summit_dir="outputs/2a_peaks/polyA_reads/",
        bam_plus="outputs/2b_polya/before_genomicAfiltering/bam_files/polyAreads_beforegenomicAfiltering_plus_nowrongstrand.bam",
        bam_minus="outputs/2b_polya/before_genomicAfiltering/bam_files/polyAreads_beforegenomicAfiltering_minus_nowrongstrand.bam",
        work_dir=config["work_dir"]
    params:
        binsize=1,
        fragmentsize=200,
        ncores=config["ncores"],
    output:
        bw_plus="outputs/ucsc_tracks/bigwig_tracks/1bp/polyA_reads/polyAreads_beforegenomicAfiltering_plus_nowrongstrand_ext.bw",
        bw_minus="outputs/ucsc_tracks/bigwig_tracks/1bp/polyA_reads/polyAreads_beforegenomicAfiltering_minus_nowrongstrand_ext.bw",
        bw_plus_raw="outputs/ucsc_tracks/bigwig_tracks/1bp/polyA_reads/polyAreads_beforegenomicAfiltering_plus_nowrongstrand_raw.bw",
        bw_minus_raw="outputs/ucsc_tracks/bigwig_tracks/1bp/polyA_reads/polyAreads_beforegenomicAfiltering_minus_nowrongstrand_raw.bw",
        bg_plus="outputs/2c1_alignment_signals/polyA_reads/polyAreads_beforegenomicAfiltering_plus_nowrongstrand.bedGraph",
        bg_minus="outputs/2c1_alignment_signals/polyA_reads/polyAreads_beforegenomicAfiltering_minus_nowrongstrand.bedGraph",
        bg_dir=directory("outputs/2c1_alignment_signals/polyA_reads/")
    conda:
        "../envs/module_2.yaml"
    log:
        "logs/2c1_bam_coverage/polyA_reads.log" 
    shell: 
        """
        bamCoverage -b {input.bam_plus} -o {output.bw_plus} -p {params.ncores} -v --binSize {params.binsize} --extendReads {params.fragmentsize} &> {log}
        bigWigToBedGraph {output.bw_plus} {output.bg_plus} &>> {log}
        echo Also creating raw big Wig track for track visualization... &>> {log}
        bamCoverage -b {input.bam_plus} -o {output.bw_plus_raw} -p {params.ncores} -v --binSize {params.binsize} &>> {log}

        bamCoverage -b {input.bam_minus} -o {output.bw_minus} -p {params.ncores} -v --binSize {params.binsize} --extendReads {params.fragmentsize} &>> {log}
        bigWigToBedGraph {output.bw_minus} {output.bg_minus} &>> {log}
        echo Also creating raw big Wig track for track visualization... &>> {log}
        bamCoverage -b {input.bam_minus} -o {output.bw_minus_raw} -p {params.ncores} -v --binSize {params.binsize} &>> {log}

        cd {input.summit_dir}
        plus_summit=$(ls *plus*peaks.xls)
        grep -v '^#\|^[[:space:]]*$' ${{plus_summit}} > ${{plus_summit%_peaks.xls}}_peaks_summits.txt 2>> {input.work_dir}/{log}
        minus_summit=$(ls *minus*peaks.xls) 
        grep -v '^#\|^[[:space:]]*$' ${{minus_summit}} > ${{minus_summit%_peaks.xls}}_peaks_summits.txt 2>> {input.work_dir}/{log}
        """

rule split_peak_ref_all_filtered_reads:
    input: 
        summit_dir="outputs/2a_peaks/all_filtered_reads/",
        bg_dir="outputs/2c1_alignment_signals/all_filtered_reads/",
        work_dir=config["work_dir"],
    params:
        compartment=config["compartment"],
        extsize=90,
        chrs=config["chrs"],
        ncores=config["ncores"],
        direction=3,
        goldmine_path=config["goldminepath"]
    output:
        directory("outputs/2c2_split_peaks/all_filtered_reads/"),
    conda:
        "../envs/module_2.yaml"
    log:
        "logs/2c2_split_peak_ref/all_filtered_reads.log" 
    shell: 
        """
        fprefix=/{params.compartment}_all_filtered_reads
        Rscript {input.work_dir}/scripts/2_peak_ref/2c_split_peaks.R -n {params.ncores} -p all_filtered_reads -b {input.work_dir}/{input.bg_dir}/ -s {input.work_dir}/{input.summit_dir}/ -m {params.goldmine_path} -f $fprefix -o {input.work_dir}/{output}/ -e {params.direction},{params.extsize} -c {params.chrs} &> {log}
        """

rule split_peak_ref_polyA_reads:
    input: 
        summit_dir="outputs/2a_peaks/polyA_reads/",
        bg_dir="outputs/2c1_alignment_signals/polyA_reads/",
        work_dir=config["work_dir"]
    params:
        compartment=config["compartment"],
        extsize=90,
        chrs=config["chrs"],
        ncores=config["ncores"],
        direction=5,
        goldmine_path=config["goldminepath"]
    output:
        directory("outputs/2c2_split_peaks/polyA_reads/"),
    conda:
        "../envs/module_2.yaml"
    log:
        "logs/2c2_split_peak_ref/polyA_reads.log" 
    shell: 
        """
        fprefix=/{params.compartment}_polyA_reads
        Rscript {input.work_dir}/scripts/2_peak_ref/2c_split_peaks.R -n {params.ncores} -p polyA_reads -b {input.work_dir}/{input.bg_dir}/ -s {input.work_dir}/{input.summit_dir}/ -m {params.goldmine_path} -f $fprefix -o {input.work_dir}/{output}/ -e {params.direction},{params.extsize} -c {params.chrs} &> {log}
        """

rule peak_ref_all_filtered_reads:
    input:
        macs_dir="outputs/2c2_split_peaks/all_filtered_reads/",
        polyA_dir="outputs/2b_polya/after_genomicAfiltering/",
        work_dir=config["work_dir"]
    params:
        min_polyA=3,
        chrs=config["chrs"],
    output:
        directory("outputs/2d_intersect_run/all_filtered_reads/"),
    conda:
        "../envs/module_2.yaml"
    log:
        "logs/2d_peak_ref/all_filtered_reads.log" 
    shell: 
        """
        if [ ! -d {input.work_dir}/{output} ]; then
            mkdir -p {input.work_dir}/{output}
        fi

        echo Retain peaks on plus strand with at least {params.min_polyA} polyA junction reads/sites... &> {log}
        bash {input.work_dir}/scripts/2_peak_ref/2d_create_peak_ref.sh -m {input.work_dir}/{input.macs_dir} -p {input.work_dir}/{input.polyA_dir} -o {input.work_dir}/{output} -s plus -c {params.chrs} -n {params.min_polyA} -t {input.work_dir}/scripts/2_peak_ref/ &>> {log}

        echo Retain peaks on minus strand with at least {params.min_polyA} polyA junction reads/sites... &>> {log}
        bash {input.work_dir}/scripts/2_peak_ref/2d_create_peak_ref.sh -m {input.work_dir}/{input.macs_dir} -p {input.work_dir}/{input.polyA_dir} -o {input.work_dir}/{output} -s minus -c {params.chrs} -n {params.min_polyA} -t {input.work_dir}/scripts/2_peak_ref/ &>> {log}
        """ 

rule peak_ref_polyA_reads:
    input:
        macs_dir="outputs/2c2_split_peaks/polyA_reads/",
        work_dir=config["work_dir"]
    params:
        min_polyA=1,
        chrs=config["chrs"],
        polyA_dir="outputs/2b_polya/before_genomicAfiltering/"
    output:
        directory("outputs/2d_intersect_run/polyA_reads/"),
    conda:
        "../envs/module_2.yaml"
    log:
        "logs/2d_peak_ref/polyA_reads.log" 
    shell: 
        """
        if [ ! -d {input.work_dir}/{output} ]; then
            mkdir -p {input.work_dir}/{output}
        fi

        echo Retain peaks on plus strand with at least {params.min_polyA} polyA junction reads/sites... &> {log}
        bash {input.work_dir}/scripts/2_peak_ref/2d_create_peak_ref.sh -m {input.work_dir}/{input.macs_dir} -p {input.work_dir}/{params.polyA_dir} -o {input.work_dir}/{output} -s plus -c {params.chrs} -n {params.min_polyA} -t {input.work_dir}/scripts/2_peak_ref/ &>> {log}

        echo Retain peaks on minus strand with at least {params.min_polyA} polyA junction reads/sites... &>> {log}
        bash {input.work_dir}/scripts/2_peak_ref/2d_create_peak_ref.sh -m {input.work_dir}/{input.macs_dir} -p {input.work_dir}/{params.polyA_dir} -o {input.work_dir}/{output} -s minus -c {params.chrs} -n {params.min_polyA} -t {input.work_dir}/scripts/2_peak_ref/ &>> {log}
        """ 

      
