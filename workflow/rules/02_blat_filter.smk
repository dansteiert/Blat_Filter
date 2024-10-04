rule scatter_maf:
    input:
        maf = os.path.join(work_dir, symlink_dir, maf_dir, "{sample}.maf")
    params:
        scatter_count = scatter_count
    output:
        maf_tmp = expand(os.path.join(work_dir, tmp_blat_dir, maf_dir,  "{{sample}}", "{{sample}}_{scatter}.maf"), scatter=[str(i) for i in range(0, scatter_count)])
    threads: 1
    container: None
    log: os.path.join(work_dir, log_dir, "{sample}_scatter_maf.log")
    resources:
        mem_mb=4000,
        runtime=60,
        nodes=1,
    run:
        ## Open Maf file and find header rows
        header_rows = []
        f = open(input.maf, "r")
        header_idx = 0
        while True:
            line = f.readline()
            if line.startswith("Hugo_Symbol"):
                # header_rows.append(line)
                break
            header_idx += 1
            header_rows.append(line)
            if header_idx > 2000:
                break
        if header_idx > 2000:
            print("Over 2000")
            print(header_rows)
        # Calculate # of Chunks
        with open(input.maf, "rb") as f:
            num_lines = sum(1 for _ in f)
        num_lines -= header_idx
        chunksize = np.ceil(num_lines/params.scatter_count)
        ## Writeout maf in smaller chunks
        header_string = "".join(header_rows)
        os.makedirs(os.path.split(output.maf_tmp[0])[0], exist_ok=True)
        scatter_index = 0
        for df_maf in pd.read_csv(input.maf, sep="\t", skiprows=header_idx, low_memory=False, chunksize=chunksize):
            df_out = df_maf.to_csv(index=False, sep='\t')
            with open(os.path.join(os.path.split(output.maf_tmp[0])[0], f"{wildcards.sample}_{scatter_index}.maf"), 'w') as f:
                f.write(header_string)
                f.write(df_out)
            scatter_index +=1




rule filter_maf:
    input:
        maf=os.path.join(work_dir, tmp_blat_dir, maf_dir, "{sample}", "{sample}_{scatter}.maf"),
        bam=os.path.join(work_dir, symlink_dir, alignment_dir, "{sample}.bam"),
        insert_size=os.path.join(work_dir, symlink_dir, insert_dir, "{sample}_insert-size.txt"),
        database=os.path.join(work_dir, ref_data_dir,  os.path.split(genome)[1].replace("fa", "2bit")),
        occ=os.path.join(work_dir, ref_data_dir, "11.occ"),
    params:
        blat_binary =blat_binary,
        output_prefix = os.path.join(work_dir, tmp_blat_dir, filtered_dir, "{sample}", "{sample}_{scatter}"),
        blat_script = os.path.join(blat_filter_script_dir, "src", "blat_filter.py")
    output:
        maf = os.path.join(work_dir, tmp_blat_dir, filtered_dir, "{sample}", "{sample}_{scatter}.blat.passed.maf"),
        maf_all = os.path.join(work_dir, tmp_blat_dir, filtered_dir, "{sample}", "{sample}_{scatter}.blat.all.maf"),
        maf_rejected = os.path.join(work_dir, tmp_blat_dir, filtered_dir, "{sample}", "{sample}_{scatter}.blat.rejected.maf"),
    threads: 1
    conda: "../envs/Filter_Blat.yaml"
    # container: blat_filter_container
    log: os.path.join(work_dir, log_dir, "{sample}_{scatter}_blat_filter.log")
    resources:
        mem_mb=4000,
        runtime=60*24 * 1,
        nodes=1,
    script:
        "../scripts/blat_filter.py"
        ## PairName
        # "echo $SHELL; cat /proc/$$/cmdline; python -V; "
        #"python /app/blat_filter.py " #{params.blat_script} "
        #"--blat_executable /app/blat " #"{params.blat_binary} "
        #"--bam {input.bam} "
        #"--maf {input.maf} "
        #"--database {input.database} "
        #"--ooc {input.occ} "
        #"--out_name {params.output_prefix} " 
        #"--maf_chr_col Chromosome "
        #"--maf_chr_col_to_score_reads Chromosome " 
        #"--maf_start_pos_col Start_Position "
        #"--maf_start_pos_col_to_score_reads Start_Position "
        #"--maf_end_pos_col End_Position " 
        #"--maf_t_ref_count_col t_ref_count "
        #"--maf_t_alt_count_col t_alt_count "
        #"--set_maf_col_manually " 
        #"--insert_size_metrics {input.insert_size} "
        # "&> {log}"



rule count_variants_maf:
    input:
        maf = os.path.join(work_dir, tmp_blat_dir, filtered_dir, "{sample}", "{sample}_{scatter}.blat.{group}.maf"),
    output:
        txt = os.path.join(work_dir, tmp_blat_dir, filtered_dir, "{sample}", "{sample}_{scatter}.blat.{group}_mutation_counts.txt"),
    params:
        count_variants_script = os.path.join(blat_filter_script_dir, "src", "count_variants.py")
    threads: 1
    # container: blat_filter_container
    conda: "../envs/Filter_Blat.yaml"
    log: os.path.join(work_dir, log_dir, "{sample}_{scatter}_{group}_mutation_counts.log")
    resources:
        mem_mb=4000,
        runtime=20,
        nodes=1,
    shell:
        ## PairName
        #"python /app/count_variants.py "#{params.count_variants_script}"
        #"{input.maf} {output.maf} "
        # "&> {log}"
        "../scripts/count_variants.py"





rule gather_maf:
    input:
        mafs = expand(os.path.join(work_dir, tmp_blat_dir, filtered_dir, "{{sample}}", "{{sample}}_{scatter}.blat.{{group}}.maf"), scatter=[str(i) for i in range(0, scatter_count)])
    params:
        chunksize=10000
    output:
        maf = os.path.join(work_dir, blat_final_dir, "{sample}.blat.{group}.maf")
    threads: 1
    container: None
    log: os.path.join(work_dir, log_dir, "{sample}_{group}_gather_maf.log")
    resources:
        mem_mb=4000,
        runtime=60,
        nodes=1,
    run:
        header = True        
        for maf in input.mafs:
            ## Open Maf file and find header rows
            header_rows = []
            f = open(maf, "r")
            header_idx = 0
            while True:
                line = f.readline()
                if line.startswith("Hugo_Symbol"):
                    header_rows.append(line)
                    break
                header_idx += 1
                header_rows.append(line)
                if header_idx > 2000:
                    break
            if header_idx > 2000:
                print("Over 2000")
                print(header_rows)
            header_string = "".join(header_rows)
            os.makedirs(os.path.split(output.maf_tmp)[0], exist_ok=True)
            for df_maf in pd.read_csv(maf, sep="\t", skiprows=header_idx, low_memory=False, chunksize=params.chunksize):
                df_out = df_maf.to_csv(index=False, sep='\t')
                if header:
                    with open(output.maf, 'w') as f:
                        f.write(header_string)
                        f.write(df_out)
                    header=False
                else:
                    df_maf.to_csv(maf, index=False, sep="\t", header=header, mode="a")
            
            
