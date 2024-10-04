rule create_links:
    params:
        maf_source_dir=maf_source_dir,
        bam_source_dir=bam_source_dir,
        insert_size_source_dir=insert_size_source_dir,
        alignment_dir =alignment_dir,
        maf_dir = maf_dir,
        insert_dir=insert_dir,
        work_dir= work_dir,
        symlink_dir = symlink_dir
    output:
        bams = expand(os.path.join(work_dir, symlink_dir, alignment_dir, "{sample}.bam"), sample=sample_list),
        bais = expand(os.path.join(work_dir, symlink_dir, alignment_dir, "{sample}.bam.bai"), sample=sample_list),
        mafs = expand(os.path.join(work_dir, symlink_dir, maf_dir, "{sample}.maf"), sample=sample_list),
        insert_size = expand(os.path.join(work_dir, symlink_dir, insert_dir, "{sample}_insert-size.txt"), sample=sample_list),
    threads: 1
    container: None
    log: os.path.join(work_dir, log_dir, "Create_Symlinks.log")
    resources:
        mem_mb=1000,
        runtime=20,
        nodes=1,
    run:
        bam_target_dir = os.path.join(params.work_dir, params.symlink_dir, params.alignment_dir)
        os.makedirs(bam_target_dir, exist_ok=True)
        maf_target_dir = os.path.join(params.work_dir, params.symlink_dir, params.maf_dir)
        os.makedirs(maf_target_dir, exist_ok=True)
        insert_size_target_dir = os.path.join(params.work_dir, params.symlink_dir, params.insert_dir)
        os.makedirs(insert_size_target_dir, exist_ok=True)
        # BAM
        for bam in os.scandir(params.bam_source_dir):
            if not (bam.name.endswith(".bam") or bam.name.endswith(".bai")): continue
            sample_id = bam.name.split("_")[0]
            tumor_type = bam.name.split("_")[1].split(".")[0]
            bam_source = os.readlink(bam.path)
            if bam.name.endswith(".bam"):
                bam_target = os.path.join(bam_target_dir, f"{sample_id}_{tumor_type}.bam")
            elif bam.name.endswith(".bai"):
                bam_target = os.path.join(bam_target_dir, f"{sample_id}_{tumor_type}.bam.bai")
            if os.path.exists(bam_target): os.remove(bam_target)
            os.symlink(bam_source, bam_target)
        # MAF
        for maf in os.scandir(params.maf_source_dir):
            if not (maf.name.endswith("maf")): continue
            maf_source = os.readlink(maf.path)
            sample_id = maf.name.split("_")[0]
            tumor_type = maf.name.split("_")[1].split(".")[0]
            maf_target = os.path.join(maf_target_dir, f"{sample_id}_{tumor_type}.maf")
            if os.path.exists(maf_target): os.remove(maf_target)
            os.symlink(maf_source, maf_target)
        # Insert Size
        for insert in os.scandir(params.insert_size_source_dir):
            # print(insert.path)
            if not insert.name.endswith("insert-size.txt"): continue
            insert_source = os.readlink(insert.path)
            sample_id = insert.name.split("_")[0]
            tumor_type = insert.name.split("_")[1].split(".")[0]
            insert_target = os.path.join(insert_size_target_dir, f"{sample_id}_{tumor_type}_insert-size.txt")
            # print(insert_target.path)
            # print(t)
            if os.path.exists(insert_target): os.remove(insert_target)
            os.symlink(insert_source, insert_target)

rule gen_occ:
    input:
        genome=genome
    params:
        blat_binary=blat_binary,
        blat_temp_out = os.path.join(work_dir, ref_data_dir, os.path.split(genome)[1].replace("fa", "psl")),
        output_directory = os.path.join(work_dir, ref_data_dir),
        file_name = "11.occ",
    output:
        os.path.join(work_dir, ref_data_dir, "11.occ")
    threads: 1
    container: None
    log: os.path.join(work_dir, log_dir, "Gen_Occ.log")
    resources:
        mem_mb=4000,
        runtime=20,
        nodes=1,
    shell:
        "cd {params.output_directory}; "
        "{params.blat_binary} "
        "{input.genome} " # Database
        "{input.genome} " # Query
        "-makeOoc={params.file_name} "
        "-t=dna -q=dna "
        "{params.blat_temp_out} "
        "&> {log}"

rule gen_2bit:
    input:
        genome=genome
    params:
        faToTwoBit_binary =faToTwoBit_binary,
    output:
        os.path.join(work_dir, ref_data_dir,  os.path.split(genome)[1].replace("fa", "2bit"))
    threads: 1
    container: None
    log: os.path.join(work_dir, log_dir, "Gen_2Bit_genome.log")
    resources:
        mem_mb=4000,
        runtime=20,
        nodes=1,
    shell:
        "{params.faToTwoBit_binary} "
        "{input.genome} "
        "{output} "
        "&> {log}"
