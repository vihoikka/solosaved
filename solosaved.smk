'''
A snakemake pipeline to search for solosaveds in phage, plasmid and Uniprot proteins.
'''

project = "run_2025_1"
base_path = "/home/vhoikkal/scratch/private/runs/solosaved" + "/" + project

### INPUTS ###

# Databases
carfsaved_hmm_db_path = "/home/vhoikkal/projects/uosa/Ville_Hoikkala/old_cluster/01_transfer/carfsaved/carfsaved.hmm"

# Uniprot
uniprot_db_path = "/mnt/shared/datasets/databases/uniprot/uniref/uniref90/uniref90.fasta.gz"

# Plasmids
plsdb_proteins = "/home/vhoikkal/projects/uosa/Ville_Hoikkala/old_cluster/01_transfer/plsdb/2024_05_31/plsdb_proteins_2024_05_31.faa"

# Phages
millard_fa = "/home/vhoikkal/projects/uosa/Ville_Hoikkala/old_cluster/01_transfer/millard_phages/Mar2025/2Mar2025_genomes.fa"
millard_proteins = "/home/vhoikkal/projects/uosa/Ville_Hoikkala/old_cluster/01_transfer/millard_phages/Mar2025/2Mar2025_vConTACT2_proteins.faa"
millard_metadata = "/home/vhoikkal/projects/uosa/Ville_Hoikkala/old_cluster/01_transfer/millard_phages/Mar2025/2Mar2025_data.tsv"

### PARAMETERS ###
solosaved_max_length = 250
solosaved_e_value = 1e-8

thread_hogger = 60
thread_med = 20
thread_small = 5
thread_min = 1

### CONFIG ###
#PLASMIDS = [os.path.basename(fasta).replace('.fna', '') for fasta in glob.glob(input_plasmid_nt_split_folder + "/*.fna")]

### PATHS ###
# Uniprot
output_folder_uniprot_subset = base_path + "/011_uniprot_small_proteins"
output_folder_uniprot_hmm = base_path + "/012_uniprot_hmm"

# Plasmids
output_folder_plasmids = base_path + "/021_annotated_plasmids"
output_folder_plasmids_hmm = base_path + "/022_plasmids_hmm"
output_folder_plasmid_locus_viz = base_path + "/023_locus_viz"

# Phages
output_folder_phage_genomes = base_path + "/031_phage_genomes"
output_folder_phage_hmm = base_path + "/032_phage_hmm"


### RULES ###

rule all:
    input: output_folder_uniprot_hmm + "/uniprot_solosaved_out.tsv"

### UNIPROT PROTEINS ###

rule uniprot_subset_small_proteins:
    '''
    Subsetting the uniprot db to only include proteins smaller than the max length. Using seqkit
    '''
    input: uniprot_db_path
    output: output_folder_uniprot_subset + "/uniprot_small.faa"
    params: max_length = solosaved_max_length
    threads: thread_hogger
    conda:
        "envs/hmmer.yaml"
    shell:
        '''
        seqkit seq -j {threads} -M {params.max_length} {input} > {output}
        '''
    

rule uniprot_hmmscan_solosaved:
    '''
    Uses carfsaved hmm database to search for solosaveds in small uniprot proteins.
    '''
    input:
        db = carfsaved_hmm_db_path,
        query = rules.uniprot_subset_small_proteins.output
    output:
        hmm_rows = output_folder_uniprot_hmm + "/uniprot_solosaved_out.tsv",
        contig_proteins = output_folder_uniprot_hmm + "/uniprot_solosaved_proteins.faa",
        rows_temp = output_folder_uniprot_hmm + "/temp_rows.tsv"
    params:
        max_length = solosaved_max_length,
        e_value = solosaved_e_value,
        temp_hmm = output_folder_uniprot_hmm + "/temp_hmm.tsv",
    threads: thread_hogger
    conda:
        "envs/hmmer.yaml"
    log:
        out = output_folder_uniprot_hmm + "/logs/solosaved_hmmscan.out",
        err = output_folder_uniprot_hmm + "/logs/solosaved_hmmscan.err"
    shell:
        '''
        hmmscan --domtblout {params.temp_hmm} --cpu {threads} -E {params.e_value} {input.db} {input.query} &> /dev/null
        echo "Removing commented rows" >> {log.out}
        grep -v "^#" {params.temp_hmm} > {output.rows_temp} ||:
        echo "Writing header" >> {log.out}
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        echo "Checking if hits were found" >> {log.out}
        if [ -s {output.rows_temp} ]; then 
            echo "Hits found" >> {log.out}
            cat {output.rows_temp} >> {output.hmm_rows}
        else
            echo "No hits found" >> {log.out}
            touch {output.hmm_rows}
        fi
        '''


### PLASMIDS ###


# rule annotate_plasmids:
#     '''
#     Uses Prokka to annotate input plasmids
#     '''
#     input:
#         plasmid_fasta = input_plasmid_nt_split_folder + "/{plasmid}.fna"
#     output:
#         proteins = output_folder_plasmids + "/{plasmid}/{plasmid}.faa",
#         nt = output_folder_plasmids + "/{plasmid}/{plasmid}.fna",
#         gff = output_folder_plasmids + "/{plasmid}/{plasmid}.gff",
#         gbk = output_folder_plasmids + "/{plasmid}/{plasmid}.gbk"
#     conda: "envs/prokka.yaml"
#     params:
#         outfolder = output_folder_plasmids + "/{plasmid}",
#         rn_hmm = ring_nuclease_hmm_db
#     threads: 5
#     log:
#         err = output_folder_plasmids + "/{plasmid}/logs/prokka.err",
#         out = output_folder_plasmids + "/{plasmid}/logs/prokka.out"
#     shell:
#         '''
#         prokka --outdir {params.outfolder} {input.plasmid_fasta} --prefix {wildcards.plasmid} --cpus {threads} --force 2> {log.err} 1> {log.out}
#         touch {output}
#         '''

# rule solosaved_hmmscan_plasmids:
#     '''
#     Runs hmmscan using saved HMMs against annotated plasmids
#     '''
#     input:
#         proteins = plsdb_proteins
#     output:
#         temp_rows = output_folder_plasmids_hmm + "/{plasmid}/{plasmid}_RN_temp.tsv",
#         hmm_rows = output_folder_plasmids_hmm + "/{plasmid}/{plasmid}_RN_hmm.tsv"
#     threads: thread_med
#     params:
#         saved_db = carfsaved_hmm_db_path,
#         evalue = "1e-8",
#         temp_hmm = output_folder_plasmids_hmm + "/{plasmid}/{plasmid}_RN_hmm.temp"
#     conda: "envs/hmmer.yaml"
#     threads: thread_min
#     log:
#         out = output_folder_plasmids_hmm + "/{plasmid}/logs/RN_search.out",
#         err = output_folder_plasmids_hmm + "/{plasmid}/logs/RN_search.err"
#     shell:
#         '''
#         hmmscan --domtblout {params.temp_hmm} --cpu {threads} -E {params.evalue} {params.saved_db} {input.proteins} &> /dev/null
#         echo "Removing commented rows" >> {log.out}
#         grep -v "^#" {params.temp_hmm} > {output.temp_rows} ||:
#         echo "Writing header" >> {log.out}
#         echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
#         echo "Checking if hits were found" >> {log.out}
#         if [ -s {output.temp_rows} ]; then 
#             echo "Hits found for {wildcards.plasmid}" >> {log.out}
#             cat {output.temp_rows} >> {output.hmm_rows}
#         else
#             echo "No hits found for {wildcards.plasmid}" >> {log.out}
#             touch {output.hmm_rows}
#         fi
#         '''

rule solosaved_hmmscan_plasmids:
    '''
    Runs hmmscan using saved HMMs against annotated plasmids
    '''
    input:
        proteins = plsdb_proteins
    output:
        temp_rows = output_folder_plasmids_hmm + "/plasmid_solosaved_temp.tsv",
        hmm_rows = output_folder_plasmids_hmm + "/plasmid_solosaved_temp_hmm.tsv"
    threads: thread_med
    params:
        saved_db = carfsaved_hmm_db_path,
        evalue = "1e-8",
        temp_hmm = output_folder_plasmids_hmm + "/plasmid_solosaved_temp_hmm.temp"
    conda: "envs/hmmer.yaml"
    threads: thread_hogger
    log:
        out = output_folder_plasmids_hmm + "/logs/plasmid_solosaved_hmm.out",
        err = output_folder_plasmids_hmm + "/logs/plasmid_solosaved_hmm.err"
    shell:
        '''
        hmmscan --domtblout {params.temp_hmm} --cpu {threads} -E {params.evalue} {params.saved_db} {input.proteins} &> /dev/null
        echo "Removing commented rows" >> {log.out}
        grep -v "^#" {params.temp_hmm} > {output.temp_rows} ||:
        echo "Writing header" >> {log.out}
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        echo "Checking if hits were found" >> {log.out}
        if [ -s {output.temp_rows} ]; then 
            echo "Hits found for {wildcards.plasmid}" >> {log.out}
            cat {output.temp_rows} >> {output.hmm_rows}
        else
            echo "No hits found for {wildcards.plasmid}" >> {log.out}
            touch {output.hmm_rows}
        fi
        '''

# rule plasmid_locus_visualiser:
#     '''
#     This rule visualises a given range of a gff file.
#     Using plasmid as wildcard.
#     GFF comes from Prokka output.
#     Coordinates for RN come from HMMER output.
#     '''
#     input:
#         gff = rules.annotate_plasmids.output.gff,
#         hmm_hits = rules.solosaved_hmmscan_plasmids.output.hmm_rows,
#         gbk = rules.annotate_plasmids.output.gbk
#     output:
#         done = output_folder_plasmid_locus_viz + "/{plasmid}/done"
#     conda:
#         "envs/locus_visualiser.yaml"
#     params:
#         outdir = output_folder_plasmid_locus_viz + "/{plasmid}",
#         gbk_out = output_folder_plasmid_locus_viz + "/04_locus_viz_gbk"
#     log:
#         out = output_folder_plasmid_locus_viz + "/logs/{plasmid}.out",
#         err = output_folder_plasmid_locus_viz + "/logs/{plasmid}.err",
#     threads: thread_singular
#     shell:
#         '''
#         if [ ! -d {params.gbk_out} ]; then
#             mkdir {params.gbk_out}
#         fi
#         python scripts/locus_viz_gff_hmm.py --sample {wildcards.plasmid} --output_folder {params.outdir} --prokka_gff {input.gff} --hmm_hits {input.hmm_hits} --gbk {input.gbk} --gbk_out_folder {params.gbk_out} 2> {log.err} 1> {log.out}
#         touch {output.done}
#         '''

### PHAGES ###


def aggregate_millard_phage_solosaved_analysis(wildcards):
    '''
    This function is used to aggregate the outputs of the millard_phage_RN_analysis rule.
    If not working, create a checkpoint to generate wildcards first.
    '''
    checkpoint_output = checkpoints.divide_millard_phages_to_folders.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{phage}/done.txt")).phage
    return expand(output_folder_phage_hmm + "/{phage}/{phage}_RN_hmm.tsv", phage=cvals)


checkpoint divide_millard_phages_to_folders:
    '''
    Takes in the genome Millard genomes.fa file and creates a folder and fasta for each phage in the fasta.
    The headers are in the format >phage_name description. We only want to keep the phage_name for the folder.
    '''
    output: directory(output_folder_phage_genomes)
    input:
        millard_fa = millard_fa, #nt multifasta of all phage genomes
        millard_proteins = millard_proteins #aa multifasta of all phage genones. Headers are >phage_name_runningnumber
    params:
        outdir = output_folder_phage_genomes,
    threads: thread_small
    run:
        from Bio import SeqIO
        import os
        print("Starting divide_millard_phages_to_folders")
        
        # Create a dictionary to store sequences for each phage
        phage_sequences = {}
        phage_protein_sequences = {}
        
        # Parse the millard_fa file and store sequences in the dictionary
        for record in SeqIO.parse(input.millard_fa, "fasta"):
            phage_name = record.id.split(" ")[0]
            if phage_name not in phage_sequences:
                phage_sequences[phage_name] = []
            phage_sequences[phage_name].append(record)
        
        # Parse the millard_proteins file and store sequences in the dictionary
        for record in SeqIO.parse(input.millard_proteins, "fasta"):
            phage_name = record.id.split("_")[0]
            if phage_name not in phage_protein_sequences:
                phage_protein_sequences[phage_name] = []
            phage_protein_sequences[phage_name].append(record)
        
        # Write sequences to files for each phage
        for phage_name, sequences in phage_sequences.items():
            phage_dir = f"{params.outdir}/{phage_name}"
            os.makedirs(phage_dir, exist_ok=True)
            
            # Write nucleotide sequences to file
            with open(f"{phage_dir}/{phage_name}.fna", "w") as out:
                SeqIO.write([seq for seq in sequences if seq.id.split(" ")[0] == phage_name], out, "fasta")
            
        for phage_name, sequences in phage_protein_sequences.items():
            phage_dir = f"{params.outdir}/{phage_name}"
            os.makedirs(phage_dir, exist_ok=True)

            # Write protein sequences to file
            with open(f"{phage_dir}/{phage_name}_proteins.faa", "w") as out:
                SeqIO.write([seq for seq in sequences if seq.id.split("_")[0] == phage_name], out, "fasta")
            
            # Write a "done" file in every phage folder
            with open(f"{phage_dir}/done.txt", "w") as out:
                out.write("done")

rule millard_phage_solosaved_hmm:
    '''
    Searches for ring nucleases in the Millard phage proteomes using HMMER and
    our pre-made RN hmm profiles similar to rule ring_nucleases
    '''
    input:
        phage_proteins = output_folder_phage_genomes + "/{phage}/{phage}_proteins.faa"
    output:
        temp_rows = output_folder_phage_hmm + "/{phage}/{phage}_RN_temp.tsv",
        hmm_rows = output_folder_phage_hmm + "/{phage}/{phage}_RN_hmm.tsv"
    params:
        solosaved_db = carfsaved_hmm_db_path,
        evalue = "1e-8",
        temp_hmm = output_folder_phage_hmm + "/{phage}/{phage}_RN_hmm.temp"
    conda: "envs/hmmer.yaml"
    threads: thread_small
    log:
        out = output_folder_phage_hmm + "/{phage}/logs/RN_search.out",
        err = output_folder_phage_hmm + "/{phage}/logs/RN_search.err"
    shell:
        '''
        hmmscan --domtblout {params.temp_hmm} --cpu {threads} -E {params.evalue} {params.solosaved_db} {input.phage_proteins} &> /dev/null
        echo "Removing commented rows" >> {log.out}
        grep -v "^#" {params.temp_hmm} > {output.temp_rows} ||:
        echo "Writing header" >> {log.out}
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        echo "Checking if hits were found" >> {log.out}
        if [ -s {output.temp_rows} ]; then 
            echo "Hits found for {wildcards.phage}" >> {log.out}
            cat {output.temp_rows} >> {output.hmm_rows}
        else
            echo "No hits found for {wildcards.phage}" >> {log.out}
            touch {output.hmm_rows}
        fi
        '''

rule concatenate_millard_phage_solosaved_analysis:
    '''
    Concatenates the results of the millard_phage_solosaved_analysis rule.
    '''
    input: aggregate_millard_phage_solosaved_analysis
    output: output_folder_phage_hmm + "/millard_phage_solosaved_analysis.tsv"
    shell:
        """
        cat {input} > {output}
        """