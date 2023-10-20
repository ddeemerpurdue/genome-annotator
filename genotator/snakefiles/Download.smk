import os
from datetime import datetime
tstamp = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
configfile: "../config.yaml"

rule dl_all:
    input:
        f'{config["verified_contigs"]}/complete.tkn'

# INTERNAL PHASE
rule Download_Taxa:
    params:
        taxon_id = f'{config["taxon_id"]}',
        source = f'{config["download_source"]}',
        filename = f'{config["download_filename"]}.zip'
    output:
        z_file = f'{config["download_filename"]}.zip'
    log: f'{config["log_folder"]}/Download_Taxa__{config["taxon_id"]}.{config["log_id"]}.log'
    group: "main"
    shell:
        """
        datasets download genome --assembly-source {params.source} taxon {params.taxon_id} --filename {params.filename} &> {log}
        """

checkpoint Process_NCBI_Zip:
    input:
        z_file = f'{config["download_filename"]}.zip'
    output:
        directory(f'{config["fastas"]}')
    log:
        f'{config["log_folder"]}/Process_NCBI__{config["taxon_id"]}.{config["log_id"]}.log'
    shell:
        """
        python ../scripts/parse_ncbi_data.py -i {input.z_file} -d {output} &> {log}
        """


rule Clean_Fastas:
    input:
        fasta = f'{config["fastas"]}/{{sample}}.fna'
    output:
        clean_fasta = f'{config["verified_contigs"]}/{{sample}}.fasta'
    log:
        f'{config["log_folder"]}/Clean_Fastas__-{{sample}}{config["taxon_id"]}.{config["log_id"]}.log'
    shell:
        """
        anvi-script-reformat-fasta {input.fasta} -o {output.clean_fasta} --simplify-names --seq-type NT &> {log}
        """


def grab_ncbi_fastas(wildcards):
    checkpoint_output = checkpoints.Process_NCBI_Zip.get(**wildcards).output[0]
    file_names = expand(
        f'{config["verified_contigs"]}/{{sample}}.{config["extension"]}',
        sample = glob_wildcards(os.path.join(checkpoint_output, "{sample}.fna")).sample)
    return file_names


rule Finished:
    input:
        grab_ncbi_fastas
    output:
        f'{config["verified_contigs"]}/complete.tkn'
    shell:
        """
        touch {output}
        """
