import os
from datetime import datetime
tstamp = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
print(f'In Download.smk, wd == {os.getcwd()}')
configfile: "../genotator/config.yaml"

rule dl_all:
    input:
        f'{config["verified_contigs"]}/complete.tkn',
        f'{config["checkm"]}/complete.tkn'
        # f'{config["verified_contigs"]}/{{sample}}.{config["extension"]}'

# INTERNAL PHASE
rule Download_Taxonomy:
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
        python ../genotator/scripts/parse_ncbi_data.py -i {input.z_file} -d {output} &> {log}
        """


def grab_ncbi_fastas(wildcards):
    checkpoint_output = checkpoints.Process_NCBI_Zip.get(**wildcards).output[0]
    verified_fastas = expand(
        f'{config["verified_contigs"]}/{{sample}}.{config["extension"]}',
        sample = glob_wildcards(os.path.join(checkpoint_output, "{sample}.fna")).sample)
    return verified_fastas


def grab_checkm_fastas(wildcards):
    checkpoint_output = checkpoints.Process_NCBI_Zip.get(**wildcards).output[0]
    checkm_fasta = expand(
        f'{config["checkm"]}/{{sample}}/lineage.ms',
        sample = glob_wildcards(os.path.join(checkpoint_output, "{sample}.fna")).sample)
    return checkm_fasta


rule Verify_Fastas:
    input:
        fasta = f'{config["fastas"]}/{{sample}}.fna'
    output:
        clean_fasta = f'{config["verified_contigs"]}/{{sample}}.{config["extension"]}'
    log:
        f'{config["log_folder"]}/Clean_Fastas__{{sample}}{config["taxon_id"]}.{config["log_id"]}.log'
    shell:
        """
        anvi-script-reformat-fasta {input.fasta} -o {output.clean_fasta} --simplify-names --seq-type NT &> {log}
        """


rule Run_CheckM:
    input:
        fasta = f'{config["fastas"]}/{{sample}}.fna'
    params:
        threads=config['threads'],
        out=f'{config["checkm"]}/{{sample}}/'
    output:
        out=f'{config["checkm"]}/{{sample}}/lineage.ms'
    log:
        "logs/qc/{sample}.log"
    shell:
        """
        checkm lineage_wf -t {params.threads} -x fa {input.fasta} {params.out} &> {log}
        """


rule Verify_All:
    input:
        grab_ncbi_fastas
    output:
        f'{config["verified_contigs"]}/complete.tkn'
    shell:
        """
        touch {output}
        """

rule Checkm_All:
    input:
        grab_checkm_fastas
    output:
        f'{config["checkm"]}/complete.tkn'
    shell:
        """
        touch {output}
        """




