import os
from datetime import datetime
tstamp = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
configfile: "config.yaml"

assemblies, = glob_wildcards(
    f'{config["contigs"]}/{{samples}}.{config["extension"]}')

rule all:
    input:
        # expand("GFF3-Final/{sample}.gff3", sample=assemblies),
        #expand("Annotations/Other/{dbs}.estimate_scgs", dbs=glob_wildcards(f'{config["contig_db"]}/{{value}}.db').value),

# INTERNAL PHASE
rule Download_Taxa:
    params:
        taxon_id = f'{config["taxon_id"]}',
        source = f'{config["download_source"]}'
        filename = f'{config["download_filename"]}.zip'
    output:
        final = f'{config["verified_contigs"]}/{{sample}}-VERIFIED.{config["extension"]}'
        final = f'{{}}'
    log: f'Logs/Reformat_Fasta__{{sample}}.{config["log_id"]}.log'
    group: "main"
    shell:
        """
        datasets download genome --assembly-source {params.source} taxon {params.taxon} --filename {params.filename}
        anvi-script-reformat-fasta {input.origin} -o {output.final} --simplify-names --seq-type NT &> {log}
        """

'''
unzip {output.final}
cd ncbi_dataset/data/

'''