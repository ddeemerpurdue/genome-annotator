import os
from datetime import datetime
tstamp = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
configfile: "config.yaml"

assemblies, = glob_wildcards(
    f'{config["contigs"]}/{{samples}}.{config["extension"]}')

rule all:
    input:
        expand("GFF3-Final/{sample}.gff3", sample=assemblies),
        #expand("Annotations/Other/{dbs}.estimate_scgs", dbs=glob_wildcards(f'{config["contig_db"]}/{{value}}.db').value),

# INTERNAL PHASE
rule Reformat_Fasta:
    input:
        origin = f'{config["contigs"]}/{{sample}}.{config["extension"]}'
    output:
        final = f'{config["verified_contigs"]}/{{sample}}-VERIFIED.{config["extension"]}'
    log: f'Logs/Reformat_Fasta__{{sample}}.{config["log_id"]}.log'
    group: "main"
    shell:
        """
        anvi-script-reformat-fasta {input.origin} -o {output.final} --simplify-names --seq-type NT &> {log}
        """

rule Create_Database:
    input:
        gen = f'{config["verified_contigs"]}/{{sample}}-VERIFIED.{config["extension"]}'
    output:
        out = f'{config["contig_db"]}/{{sample}}.db'
    params:
        title = f'{config["database_name"]}',
    log: f'Logs/Create_Database__{{sample}}.{config["log_id"]}.log'
    threads: config["threads"]
    group: "main"
    shell:
        """
        anvi-gen-contigs-database -T {threads} -f {input.gen} -n '{wildcards.sample} {params.title}' -o {output.out} &> {log}
        """

rule Run_HMMs:
    input:
        one = ancient(f'{config["contig_db"]}/{{sample}}.db'),
    output:
        touch("Annotations/HMMs/{sample}.run_hmm")
    log: f'Logs/Run_HMMs__{{sample}}.{config["log_id"]}.log'
    threads: config["threads"]
    group: "main"
    priority: 25
    shell:
        """
        anvi-run-hmms -c {input.one} -T {threads} --just-do-it --also-scan-trnas &> {log}
        """

rule Run_SCG_Taxonomy:
    input:
        db = ancient(f'{config["contig_db"]}/{{sample}}.db')
    output:
        estimate_out = "Annotations/Other/{sample}.estimate_scgs",
    log:
        f'Logs/Run_SCGs__{{sample}}.{config["log_id"]}.log'
    priority: 90
    threads: config["threads"]
    shell:
        """
        anvi-run-scg-taxonomy -c {input.db} -T {threads} &> {log}
        anvi-estimate-scg-taxonomy -c {input.db} --output-file {output.estimate_out} -T {threads} --metagenome-mode
        """

rule Run_COGs:
    input:
        ancient(f'{config["contig_db"]}/{{sample}}.db')
    output:
        touch("Annotations/Cog/{sample}.cog")
    log: f'Logs/Run_COGs__{{sample}}.{config["log_id"]}.log'
    threads: config["threads"]
    group: "main"
    priority: 25
    shell:
        """
        anvi-run-ncbi-cogs -T {threads} -c {input} &> {log}
        """

rule Run_Kegg:
    input:
        ancient(f'{config["contig_db"]}/{{sample}}.db')
    output:
        touch("Annotations/Kegg/{sample}.kegg")
    log: f'Logs/Run_Kegg__{{sample}}.{config["log_id"]}.log'
    threads: config["threads"]
    group: "main"
    priority: 25
    shell:
        """
        anvi-run-kegg-kofams -c {input} -T {threads} --just-do-it &> {log}
        """

rule Run_PFams:
    input:
        ancient(f'{config["contig_db"]}/{{sample}}.db')
    output:
        touch("Annotations/Pfam/{sample}.pfam")
    log: f'Logs/Run_PFams__{{sample}}.{config["log_id"]}.log'
    threads: config["threads"]
    group: "main"
    priority: 25
    shell:
        """
        anvi-run-pfams -c {input} -T {threads} &> {log}
        """

# EXPORTING PHASE
rule Export_FAA:
    input:
        ancient(f'{config["contig_db"]}/{{sample}}.db')
    output:
        "FAAs/{sample}.faa"
    log: f'Logs/Export_FAA__{{sample}}.{config["log_id"]}.log'
    group: "main"
    priority: 15
    shell:
        """
        anvi-get-sequences-for-gene-calls -c {input} --get-aa-sequences -o {output} &> {log}
        """

rule Export_Gene_Calls:
    input:
        ancient(f'{config["contig_db"]}/{{sample}}.db')
    output:
        "GeneCalls/{sample}.gff"
    log: f'Logs/Export_Gene_Calls__{{sample}}.{config["log_id"]}.log'
    group: "main"
    priority: 15
    shell:
        """
        anvi-export-gene-calls -c {input} -o {output} --gene-caller prodigal &> {log}
        """

# EXTERNAL ANNOTATION
rule RAST_Run:
    input:
        "FAAs/{sample}.faa"
    output:
        "Annotations/RAST/{sample}-RAST-FAA.txt"
    log: f'Logs/RAST_Run__{{sample}}.{config["log_id"]}.log'
    group: "main"
    priority: 10
    shell:
        """
        svr_assign_using_figfams < {input} > {output}
        """

rule RAST_Reformat:
    input:
        "Annotations/RAST/{sample}-RAST-FAA.txt"
    output:
        "Annotations/RAST/{sample}-RAST-FAA.txt.anvio"
    log: f'Logs/RAST_Reformat__{{sample}}.{config["log_id"]}.log'
    group: "main"
    priority: 10
    shell:
        """
        rast-table.py {input} {output} &> {log}
        """

rule RAST_Import:
    input:
        db = ancient(f'{config["contig_db"]}/{{sample}}.db'),
        imp = "Annotations/RAST/{sample}-RAST-FAA.txt.anvio"
    output:
        touch("Annotations/RAST/{sample}.rast_added")
    log: f'Logs/RAST_Import__{{sample}}.{config["log_id"]}.log'
    group: "main"
    priority: 10
    shell:
        """
        anvi-import-functions -c {input.db} -i {input.imp} &> {log}
        """

rule CAZyme_Run:
    input:
        fasta = f'{config["verified_contigs"]}/{{sample}}-VERIFIED.{config["extension"]}'
    params:
        outdir = f'Annotations/DBCan4/{{sample}}/',
        outpref = f'{{sample}}_',
        db = f'{config["cazyme_db"]}'
    output:
        final = f'Annotations/DBCan4/{{sample}}/{{sample}}_dbsub.out'
    log: f'Logs/CAZyme_Run__{{sample}}.{config["log_id"]}.log'
    threads: config["threads"]
    shell:
        """
        run_dbcan {input.fasta} prok --out_dir {params.outdir} --tf_cpu {threads} --out_pre {params.outpref} \
        --db_dir {params.db} -c cluster
        """

rule CAZyme_Reformat:
    input:
        f'Annotations/DBCan4/{{sample}}/{{sample}}_dbsub.out'
    params:
        mapper = f'{config["cazyme_mapper"]}'
    output:
        f'Annotations/DBCan4/{{sample}}/results.anvio.tbl'
    log: f'Logs/CAZyme_Reformat__{{sample}}.{config["log_id"]}.log'
    shell:
        """
        python scripts/cazyme-table.py {input} {params.mapper} {output} &> {log}
        """

rule CAZyme_Import:
    input:
        db = ancient(f'{config["contig_db"]}/{{sample}}.db'),
        imp = f'Annotations/DBCan4/{{sampe}}/results.anvio.tbl'
    output:
        touch("Annotations/DBCan4/{sample}/{sample}.tigr_added")
    log: f'Logs/CAZyme_Import__{{sample}}.{config["log_id"]}.log'
    group: "main"
    priority: 10
    shell:
        """
        anvi-import-functions -c {input.db} -i {input.imp}  &> {log}
        """

rule TigrFam_Run:
    input:
        "FAAs/{sample}.faa"
    params:
        db = f'{config["tigrfam_db"]}'
    output:
        one = "Annotations/TigrFamResults/{sample}.hmmer.TIGR.hmm",
        two = "Annotations/TigrFamResults/{sample}.hmmer.TIGR.tbl"
    log: f'Logs/TigrFam_Run__{{sample}}.{config["log_id"]}.log'
    threads: config["threads"]
    group: "main"
    priority: 10
    shell:
        """
        hmmsearch -o {output.one} --tblout {output.two} --cpu {threads} {params.db} {input} &> {log}
        """

rule TigrFam_Reformat:
    input:
        "Annotations/TigrFamResults/{sample}.hmmer.TIGR.tbl"
    params:
        tfam = f'{config["tigrfam_roles"]}'
    output:
        "Annotations/TigrFamResults/{sample}.hmmer.TIGR.anvio.tbl"
    log: f'Logs/TigrFam_Reformat__{{sample}}.{config["log_id"]}.log'
    group: "main"
    priority: 10
    shell:
        """
        scripts/tigrfam-table.py {input} {params.tfam} {output} &> {log}
        """

rule TigrFam_Import:
    input:
        db = ancient(f'{config["contig_db"]}/{{sample}}.db'),
        imp = "Annotations/TigrFamResults/{sample}.hmmer.TIGR.anvio.tbl"
    output:
        touch("Annotations/TigrFamResults/{sample}.tigr_added")
    log: f'Logs/TigrFam_Import__{{sample}}.{config["log_id"]}.log'
    group: "main"
    priority: 10
    shell:
        """
        anvi-import-functions -c {input.db} -i {input.imp}  &> {log}
        """

# EXPORTING ANNOTATIONS
rule Export_Annotations:
    input:
        db = ancient(f'{config["contig_db"]}/{{sample}}.db'),
        hmms = "Annotations/HMMs/{sample}.run_hmm",
        cogs = "Annotations/Cog/{sample}.cog",
        kegg = "Annotations/Kegg/{sample}.kegg",
        pfams = "Annotations/Pfam/{sample}.pfam",
        figfams = "Annotations/RAST/{sample}.rast_added",
        tigrfams = "Annotations/TigrFamResults/{sample}.tigr_added",
    output:
        out = "Annotations/Annotations-Exported/{sample}.gff3"
    log: f'Logs/Export_Annotations__{{sample}}.{config["log_id"]}.log'
    params:
        annotations = "FigFams,KEGG_Module,COG20_PATHWAY,TIGRFAM,KOfam,KEGG_Class,COG20_FUNCTION,Pfam,COG20_CATEGORY"
    group: "exporting"
    priority: 0
    shell:
        """
        anvi-export-functions -c {input.db} --annotation-sources {params.annotations} -o {output.out} &> {log}
        """

rule Reformat_Gff3:
    input:
        anvio_annotations = "Annotations/Annotations-Exported/{sample}.gff3",
        gene_calls = "GeneCalls/{sample}.gff",
    output:
        gff3_final = "GFF3-Final/{sample}.gff3"
    log: f'Logs/Reformat_Gff3__{{sample}}.{config["log_id"]}.log'
    group: "exporting"
    priority: 0
    shell:
        """
        scripts/combineFunctionsAndGeneCalls.py {input.gene_calls} {input.anvio_annotations} {output.gff3_final} &> {log}
        """
