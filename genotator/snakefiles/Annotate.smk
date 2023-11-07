import os
from datetime import datetime
tstamp = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
# configfile: "../config.yaml"

# include: "Filter.smk"
# module filter_workflow:
#     snakefile: "Filter.smk"

# use rule * from filter_workflow as other_*
# use rule Place_Final_Genomes from filter_workflow as other_Place_Final_Genomes with:
#     output: directory(f'{config["verified_contigs"]}')


samples, = glob_wildcards(
    f'{config["filtered_contigs"]}/{{samples}}.{config["extension"]}')

rule all:
    input:
        expand(f'{config["gff_final"]}/{{sample}}.gff3', sample=samples)

rule Create_Database:
    input:
        gen = f'{config["filtered_contigs"]}/{{sample}}.{config["extension"]}'
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
        touch(f'{config["annotations"]}/HMMs/{{sample}}.run_hmm')
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
        estimate_out = f'{config["annotations"]}/Other/{{sample}}.estimate_scgs',
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
        touch(f'{config["annotations"]}/Cog/{{sample}}.cog')
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
        touch(f'{config["annotations"]}/Kegg/{{sample}}.kegg')
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
        touch(f'{config["annotations"]}/Pfam/{{sample}}.pfam')
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
        f'{config["faas"]}/{{sample}}.faa'
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
        f'{config["gene_calls"]}/{{sample}}.gff'
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
        f'{config["faas"]}/{{sample}}.faa'
    output:
        f'{config["annotations"]}/RAST/{{sample}}-RAST-FAA.txt'
    log: f'Logs/RAST_Run__{{sample}}.{config["log_id"]}.log'
    group: "main"
    priority: 10
    shell:
        """
        svr_assign_using_figfams < {input} > {output}
        """

rule RAST_Reformat:
    input:
        f'{config["annotations"]}/RAST/{{sample}}-RAST-FAA.txt'
    output:
        f'{config["annotations"]}/RAST/{{sample}}-RAST-FAA.txt.anvio'
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
        imp = f'{config["annotations"]}/RAST/{{sample}}-RAST-FAA.txt.anvio'
    output:
        touch(f'{config["annotations"]}/RAST/{{sample}}.rast_added')
    log: f'Logs/RAST_Import__{{sample}}.{config["log_id"]}.log'
    group: "main"
    priority: 10
    shell:
        """
        anvi-import-functions -c {input.db} -i {input.imp} &> {log}
        """

# rule CAZyme_Run:
#     input:
#         fasta = f'{config["filtered_contigs"]}/{{sample}}-VERIFIED.{config["extension"]}'
#     params:
#         outdir = f'{config["annotations"]}/DBCan4/{{sample}}/',
#         outpref = f'{{sample}}_',
#         db = f'{config["cazyme_db"]}'
#     output:
#         final = f'{config["annotations"]}/DBCan4/{{sample}}/{{sample}}_dbsub.out'
#     log: f'Logs/CAZyme_Run__{{sample}}.{config["log_id"]}.log'
#     threads: config["threads"]
#     shell:
#         """
#         run_dbcan {input.fasta} prok --out_dir {params.outdir} --tf_cpu {threads} --out_pre {params.outpref} \
#         --db_dir {params.db} -c cluster
#         """

rule CAZyme_Run:
    input:
        faa = f'{config["faas"]}/{{sample}}.faa',
        gene_calls = f'{config["gene_calls"]}/{{sample}}.gff'
    params:
        outdir = f'{config["annotations"]}/DBCan4/{{sample}}/',
        outpref = f'{{sample}}_',
        db = f'{config["cazyme_db"]}'
    output:
        final = f'{config["annotations"]}/DBCan4/{{sample}}/{{sample}}_dbsub.out'
    log: f'Logs/CAZyme_Run__{{sample}}.{config["log_id"]}.log'
    threads: config["threads"]
    shell:
        """
        run_dbcan {input.faa} meta --out_dir {params.outdir} --tf_cpu {threads} --out_pre {params.outpref} \
        --db_dir {params.db} -c {input.gene_calls} &> {log}
        """

rule CAZyme_Reformat:
    input:
        f'{config["annotations"]}/DBCan4/{{sample}}/{{sample}}_dbsub.out'
    params:
        mapper = f'{config["cazyme_mapper"]}'
    output:
        f'{config["annotations"]}/DBCan4/{{sample}}/results.anvio.tbl'
    log: f'Logs/CAZyme_Reformat__{{sample}}.{config["log_id"]}.log'
    shell:
        """
        python scripts/cazyme-table.py {input} {params.mapper} {output} &> {log}
        """

rule CAZyme_Import:
    input:
        db = ancient(f'{config["contig_db"]}/{{sample}}.db'),
        imp = f'{config["annotations"]}/DBCan4/{{sample}}/results.anvio.tbl'
    output:
        touch(
            f'{config["annotations"]}/DBCan4/{{sample}}/{{sample}}.cazy_added')
    log: f'Logs/CAZyme_Import__{{sample}}.{config["log_id"]}.log'
    group: "main"
    priority: 10
    shell:
        """
        anvi-import-functions -c {input.db} -i {input.imp}  &> {log}
        """

rule TigrFam_Run:
    input:
        f'{config["faas"]}/{{sample}}.faa'
    params:
        db = f'{config["tigrfam_db"]}'
    output:
        one = f'{config["annotations"]}/TigrFamResults/{{sample}}.hmmer.TIGR.hmm',
        two = f'{config["annotations"]}/TigrFamResults/{{sample}}.hmmer.TIGR.tbl'
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
        f'{config["annotations"]}/TigrFamResults/{{sample}}.hmmer.TIGR.tbl'
    params:
        tfam = f'{config["tigrfam_roles"]}'
    output:
        f'{config["annotations"]}/TigrFamResults/{{sample}}.hmmer.TIGR.anvio.tbl'
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
        imp = f'{config["annotations"]}/TigrFamResults/{{sample}}.hmmer.TIGR.anvio.tbl'
    output:
        touch(f'{config["annotations"]}/TigrFamResults/{{sample}}.tigr_added')
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
        hmms = f'{config["annotations"]}/HMMs/{{sample}}.run_hmm',
        cogs = f'{config["annotations"]}/Cog/{{sample}}.cog',
        kegg = f'{config["annotations"]}/Kegg/{{sample}}.kegg',
        pfams = f'{config["annotations"]}/Pfam/{{sample}}.pfam',
        # figfams = f'{config["annotations"]}/RAST/{{sample}}.rast_added',
        tigrfams = f'{config["annotations"]}/TigrFamResults/{{sample}}.tigr_added',
        #cazys = f'{config["annotations"]}/DBCan4/{{sample}}/{{sample}}.cazy_added'
    output:
        out = f'{config["annotations"]}/Annotations-Exported/{{sample}}.gff3'
    log: f'Logs/Export_Annotations__{{sample}}.{config["log_id"]}.log'
    params:
        # annotations = "FigFams,KEGG_Module,COG20_PATHWAY,TIGRFAM,CAZyme,KOfam,KEGG_Class,COG20_FUNCTION,Pfam,COG20_CATEGORY"
        annotations = "KEGG_Module,COG20_PATHWAY,TIGRFAM,KOfam,KEGG_Class,COG20_FUNCTION,Pfam,COG20_CATEGORY"
    group: "exporting"
    priority: 0
    shell:
        """
        anvi-export-functions -c {input.db} --annotation-sources {params.annotations} -o {output.out} &> {log}
        """

rule Reformat_Gff3:
    input:
        anvio_annotations = f'{config["annotations"]}/Annotations-Exported/{{sample}}.gff3',
        gene_calls = f'{config["gene_calls"]}/{{sample}}.gff',
    output:
        gff3_final = f'{config["gff_final"]}/{{sample}}.gff3'
    log: f'Logs/Reformat_Gff3__{{sample}}.{config["log_id"]}.log'
    group: "exporting"
    priority: 0
    shell:
        """
        scripts/combineFunctionsAndGeneCalls.py {input.gene_calls} {input.anvio_annotations} {output.gff3_final} &> {log}
        """
