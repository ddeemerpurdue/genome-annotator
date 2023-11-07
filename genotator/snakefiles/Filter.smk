import os
import shutil
print(f'Current working directory: {os.getcwd()}')
from datetime import datetime
tstamp = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")

include: "./Download.smk"

# ~/Documents/bioinf-programs/FastANI/build/
# For FastANI

rule filter_all:
    input:
        a=f'{config["verified_contigs"]}/complete.tkn',
        b=f'{config["checkm"]}/complete.tkn',
        c=f'{config["fastani"]}/{config["taxon_id"]}-Filtered.txt',
        d=f'{config["final_contigs"]}/complete.tkn'
        # b="CheckM/{sample}/lineage.ms"


rule Create_Qlist:
    input:
        f'{config["verified_contigs"]}/complete.tkn'
    params:
        extension = f'{config["extension"]}',
        direct = f'{config["verified_contigs"]}'
    output:
        f'{config["fastani"]}/{config["taxon_id"]}.txt'
    log: f'{config["log_folder"]}/Create_Qlist__{config["taxon_id"]}.{config["log_id"]}.log'
    group: "main"
    shell:
        """
        python ../scripts/create_qlist.py -d {params.direct} -e {params.extension} -o {output} &> {log}
        """


rule Run_FastAni:
    input: f'{config["fastani"]}/{config["taxon_id"]}.txt'
    params:
        threads = f'{config["threads"]}'
    output: f'{config["fastani"]}/{config["taxon_id"]}-FastANI.txt'
    log: f'{config["log_folder"]}/Run_FastAni__{config["taxon_id"]}.{config["log_id"]}.log'
    shell:
        """
        fastANI -t {params.threads} --rl {input} --ql {input} -o {output} &> {log}
        """


rule Filter_Ani:
    input:
        f'{config["fastani"]}/{config["taxon_id"]}-FastANI.txt'
    params:
        threshold = f'{config["fastani_thresh"]}'
    output:
        ani = f'{config["fastani"]}/{config["taxon_id"]}-Filtered.txt',
        genome_list = f'{config["fastani"]}/{config["taxon_id"]}-Genomes.txt',
    log:
        f'{config["log_folder"]}/Filter_Ani__{config["taxon_id"]}.{config["log_id"]}.log'
    shell:
        """
        python ../scripts/parseFastAniNonRedundant.py -i {input} -t {params.threshold} -o {output.ani} &> {log}
        cut -f1 {output.ani} > {output.genome_list}
        """


checkpoint Place_Genomes:
    input: f'{config["fastani"]}/{config["taxon_id"]}-Genomes.txt'
    params:
        gen_dir = directory(f'{config["verified_contigs"]}')
    output: directory(f'{config["filtered_contigs"]}')
    run:
        os.makedirs(output[0], exist_ok=True)
        with open(input[0], 'r') as infile:
            for line in infile.readlines():
                genome = line.strip()
                src = os.path.join(params.gen_dir, genome)
                dest = os.path.join(output[0], genome)
                shutil.copyfile(src, dest)
    

def grab_genomes(wildcards):
    checkpoint_output = checkpoints.Place_Genomes.get(**wildcards).output[0]
    file_names = expand(
        f'{config["final_contigs"]}/{{sample}}.tkn',
        sample = glob_wildcards(os.path.join(checkpoint_output, "{sample}.fasta")).sample)
    return file_names


rule Verify_Filtered:
    input:
        f'{config["filtered_contigs"]}/{{sample}}.{config["extension"]}'
    output:
        f'{config["final_contigs"]}/{{sample}}.tkn'
    log:
        f'{config["log_folder"]}/Verify_Filtered__{{sample}}{config["taxon_id"]}.{config["log_id"]}.log'
    shell:
        """
        touch {output}
        """


rule Filter_All:
    input: grab_genomes
    output: f'{config["final_contigs"]}/complete.tkn'
    shell:
        """
        touch {output}
        """