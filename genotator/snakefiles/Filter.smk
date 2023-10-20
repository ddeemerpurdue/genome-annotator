import os
import shutil
print(f'Current working directory: {os.getcwd()}')
from datetime import datetime
tstamp = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
configfile: "../config.yaml"

include: "./Download.smk"

# ~/Documents/bioinf-programs/FastANI/build/
# For FastANI

rule filter_all:
    input:
        f'{config["verified_contigs"]}/complete.tkn'


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
        genome_list = f'{config["fastani"]}/{config["taxon_id"]}-Genomes.txt'
    log:
        f'{config["log_folder"]}/Filter_Ani__{config["taxon_id"]}.{config["log_id"]}.log'
    shell:
        """
        python ../scripts/parseFastAniNonRedundant.py -i {input} -t {params.threshold} -o {output.ani} &> {log}
        cut -f1 {output.ani} > {output.genome_list}
        """


rule Place_Final_Genomes:
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
    
