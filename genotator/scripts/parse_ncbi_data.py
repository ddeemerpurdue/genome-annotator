'''
Python script that parses NCBI zip file from the
Datasets CLI
'''

import glob
import os
import pathlib
import shutil
import zipfile

import argparse


def unzip_data(file, base_location):
    # if os.path.isdir(base_location):
        # print('Base location exists...skipping unzipping')
        # return 0

    with zipfile.ZipFile(file, 'r') as zip_ref:
        zip_ref.extractall(os.path.join(base_location, "ncbi_data"))
    return 0


def crawl_ncbi(base_location):
    fastas = {}
    ncbi_base = os.path.join(base_location, "ncbi_data/ncbi_dataset/data/")
    if os.path.isdir(ncbi_base):
        folders = glob.glob(f'{ncbi_base}/GCF_*')
        for folder in folders:
            bn_folder = os.path.basename(folder)
            genome_prefix = bn_folder.replace('_', '').split('.')[0]
            genome_version = bn_folder.replace('_', '').split('.')[1]  # Add versioning
            for file in glob.glob(f'{folder}/*.fna'):
                if file.startswith('cds'):
                    pass
                else:
                    bn = os.path.basename(file)
                    if bn in fastas:
                        if int(genome_version) > int(fastas[bn][2]):
                            fastas[bn] = (bn, file, genome_version)
                    else:
                        fastas[bn] = (bn, file, genome_version)

    return fastas


def main(args):
    input_file = pathlib.Path(args.Input)
    base_location = input_file.parents[0]
    print(f'{input_file=}')
    print(f'{base_location=}')
    unzip_data(args.Input, base_location)
    fastas = crawl_ncbi(base_location)
    if len(fastas) == 0:
        print('No fasta files found...')
        return 1

    print(f'Found {len(fastas)} fasta files.')
    fasta_locker = args.Directory
    
    print(f'Fasta locker: {fasta_locker}')

    if not os.path.exists(fasta_locker):
        print(f'Making directory: {fasta_locker}')
        os.mkdir(fasta_locker)
    else:
        print(f'Skipping making the directory')
    for fasta in fastas:
        data = fastas[fasta]
        dest = os.path.join(fasta_locker, data[0])
        # dest = f"{fasta_locker}/"
        print(f'Copying:\n{data[1]}\n--->\n{dest}\n\n')
        shutil.copy(data[1], dest)
    return 0



def parse_args():
    parser = argparse.ArgumentParser(description="NCBI Parser")
    parser.add_argument("-i", "--Input", help="Input .zip file", required=True)
    parser.add_argument("-d", "--Directory", help="Fasta storage locker",
                        required=True)
    return parser


if __name__ == "__main__":
    parser = parse_args()
    args = parser.parse_args()
    main(args)