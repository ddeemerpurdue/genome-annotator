import glob

def main():
    folders = glob.glob("GCF_*")
    for folder in folders:
        genome_prefix = folder.replace('_', '').split('.')[0]
        for file in glob.glob(f'{folder}/*.fna'):
            if file.startswith('cds'):
                pass
            else:
                genome_loc = f'{folder}/{file}'