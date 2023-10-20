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


def main(directory: str, extension: str, output: str) -> int:
    assert os.path.isdir(directory), "Invalid directory"
    files = glob.glob(f"{directory}/*.{extension}")
    assert len(files) > 1, "Only 1 sample"

    out_lines = '\n'.join(files)
    with open(output, 'w') as out:
        out.write(out_lines)
    return 0



def parse_args():
    parser = argparse.ArgumentParser(description="NCBI Parser")
    parser.add_argument("-d", "--Directory", help="Fasta storage locker", required=True)
    parser.add_argument("-e", "--Extension", help="Sequence extension", required=True)
    parser.add_argument("-o", "--Output", help="Output file name", required=True)
    return parser


if __name__ == "__main__":
    parser = parse_args()
    args = parser.parse_args()
    directory, extension, output = args.Directory, args.Extension, args.Output
    main(directory, extension, output)