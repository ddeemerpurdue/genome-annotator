#!/usr/bin/env python

import argparse
import os
import subprocess
import sys


snakefile_key = {
    "download": "../genotator/snakefiles/Download.smk",
    "filter": "../genotator/snakefiles/Filter.smk",
    "annotate": "../genotator/snakefiles/Annotate.smk"
}


def build_cli(command, snakefile):
    default_config = os.path.join(os.getcwd(), "../genotator/config.yaml")
    cli_prefix = ["snakemake", "-s"]
    cli_suffix = command.copy()
    if not "--configfile" in cli_suffix:
        cli_suffix.extend(["--configfile", default_config])
    cli_suffix[0] = snakefile
    return cli_prefix + cli_suffix

def main(command):
    '''
    Todo
    '''
    help_flags = ['--help', '-h', '-H']
    list_flags = ['--list', '-L', '--List']
    if len(command) < 1:
        print("Placeholder for Genotator main help", '\n')
    else:
        program = command[0]
        if program in help_flags:
            print('Help flag was raise. TODO')
            return 0
        if program in list_flags:
            print('List flag was raised. List all options TODO')
            return 0

        if program in snakefile_key:
            snakefile = snakefile_key[program]
            full_cli = build_cli(command, snakefile)
    
            print(f'Running {full_cli=}')
            subprocess.run(full_cli, check=True)

            return 0

        return 1

if __name__ == '__main__':
    arguments = sys.argv[1:]
    main(arguments)