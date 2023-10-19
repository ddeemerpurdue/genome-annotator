#!/usr/bin/env python3
import sys


def main(results, translate, output):
    source = 'TIGRFAM'

    translater = {}
    with open(translate) as tran:
        line = tran.readline()
        while line:
            values = line.strip().split('\t')
            translater[values[0]] = values[1]
            line = tran.readline()

    with open(output, 'w') as out:
        header = f"gene_callers_id\tsource\taccession\tfunction\te_value\n"
        out.write(header)
        with open(results) as file:
            line = file.readline()
            while line:
                if line.startswith('#'):
                    pass
                else:
                    values = line.split()
                    gene_id = values[0]
                    tigr = values[2]
                    assert tigr.startswith('TIGR'), 'Invalid Tigr name!'
                    e_val = values[4]
                    function = translater[tigr]
                    writeline = f"{gene_id}\t{source}\t{tigr}\t{function}\t{e_val}\n"
                    out.write(writeline)
                line = file.readline()


if __name__ == '__main__':
    results = sys.argv[1]
    translate = sys.argv[2]
    output = sys.argv[3]
    main(results, translate, output)
