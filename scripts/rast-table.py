#!/usr/bin/env python3
import sys


def main(results, output):
    '''
    '''
    pseudo_acc = {}
    cnt = 0
    source = 'FigFams'

    with open(output, 'w') as out:
        header = f"gene_callers_id\tsource\taccession\tfunction\te_value\n"
        out.write(header)
        with open(results, encoding='cp1252') as file:
            line = file.readline()
            while line:
                values = line.strip().split('\t')
                gene_id = values[1]
                function = values[2]
                e_val = '1.0e-50'

                if function in pseudo_acc:
                    pass
                else:
                    cnt += 1
                    pseudo_acc[function] = f"FigFAM{str(cnt)}"

                writeline = f"{gene_id}\t{source}\t{pseudo_acc[function]}\t{function}\t{e_val}\n"
                out.write(writeline)
                line = file.readline()


results = sys.argv[1]
output = sys.argv[2]

if __name__ == '__main__':
    main(results, output)
