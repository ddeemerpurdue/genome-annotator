#!/usr/bin/env python3
'''
map file is: fam-substrate-mapping-08252022.tsv
results file is: GCF000154385_dbsub.out
'''
import sys


def create_map(map_file):
    return_obj = {}
    print(f'Analyzing {map_file}')
    with open(map_file) as file:
        line = file.readline()
        line = file.readline()  # Skip header
        while line:
            values = line.strip().split('\t')
            if len(values) == 5:
                substr, sub, family, name, ec_num = values
            elif len(values) == 4:
                substr, sub, family, name = values
            else:
                print(f'Weird number of values...{len(values)}')
                print(f'{values}')
                line = file.readline()
                continue
            return_obj.setdefault(family, "")
            if return_obj[family] == "":
                return_obj[family] = name
            else:
                return_obj[family] = ';'.join([return_obj[family], name])

            line = file.readline()

    return return_obj


def main(results, output, map_file):
    mapper = create_map(map_file)
    source = 'CAZyme'

    with open(output, 'w') as out:
        header = f"gene_callers_id\tsource\taccession\tfunction\te_value\n"
        out.write(header)
        with open(results) as file:
            line = file.readline()
            line = file.readline()  # Skip header
            while line:
                values = line.split('\t')
                gene_id = values[5]
                cazyme_sub = values[0]
                cazyme = cazyme_sub.split('_')[0]
                diamond = values[4]
                e_val = values[7]
                function = mapper.get(cazyme)
                if not diamond == '-':
                    writeline = f"{gene_id}\t{source}\t{cazyme_sub}\t{function}\t{e_val}\n"
                    out.write(writeline)
                    line = file.readline()


if __name__ == '__main__':

    results = sys.argv[1]
    output = sys.argv[2]
    map_file = sys.argv[3]
    main(results, output, map_file)
