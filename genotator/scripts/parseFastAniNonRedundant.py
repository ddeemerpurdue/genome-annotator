import argparse
from collections import defaultdict
import os


def parse_fastani(in_file, threshold, basename=True):
    matches = []

    with open(in_file) as file:
        line = file.readline()
        while line:
            ref, quer, ani, mat, tot = line.strip().split('\t')
            if basename:
                ref = os.path.basename(ref)
                quer = os.path.basename(quer)
            if ref == quer:
                line = file.readline()
                continue
            if float(ani) >= threshold:
                matches.append((ref, quer))

            line = file.readline()
    print(matches)
    print(len(matches))
    return matches


def dfs(adj_list, visited, vertex, result, key):
    visited.add(vertex)
    result[key].append(vertex)
    for neighbor in adj_list[vertex]:
        if neighbor not in visited:
            dfs(adj_list, visited, neighbor, result, key)


def connected_components(matches):
    adj_list = defaultdict(list)
    for x, y in matches:
        adj_list[x].append(y)
        adj_list[y].append(x)

    result = defaultdict(list)
    visited = set()
    for vertex in adj_list:
        if vertex not in visited:
            dfs(adj_list, visited, vertex, result, vertex)
    return result


def write_groups(results, output):
    with open(output, 'w') as out:
        for result in results:
            values = results[result]
            out.write('\t'.join(values) + '\n')
    return 0


def parse_args():
    parser = argparse.ArgumentParser(description="Parser")
    parser.add_argument("-i", "--Input",
                        help="Input fastANI file.")
    parser.add_argument('-o', '--Output',
                        help='Output file to write to.',
                        default='fastAniParse.txt', required=False)
    parser.add_argument('-t', '--Threshold',
                        help='ANI threshold to filter redundancy',
                        default=99, required=False, type=float)
    return parser


if __name__ == '__main__':
    parser = parse_args()
    args = parser.parse_args()
    matches = parse_fastani(args.Input, args.Threshold)
    results = connected_components(matches)
    write_groups(results, args.Output)
