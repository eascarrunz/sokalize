#!/usr/bin/env python3
import re
import sys
import argparse

try:
    import numpy as np
except ImportError:
    sys.exit("""You need numpy to run this script!
Install it from https://pypi.python.org/pypi/numpy
or run pip install numpy""")

parser = argparse.ArgumentParser(description="Recode a NEXUS multistate matrix into a binary matrix")
parser.add_argument('matrix_file',
                    help='Path of NEXUS file with the character matrix to recode. The matrix cannot be interleaved',
                    type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('-o', '--output', nargs='?', help='Name for the output file', type=argparse.FileType('w'),
                    default=sys.stdout, const=sys.stdout)
parser.add_argument('-f', '--format', help='Output format. If nexus is selected, detailed recoding information will be '
                                           'written out as a comment', choices=['nexus', 'phylip', 'tnt'],
                    default='nexus')
ctype_args = parser.add_mutually_exclusive_group()
ctype_args.add_argument('-O', '--ordered', help='Force all characters to be treated as ordered', action='store_true',
                        default=False)
ctype_args.add_argument('-U', '--unordered', help='Force all characters to be treated as unordered',
                        action='store_true',
                        default=False)
parser.add_argument('-a', '--remove_autapomorphic', help='Remove autapomorphic characters after recoding',
                    action='store_true', default=False)
parser.add_argument('-i', '--remove_invariant', help='Remove invariant characters after recoding', action='store_true',
                    default=False)
parser.add_argument('-m', '--remove_missing',
                    help='Remove characters without coded data (all missing or inapplicable) after recoding',
                    action='store_true', default=False)
args = parser.parse_args()

SYMBOLS = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
MIN_SEP = 2  # Minimum white space separation for writing out data
END_PATTERN = re.compile(r'''end\s*;''', re.IGNORECASE)
symbol_map = dict(zip(SYMBOLS, [2 ** i for i in range(len(SYMBOLS))]))
nsp = 0  # Number of species
sp_names = []  # Names of species
nchar_in = 0  # Number of characters
nstates = []  # Numbers of states
ctype = []  # Character types (unordered: 0, ordered: 1)
explicit_ctype = False
current_block = 0  # Which NEXUS block is currently being parsed. 0 = none; 1 = character, 2 = assumptions


def read_symbol(x):
    if x == "?":
        y = [0, 2]
    elif x == "-":
        y = [0, 1]
    else:
        y = [symbol_map[x], 0]
    return y


if __name__ == "__main__":
    # Store the raw data in memory as strings (data0)
    data0 = []
    in_matrix = False
    file_pos = 1
    # Parse NEXUS file
    for line in args.matrix_file:
        ORD_PATTERN = re.compile(r'''ord:\s*(\d+-*\s*)+[,;]''')
        UNORD_PATTERN = re.compile(r'''unord:\s*(\d+-*\s*)+[,;]''')
        if current_block == 0:  # Not in a block / in unsupported block
            if re.match(r'''^\s*begin\s+characters\s*;|^\s*begin\s+data\s*;''', line, re.IGNORECASE):
                current_block = 1
            elif re.match(r'''^\s*begin\s+assumptions\s*;''', line, re.IGNORECASE):
                current_block = 2

        elif current_block == 1:  # Data matrix block
            cursor = 1
            if re.match(r'''\s*matrix\s*''', line, re.IGNORECASE):
                in_matrix = True
                continue
            elif in_matrix:
                if re.match(r'''^\s*;\s*''', line):
                    in_matrix = False
                    continue
                if re.match(r'''^\s*$''', line):
                    continue
                try:
                    re.match(r'''^\s*[_0-9A-z\'\-.,/]+\s+[0-9A-z?\-\s{(\[})\]]+$''', line)
                except:
                    sys.stderr.write("Invalid format in line " + str(file_pos) + ":\n" + line)
                foo = re.search(r'''(^\s*[^ ]+)\s+(.+)''', line).groups()
                sp_names.append(foo[0].strip())
                raw_str = foo[1].strip()
                # Count characters
                in_poly = False
                count_char = 0
                count_space = 0
                for p in raw_str:
                    if p == "[" or p == "(" or p == "{":
                        in_poly = True
                        count_char += 1
                    elif p == "]" or p == ")" or p == "}":
                        in_poly = False
                    elif p == " ":
                        count_space += 1
                        continue
                    else:
                        count_char += 0 if in_poly else 1
                if cursor == 1:
                    nchar_in = count_char
                    data0.append(raw_str)
                else:
                    try:
                        nchar_in == count_char
                    except:
                        "Line {0} does not have the same number of characters as the preceding lines".format(
                            str(cursor))
                    data0.append(raw_str)
                cursor += 1
            elif END_PATTERN.match(line):
                current_block = 0
            if not explicit_ctype:
                ctype = [1] * nchar_in if args.ordered else [0] * nchar_in
        elif current_block == 2 and not args.unordered and not args.ordered:  # Assumptions block
            ctype_statement = re.search(r'''^\s*.*[un]?ord\s*:.*;''', line, re.IGNORECASE)
            if ctype_statement:
                explicit_ctype = True
                for part in ctype_statement.group().split(','):
                    char_index_aggregator = []
                    selection_statement = re.findall(r'''(\d+/*-*\d*)''', part)
                    for i in range(len(selection_statement)):
                        selection_statement[i] = re.sub(r'''/\d+''', '', selection_statement[i])
                        if re.search(r'''\s*-\s*''', selection_statement[i]):
                            selection_statement_parts = re.findall(r'''(\d+)''', part)
                            for x in range(int(selection_statement_parts[0]),
                                           int(selection_statement_parts[1]) + 1):
                                char_index_aggregator.append(x)
                        else:
                            char_index_aggregator.append(int(selection_statement[i]))
                    if re.match(r'''\s*ord''', part, re.IGNORECASE):
                        for x in char_index_aggregator:
                            ctype[int(x) - 1] = 1
                    elif re.match(r'''\s*unord''', part, re.IGNORECASE):
                        for x in char_index_aggregator:
                            ctype[int(x) - 1] = 0
        nsp = len(sp_names)
        file_pos += 1

    # Encode character state information with the symbol_map into data1
    data1 = np.zeros(shape=(nsp, nchar_in), dtype='int32')
    obs_type = np.zeros(shape=(nsp, nchar_in), dtype='int8')
    gor = False
    for sp, line in zip(range(nsp), data0):
        in_poly = False
        states_current = 0
        cursor = 0
        for p in line:
            if cursor == 13 and states_current > 3:
                gor = True
            if p == "[" or p == '(' or p == '{':
                in_poly = True
            elif p == ']' or p == ')' or p == '}':
                in_poly = False
                data1[sp, cursor] = states_current
                cursor += 1
                states_current = 0
            elif p == " ":
                continue
            else:
                if in_poly:
                    foo = read_symbol(p)
                    states_current += foo[0]
                    obs_type[sp, cursor] = foo[1]
                else:
                    states_current = 0
                    foo = read_symbol(p)
                    data1[sp, cursor] = foo[0]
                    obs_type[sp, cursor] = foo[1]
                    cursor += 1
    del data0
    nstates = np.bitwise_or.reduce(data1, 0).tolist()
    nstates = [x.bit_length() for x in nstates]
    # Identify binary characters
    is_binary = [False if x > 2 else True for x in nstates]

    # Check for polymorphism in binary characters, and recode it as '?'
    for i in range(nsp):
        for j in [y for y, x in enumerate(is_binary) if x]:
            if data1[i, j] == 3:
                data1[i, j] = 0
                obs_type[i, j] = 2
    # # Identify cells with observed data (not missing or inapplicable)
    data2 = []
    for nstate, char in zip(nstates, range(nchar_in)):
        if nstate > 2:
            tmp = (data1[:, char][:, np.newaxis] >> np.arange(nstate)[::-1] & 1)
        elif nstate == 2:
            tmp = (data1[:, char][:, np.newaxis] >> np.arange(nstate)[::-1] & 1)[:, 0]
        else:
            if data1[:, char].any():
                tmp = (data1[:, char][:, np.newaxis] >> np.arange(nstate)[::-1] & 0)[:, 0]
            else:
                tmp = np.zeros((nsp, 1), dtype='int8')
        data2.append(tmp)

    # data2 now contains a list of arrays, each array representing a character with one column per state
    # character arrays in data2 will be referred as char2 and each column
    # Character analysis
    outmat = np.zeros(shape=(nsp, 0), dtype="int8").astype(str)
    char_counter = 0
    char_report = ''
    for i, char2, nstate in zip(range(nchar_in), data2, nstates):
        add_char = True
        if nstates[i] <= 2:
            char_report += "Character " + str(i + 1) + " (binary) "
            tmp = char2[obs_type[:, i] == 0]
            # args.output.write(str(len(tmp)) + "\n")
            state_freqs = np.unique(tmp, return_counts=True)[1]
            all_missing = False
            if len(tmp) == 0:
                char_report += "has no observed states."
                all_missing = True
                if args.remove_missing:
                    add_char = False
            elif len(state_freqs) == 1:
                char_report += "is invariant."
                if args.remove_invariant:
                    add_char = False
            elif state_freqs[0] < 2 or state_freqs[1] < 2:
                char_report += "is autapomorphic."
                if args.remove_autapomorphic:
                    add_char = False
            else:
                char_report += "is parsimony-informative."
            if all_missing:
                bar = np.zeros((nsp, 1), dtype="int8").astype(str)
            else:
                bar = char2.astype(str)
            bar[obs_type[:, i] == 1] = "-"
            bar[obs_type[:, i] == 2] = "?"
            if add_char:
                outmat = np.c_[outmat, bar]
                char_counter += 1
                char_report += " Added to the recoded matrix as character " + str(char_counter) + "\n"
            else:
                char_report += " Dropped from the recoded matrix.\n"
                add_char = True
        else:
            if ctype[i]:
                for lin in range(char2.shape[0]):
                    msb_index = np.nonzero(char2[lin,])[0]
                    msb_index = None if len(msb_index) == 0 else msb_index[0]
                    if msb_index is not None:
                        char2[lin, msb_index:] = 1
            for char3 in range(nstate):
                tmp = char2[obs_type[:, i] == 0, char3]
                state_freqs = np.unique(tmp, return_counts=True)[1]
                char_report += "Character " + str(i + 1) + " state " + str(char3) + " "
                if len(state_freqs) == 1:
                    char_report += "is invariant."
                    if args.remove_invariant:
                        add_char = False
                elif state_freqs[0] < 2 or state_freqs[1] < 2:
                    char_report += "is autapomorphic."
                    if args.remove_autapomorphic:
                        add_char = False
                else:
                    char_report += "is parsimony-informative."
                    add_char = True
                bar = char2[:, char3].astype(str)
                bar[obs_type[:, i] == 1] = "-"
                bar[obs_type[:, i] == 2] = "?"
                if add_char:
                    outmat = np.c_[outmat, bar]
                    char_counter += 1
                    char_report += " Added to the recoded matrix as character " + str(char_counter) + "\n"
                else:
                    char_report += " Dropped from the recoded matrix.\n"
                    add_char = True
    nchar_out = outmat.shape[1]
    indent = 2 * " "
    foo = max([len(x) for x in sp_names]) + MIN_SEP
    if args.format == 'nexus':
        args.output.write(
            '#NEXUS\nBEGIN TAXA;\n' + indent + "DIMENSIONS NTAX=" + str(nsp) + ";\n" + indent + "TAXLABELS\n")
        for sp in sp_names:
            args.output.write(indent + indent + sp + '\n')
        args.output.write(indent + ';\nEND;\n')
        args.output.write(
            "BEGIN CHARACTERS;\n" + indent + "DIMENSIONS NCHAR=" + str(nchar_out) + " NTAX=" + str(nsp) + ";\n")
        args.output.write(indent + "FORMAT DATATYPE=STANDARD MISSING=? GAP=- SYMBOLS=\"0 1\";\n" + indent + "MATRIX\n")
        for i in range(nsp):
            padding = foo - len(sp_names[i])
            args.output.write(indent + sp_names[i] + " " * padding + "".join(outmat[i]) + "\n")
        args.output.write(indent + ";\nEND;\nBEGIN RECODING_REPORT;\n")
        args.output.write(indent + "Recoding a multistate matrix into binary:\n")
        if args.unordered:
            args.output.write(indent + "All characters were treated as unordered.\n")
        elif args.ordered:
            args.output.write(indent + "All characters were treated as ordered.\n")
        else:
            args.output.write(
                indent + "Characters were treated as ordered or unordered as specified in the ASSUMPTIONS block.\n")
        for line in char_report.splitlines():
            args.output.write(indent + line + "\n")
        args.output.write("END;")
    elif args.format == 'tnt':
        args.output.write("XREAD\n")
        args.output.write(str(nchar_out) + " " + str(nsp) + "\n")
        for i in range(nsp):
            padding = foo - len(sp_names[i])
            args.output.write(sp_names[i] + " " * padding + "".join(outmat[i]) + "\n")
        args.output.write(";\nPROC / ;\n")
    elif args.format == 'phylip':
        args.output.write(str(nsp) + " " + str(nchar_out) + "\n")
        for i in range(nsp):
            padding = foo - len(sp_names[i])
            args.output.write(sp_names[i] + " " * padding + "".join(outmat[i]) + "\n")
