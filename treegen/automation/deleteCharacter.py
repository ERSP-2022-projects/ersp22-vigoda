import argparse
import random
# parses a nexus file for matrix, then for each line in matrix, replaces characters with '?' with probability 
def deleteCharacter(infilename, outfilename, probability):
    with open(infilename, 'r') as infile:
        with open(outfilename, 'w') as outfile:
            matrix = False
            for line in infile:
                if (line == 'matrix\n'):
                    matrix = True
                    outfile.write(line)
                elif (line == ';\n'):
                    matrix = False
                    outfile.write(line)
                elif (not matrix):
                    outfile.write(line)
                else:
                    charlist = []
                    for char in line:
                        charlist.append(char)
                        if char in ['A', 'T', 'G', 'C']:
                            replace = random.choices([True, False], weights=[probability, 1.0 - probability], k = 1)[0]
                            if (replace):
                                charlist[-1] = '?'
                    outfile.write(''.join(charlist))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', type = str)
    parser.add_argument('-o', '--outfile', type = str)
    parser.add_argument('-p', '--probability', type = float)
    args = parser.parse_args()
    infile = getattr(args, 'infile')
    outfile = getattr(args, 'outfile')
    prob = getattr(args, 'probability')
    deleteCharacter(infile, outfile, prob)


if __name__ == "__main__":
    main()

