import os
import subprocess
import shutil
import re
import numpy as np
import time
from contextlib import contextmanager
from datetime import datetime
import random
import dendropy
from dendropy.calculate import treecompare
@contextmanager
def cwd(path):
    oldpwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(oldpwd)



def print_timestamped(*args, **kwargs):
    now = datetime.now()
    timestamp = now.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}]", *args, **kwargs)
'''
Generates a single tree topology and the accompanying character sequence by calling treegen.out. Results in a .nex 
file corresponding to the character sequence and a .txt file corresponding to the generated tree topology in Newick format.
The generating tree is also included in a 'trees' block at the end of the NEXUS file
Args:
    species (int): number of species in tree
    seqlen (int): number of characters to include in generated genetic sequence
    b_length (float): float between 0 and 1 indicating the desired branch length in the tree
        (0.1 by default, which is approximately equal to p_mutate of ~0.082 using default JC69 rate matrix in NewickTree)
    mutation_model (string): type of mutation model to apply to treegen.out (not yet implemented)
    seed (int): seed to generate tree with, randomly set by default
    results_directory (string): name of directory within files to store nexus file in

Returns:
    (string): path to generated nexus file with genetic character sequences
'''
def treegen(species=10, seqlen=1000, b_length=0.2,
            mutation_model="jc69", seed=None, results_directory = None):
    # set the path to the C++ executable and the command-line arguments
    executable = r"./treegen.out"
    treegen_directory = r"./treegen/"

    if (results_directory == None):
        relative_results = ""
        results_directory = r"./treegen/results"

    else:
        relative_results = results_directory
        results_directory = os.path.join("./treegen/results", results_directory)
    if (not os.path.isdir(results_directory)):
        os.makedirs(results_directory)
    
    # run the C++ executable with the specified arguments
    command = [executable, 
                f"species={species}",
                f"seqlen={seqlen}", 
                f"b_length={b_length}",
                f"mutation_model={mutation_model}",
                f"filepath={relative_results}"]
    if seed is not None:
        command.append(f"seed={seed}")
        
    #change directories and run 
    with cwd(treegen_directory):
        process = subprocess.run(
            command, capture_output=True)
        
    stdout = (process.stdout).decode()
    stderr = (process.stderr).decode()
        
    # get the prefix of the newly created NEXUS file
    prefix = (stdout).split('\n')[-2].split()[-1]

    for filename in os.listdir(results_directory):
        if filename.startswith(prefix) and filename.endswith('.nex'):
            return os.path.join(results_directory, filename)

    raise FileNotFoundError(
        f"The file with prefix: '{prefix}' does not exist in results directory")

def generate_Mixture(filepath1 = None,filepath2=None,
                     species=10, seqlen=1000, b_length=0.2,
                    mutation_model="jc69", seed=None):
    pass
    # if filepath1 != None:
        
    

'''
Generate NEXUS file containing the following:
    - (data block) character sequence data to analyze
    - (trees block) original generating tree(s)
    - (mrbayes block) script to run MrBayes with given parameters
Args:
    nexus_filepath (string): path to NEXUS file to source character sequence data from
    target_directory (string): directory to store MrBayes script and analysis results in
    output_filepath (string): filepath to store final MrBayes script in (this is actually a file)

'''
def generate_mrbayes_script(nexus_filepath, target_directory=None,
                            output_filepath=None, ngen=5000,
                            samplefreq=100, nchains=1, nruns=2,
                            printfreq=100, diagnfreq=100,
                            stopval=0.00001, burnin=250, nst=6,
                            rates="gamma", addToNexusFile=True):
    input_filename = os.path.basename(nexus_filepath)
    prefix = input_filename.split("_")[0]

    # Create the directory inside the analysis folder
    if target_directory == None:
        analysis_dir = os.path.normpath(
            "./treeanalysis/analysis/misc")
    else:
        analysis_dir = target_directory
    analysis_dir = os.path.join(os.path.normpath(analysis_dir), 'data')
    os.makedirs(analysis_dir, exist_ok=True)
    src_dir = os.path.dirname(nexus_filepath)
    dest_dir = analysis_dir
    if (not os.path.exists(dest_dir)):
        os.makedirs(dest_dir)
    for file_name in os.listdir(src_dir):
        # Check if the file starts with the prefix "input_filename"
        if file_name.startswith(prefix):
            # If the file starts with the prefix, construct the source and destination file paths
            src_path = os.path.join(src_dir, file_name)
            dest_path = os.path.join(dest_dir, file_name)
            # Check if the file already exists in the destination directory
            if os.path.exists(dest_path):
                print_timestamped(
                    f"File '{file_name}' already exists in the destination directory, skipping...")
            else:
                # If the file does not exist in the destination directory, move it there
                shutil.move(src_path, dest_path)
    nexus_filepath = os.path.join(analysis_dir, input_filename)
    if output_filepath != None:
        output_filename = os.path.basename(output_filepath)
        output_filepath = os.path.join(analysis_dir, output_filename)
    TEMPLATE_PATH = r"./treeanalysis/mrbayes_template.nexus"

    temp_path = os.path.normpath(TEMPLATE_PATH)

    nexus_filepath = os.path.normpath(nexus_filepath)
    if output_filepath == None:
        filename = os.path.basename(nexus_filepath)
        output_filepath = os.path.join(os.path.dirname(
            nexus_filepath), f"{filename}.mrbayes")
    output_filepath = os.path.normpath(output_filepath)
    with open(temp_path, "r") as template_file:
        template = template_file.read()
    script = template
    if addToNexusFile == False:
        script = script.replace("<input_filename>", nexus_filepath)
    log_filename = os.path.splitext(os.path.basename(nexus_filepath))[0] + "_log.txt"
    script = script.replace("<log_filename>", log_filename)
    script = script.replace("<ngen>", str(ngen))
    script = script.replace("<samplefreq>", str(samplefreq))
    script = script.replace("<nchains>", str(nchains))
    script = script.replace("<nruns>", str(nruns))
    script = script.replace("<printfreq>", str(printfreq))
    script = script.replace("<diagnfreq>", str(diagnfreq))
    script = script.replace("<burnin>", str(burnin))
    script = script.replace("<stopval>", str(stopval))
    script = script.replace("<nst>", str(nst))
    script = script.replace("<rates>", str(rates))

    if addToNexusFile:
        with open(nexus_filepath, "a") as nexus_file:
            nexus_file.write(
                "\n" + script.replace("execute <input_filename>;", "") + "\n")
            output_filepath = nexus_filepath
    else:
        with open(output_filepath, "w") as script_file:
            script_file.write(script)
    return output_filepath

'''
runs MrBayes from a given script
Args:
    script_filepath (string): Path to MrBayes script to be run. Should include character data, generating tree(s), and MrBayes commands
Returns:
    (tuple):
        - (string): standard output of MrBayes process
        - (string): standard error of MrBayes process
        - (string): file path to directory containing MrBayes results
'''


def run_mrbayes(script_filepath, directory=None):
    script_filepath = os.path.normpath(os.path.abspath(script_filepath))
    assert os.path.isfile(script_filepath), f"script_filepath doesn't exist"
    
    if directory == None:
        directory = os.path.dirname(script_filepath)
    if (len(script_filepath) >= 99):
        short_filepath = os.path.relpath(script_filepath,start = directory)
        script_filepath = short_filepath
    directory = os.path.normpath(directory)
    
    #mrBayes only takes filepaths of max length 100 
    #so this replaces filepaths with a short relative filepath
        
    command = ["mb", script_filepath]    
    with cwd(directory):  
        process = subprocess.run(
            command, capture_output=True)

    return (process.stdout, process.stderr, directory)

'''
Generates/takes in a single NEXUS file with genetic character sequenes and runs MrBayes on it.
Args:
    **kwargs: see below
Returns:
    (string): filepath to MrBayes output
'''
def generateAndRun(**kwargs):

    nexus_filepath = kwargs.get('nexus_filepath', None)

    if nexus_filepath == None:
        treegen_params = {
            'species': kwargs.get('species', None),
            'seqlen': kwargs.get('seqlen', None),
            'b_length': kwargs.get('b_length', None),
            'mutation_model': kwargs.get('mutation_model', None),
            'seed': kwargs.get('seed', None),
        }
        nexus_filepath = treegen(
            **{k: v for k, v in treegen_params.items() if v is not None})

    generate_mrbayes_script_params = {
        'nexus_filepath': nexus_filepath,
        'output_filepath': kwargs.get('output_filepath', None),
        'target_directory': kwargs.get('target_directory', None),
        'ngen': kwargs.get('ngen', None),
        'samplefreq': kwargs.get('samplefreq', None),
        'stopval': kwargs.get('stopval', None),
        'nchains': kwargs.get('nchains', None),
        'nruns': kwargs.get('nruns', None),
        'printfreq': kwargs.get('printfreq', None),
        'diagnfreq': kwargs.get('diagnfreq', None),
        'burnin': kwargs.get('burnin', None),
        'nst': kwargs.get('nst', None),
        'rates': kwargs.get('rates', None),
        'addToNexusFile': kwargs.get('addToNexusFile', None),
    }
    script_filepath = generate_mrbayes_script(
        **{k: v for k, v in generate_mrbayes_script_params.items() if v is not None})
    print_timestamped('script filepath:', script_filepath)
    output, error, dir_path = run_mrbayes(script_filepath)
    return dir_path

# 
def CreateAndAnalyze(numSamples, test_dir_name="misc", **kwargs):
    analysis_directory = r"./treeanalysis/analysis"
    target_directory = os.path.join(analysis_directory, test_dir_name)
    if not os.path.exists(target_directory):
        os.makedirs(target_directory)
        print_timestamped(f"Created directory: {target_directory}")

    kwargs["target_directory"] = target_directory

    for i in range(numSamples):
        generateAndRun(**kwargs)
        print_timestamped(f"Analyzing sample {i+1}")
    print_timestamped(f"Analyzing {numSamples} samples completed")
    return target_directory


def extract_results(extract_from_directory):
    output_directory = extract_from_directory
    extract_from_directory = os.path.join(os.path.normpath(extract_from_directory), 'data')
    table = extract_data_from_directory(extract_from_directory, output_directory=output_directory)
    print_timestamped(f"{table=}")

# (WIP) needs to include branch length, percentage of characters removed, ntaxa, nchars, seed, swap seed, distance from correct tree
# need some way to ensure that same data points are not recorded multiple times
# primary
def extract_data_from_directory(target_directory, output_directory = None, output_filename="_results.csv"):
    if output_directory == None:
        output_directory = target_directory
    # create an empty list to store the data
    data_list = []

    # loop over all files in the specified directory
    for file_name in os.listdir(target_directory):
        if file_name.endswith('.mcmc'):
            # prefix is filename before extension
            file_name_prefix = file_name.split('.')[0]
            # seed is stored as first sequence before any percentage markers (e.g. 1234_50.nex
            # has a seed of 1234 and a percentage (removed) of 50)
            if (file_name_prefix.find('_') != -1):
                file_name_seed = file_name_prefix[ : file_name_prefix.find('_')]
                percentage = file_name_prefix[file_name_prefix.find('_') + 1 : ]
            else:
                file_name_seed = file_name_prefix
                percentage = str(100)

            mcmc_file_path = os.path.join(target_directory, file_name)
            trprobs_file_path = os.path.join(target_directory, file_name_prefix + '.nex.trprobs')
            log_file_path = os.path.join(target_directory, file_name_prefix + '_log.txt')
            tree_file_path = os.path.join(target_directory, file_name_seed + '_tree.txt')
            # extract the relevant data and add it to the data_list
            prefix, gen_num, std_dev = extract_data_gen(mcmc_file_path)
            predicted_tree = extract_data_predicted_tree(trprobs_file_path)
            generating_tree, internal_b_length, terminal_b_length = extract_data_generating_tree(tree_file_path)
            tns = dendropy.TaxonNamespace()
            p_tree = dendropy.Tree.get(data = predicted_tree, schema='newick', taxon_namespace=tns)
            g_tree = dendropy.Tree.get(data = generating_tree, schema='newick', taxon_namespace = tns)
            #implements unweighted Robinson-Foulds distance
            distance = treecompare.symmetric_difference(p_tree, g_tree)
            seed, swapseed = extract_data_seed(log_file_path)
            data_list.append([file_name_prefix, percentage, internal_b_length, terminal_b_length, seed, swapseed, gen_num, std_dev, distance])

    # convert the data_list to a numpy array
    data_array = np.array(data_list)

    # append the .csv extension to the output filename if it's not already present
    if not output_filename.endswith('.csv'):
        output_filename += '.csv'

    # write the data_array to a CSV file in the target_directory
    output_file_path = os.path.join(output_directory, output_filename)
    header = "prefix, percent_removed, internal_length, terminal_length, seed, swapseed, gen_num, ASDSF, distance"
    np.savetxt(output_file_path, data_array,
               delimiter=',', fmt='%s', header=header)

    # return the data_array
    return data_array
# parses generating tree topology and branch lengths from _tree.txt file
def extract_data_generating_tree(file_path):
    with open(file_path, 'r') as file:
        content = file.readlines()
    tokens = ['topology = ', 'internal branch length = ', 'terminal branch length = ']
    for index, token in enumerate(tokens):
        if (content[index].find(token) == -1):
            raise Exception(file_path + " doesn't contain at least one of topology, internal branch length, or terminal branch length")
    topology = content[0][len(tokens[0]) : ].strip()
    internal_b_length = content[1][len(tokens[1]) : ].strip()
    terminal_b_length = content[2][len(tokens[2]) : ].strip()
    return topology, internal_b_length, terminal_b_length

# parses highest probability tree from .trprobs file
def extract_data_predicted_tree(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if (line.find('tree tree_1') != -1):
                return line[line.find('('):]
    return ''

# parses seed data from _log.txt file
def extract_data_seed(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if (line.find('Seed = ') != -1):
                seed = int(line[line.find('Seed = ') + 7 : ])
            if (line.find('Swapseed = ') != -1):
                swapseed = int(line[line.find('Swapseed = ') + 11: ])
            
    return seed, swapseed
# extract data from .mcmc file
# note that prefix contains just the generating seed and not any percentage values
def extract_data_gen(file_path):
    # extract the prefix number tag from the file name
    print_timestamped(file_path)
    file_path = os.path.normpath(file_path)
    prefix = re.findall(r'^(\d+)', os.path.basename(file_path))[0]

    # open the file and read its contents
    with open(file_path, 'r') as f:
        contents = f.read()

    # find the last row of the table
    last_row = contents.split('\n')[-2]

    # extract the generation number and StdDev value
    gen_num = last_row.split('\t')[0]
    std_dev = last_row.split('\t')[-1]

    return (prefix, gen_num, std_dev)

def remove_characters(infilename, outfilename, probability):
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
    return outfilename
'''

'''
def missing_character_wrapper(num_samples, percentages, test_dir_name = None, species = [10], seqlen = [10000], b_length = 0.2, character_directory=None):

    if (test_dir_name == None):
        test_dir_name = str(num_samples) + "_"
        for parameter, label in [(percentages, 'p'), (species, 'sp'), (seqlen, 'sl')]:
            test_dir_name += label + '-'
            for value in parameter:
                test_dir_name += str(value)
            test_dir_name += '_'
        test_dir_name += 'b-' + str(b_length)
    mrbayes_params = {
        "ngen": 1000000,
        "samplefreq": 100,
        "printfreq": 100,
        "diagnfreq": 100,
        "stopval": 0.01
    }
    analysis_directory = r"./treeanalysis/analysis"
    target_directory = os.path.join(analysis_directory, test_dir_name)
    sequences_generated = []
    if not os.path.exists(target_directory):
        os.makedirs(target_directory)
        print_timestamped(f"Created directory: {target_directory}")
    # generate topologies using treegen
    for length in seqlen:
        for ntaxa in species:
            for i in range(num_samples):
                character_filepath = treegen(species = ntaxa, seqlen = length, b_length = b_length, results_directory=character_directory)
                results_directory = os.path.dirname(character_filepath)
                filename_prefix = os.path.splitext(os.path.basename(character_filepath))[0]

                for percentage in percentages:
                    percentage = round(percentage, 2)
                    suffix = '_' + '{:02d}'.format(int(100 * percentage)) + '.nex'
                    output_filename = filename_prefix + suffix
                    sequences_generated.append(remove_characters(character_filepath, os.path.join(results_directory, output_filename), percentage))
    for sequence in sequences_generated:
        generateAndRun(nexus_filepath=sequence, target_directory = target_directory, **mrbayes_params)
    # for each nexus file, remove characters

    # run mrbayes on all files
    # extract data from all runs into csv
    extract_results(target_directory)

    pass
def seqlen_wrapper(test_dir_name, seqlen, numSamples):
    treegen_params = {
        "species": 20,
        "seqlen": seqlen,
    }
    mrbayes_params = {
        "ngen": 5000,
        "samplefreq": 100,
        "printfreq": 100,
        "diagnfreq": 100,
        "stopval": 0.0001
    }
    all_params = {**treegen_params, **mrbayes_params}
    test_directory = CreateAndAnalyze(
        test_dir_name=test_dir_name, numSamples=numSamples, **all_params)
    time.sleep(60)
    print_timestamped(f"Extracting Results in {test_dir_name}")
    print_timestamped(extract_results(test_directory))


def species_wrapper(test_dir_name, species, numSamples, seqlen=10000):
    treegen_params = {
        "species": species,
        "seqlen": seqlen,
    }
    mrbayes_params = {
        "ngen": 50000,
        "samplefreq": 100,
        "printfreq": 100,
        "diagnfreq": 100,
        "stopval": 0.0001
    }
    all_params = {**treegen_params, **mrbayes_params}
    test_directory = CreateAndAnalyze(
        test_dir_name=test_dir_name, numSamples=numSamples, **all_params)
    time.sleep(60)
    print_timestamped(f"Extracting Results in {test_dir_name}")
    print_timestamped(extract_results(test_directory))

'''
takes in a directory of nexus files and a list of percentages
to remove. Goes through each nexus file and creates a copy with the prefix
fix extract data to strip everything before .nex instead of using digits
Args:
    source_directory (string): directory containing nexus files with no missing characters
    output_directory (string): directory to output files with characters removed
        files will be named [original]_[percentage].nex
    percentages (list or listlike): list of floats between 0 and 1 representing what percentage of characters to remove
Returns:
    (list): list of paths to files created
'''
# def remove_characters_wrapper(source_directory, output_directory, percentages):
#     files_created = []
#     for file in os.listdir(source_directory):
#         if (os.path.splitext(file)[-1] == '.nex'):
#             prefix = os.path.splitext(file)[0]
#             for percentage in percentages:
#                 percentage = round(percentage, 2)
#                 suffix = '_' + '{:02d}'.format(int(100 * percentage))
#                 remove_characters(os.path.join(source_directory, file), os.path.join(output_directory, prefix + suffix + '.nex'), percentage)
#                 files_created.append(os.path.join(prefix + suffix))
#     return files_created


def main():
    missing_character_wrapper(num_samples=25, percentages=[0.0,0.2,0.4,0.6,0.8], species = [100], seqlen = [10000], b_length = 0.2, character_directory="test2")

def fake_main():
    seqlen_wrapper("misc", seqlen=1000, numSamples=3)


if __name__ == "__main__":
    main()