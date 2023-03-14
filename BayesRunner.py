import os
import subprocess
import shutil
import re
import numpy as np
import time
from contextlib import contextmanager
from datetime import datetime

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


def treegen(species=10, seqlen=1000, p_mutate=0.2,
            mutation_model="jc69", seed=None):
    # set the path to the C++ executable and the command-line arguments
    executable = r"./treegen.out"
    treegen_directory = r"./treegen/"
    results_directory = r"./treegen/results"
    
    # run the C++ executable with the specified arguments
    command = [executable, f"species={species}",
               f"seqlen={seqlen}", f"p_mutate={p_mutate}",
               f"mutation_model={mutation_model}"]
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
                     species=10, seqlen=1000, p_mutate=0.2,
                    mutation_model="jc69", seed=None):
    pass
    # if filepath1 != None:
        
    


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
    analysis_dir = os.path.normpath(analysis_dir)
    os.makedirs(analysis_dir, exist_ok=True)
    src_dir = os.path.dirname(nexus_filepath)
    dest_dir = analysis_dir
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
    log_filename = os.path.basename(nexus_filepath)
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


def generateAndRun(**kwargs):
    nexus_filepath = kwargs.get('nexus_filepath', None)
    if nexus_filepath == None:
        treegen_params = {
            'species': kwargs.get('species', None),
            'seqlen': kwargs.get('seqlen', None),
            'p_mutate': kwargs.get('p_mutate', None),
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
    output, error, dir_path = run_mrbayes(script_filepath)
    return dir_path


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
    extract_from_directory = os.path.normpath(extract_from_directory)
    table = extract_data_from_directory(extract_from_directory)
    print_timestamped(f"{table=}")


def extract_data_from_directory(target_directory, output_filename="_results.csv"):
    # create an empty list to store the data
    data_list = []

    # loop over all files in the specified directory
    for file_name in os.listdir(target_directory):
        if file_name.endswith('.mcmc'):
            file_path = os.path.join(target_directory, file_name)

            # extract the relevant data and add it to the data_list
            prefix, gen_num, std_dev = extract_data(file_path)
            data_list.append([prefix, gen_num, std_dev])

    # convert the data_list to a numpy array
    data_array = np.array(data_list)

    # append the .csv extension to the output filename if it's not already present
    if not output_filename.endswith('.csv'):
        output_filename += '.csv'

    # write the data_array to a CSV file in the target_directory
    output_file_path = os.path.join(target_directory, output_filename)
    header = "prefix,generation number,ASDSF"
    np.savetxt(output_file_path, data_array,
               delimiter=',', fmt='%s', header=header)

    # return the data_array
    return data_array


def extract_data(file_path):
    # extract the prefix number tag from the file name
    file_path = os.path.normpath(file_path)
    prefix = re.findall(r'^(\d+)_', os.path.basename(file_path))[0]

    # open the file and read its contents
    with open(file_path, 'r') as f:
        contents = f.read()

    # find the last row of the table
    last_row = contents.split('\n')[-2]

    # extract the generation number and StdDev value
    gen_num = last_row.split('\t')[0]
    std_dev = last_row.split('\t')[-1]

    return (prefix, gen_num, std_dev)






















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


def main():
    seqlen_wrapper("nCharsMixingTime-20000", seqlen=20000, numSamples=10)

    print_timestamped("Program completed")

def fake_main():
    seqlen_wrapper("misc", seqlen=1000, numSamples=3)


if __name__ == "__main__":
    main()