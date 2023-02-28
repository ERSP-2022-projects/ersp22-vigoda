import os
import subprocess
import shutil


def generate_mrbayes_script(nexus_filepath, output_filepath=None, ngen=50000, samplefreq=100, nchains=1, nruns=2, printfreq=1000, diagnfreq=1000, burnin=250, nst=6, rates="gamma", addToNexusFile=True):
    input_filename = os.path.basename(nexus_filepath)
    prefix = input_filename.split("_")[0]

    # Create the directory inside the analysis folder
    analysis_dir = os.path.join(
        "C:/Users/yasha/Github/ersp22-vigoda/treeanalysis/analysis", prefix)
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
                print(
                    f"File '{file_name}' already exists in the destination directory, skipping...")
            else:
                # If the file does not exist in the destination directory, move it there
                shutil.move(src_path, dest_path)
    nexus_filepath = os.path.join(analysis_dir, input_filename)
    if output_filepath != None:
        output_filename = os.path.basename(output_filepath)
        output_filepath = os.path.join(analysis_dir, output_filename)
    TEMPLATE_PATH = r"C:\Users\yasha\Github\ersp22-vigoda\treeanalysis\mrbayes_template.nexus"

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
    if directory == None:
        directory = os.path.dirname(script_filepath)
    old_path = os.getcwd()
    os.chdir(directory)
    command = ["mb", script_filepath]
    process = subprocess.run(command, capture_output=True, text=True)

    os.chdir(old_path)
    return (process.stdout, process.stderr)


def treegen(species=10, seqlen=1000, p_mutate=0.2, mutation_model="jc69", seed=None):
    # set the path to the C++ executable and the command-line arguments
    executable_path = "./treegen.out"
    treegen_directory = "C:/Users/yasha/Github/ersp22-vigoda/treegen"
    results_directory = "C:/Users/yasha/Github/ersp22-vigoda/treegen/results"
    old_directory = os.getcwd()
    # set the working directory to the location of the C++ executable
    os.chdir(treegen_directory)

    # run the C++ executable with the specified arguments
    command = [executable_path, f"species={species}",
               f"seqlen={seqlen}", f"p_mutate={p_mutate}",
               f"mutation_model={mutation_model}"]
    if seed is not None:
        command.append(f"seed={seed}")
    process = subprocess.run(["./treegen.out"], capture_output=True, text=True)

    # get the prefix of the newly created NEXUS file
    prefix = (process.stdout).split('\n')[-2].split()[-1]

    # change the working directory back to the original location
    os.chdir(old_directory)

    for filename in os.listdir(results_directory):
        if filename.startswith(prefix) and filename.endswith('.nex'):
            return os.path.join(results_directory, filename)

    raise FileNotFoundError(
        f"The file with prefix: '{prefix}' does not exist in results directory")


def generateAndRun(**kwargs):
    treegen_params = {
        'species': kwargs['species'],
        'seqlen': kwargs['seqlen'],
        'p_mutate': kwargs['p_mutate'],
        'mutation_model': kwargs['mutation_model'],
        'seed': kwargs['seed'],
    }
    nexus_filepath = treegen(**treegen_params)
    generate_mrbayes_script_params = {
        'nexus_filepath': nexus_filepath,
        'ngen': kwargs['ngen'],
        'samplefreq': kwargs['samplefreq'],
        'nchains': kwargs['nchains'],
        'nruns': kwargs['nruns'],
        'printfreq': kwargs['printfreq'],
        'diagnfreq': kwargs['diagnfreq'],
        'burnin': kwargs['burnin'],
        'nst': kwargs['nst'],
        'rates': kwargs['rates'],
        'addToNexusFile': kwargs['addToNexusFile'],
    }

    script_filepath = generate_mrbayes_script(**generate_mrbayes_script_params)
    return run_mrbayes(script_filepath)


if __name__ == "__main__":
    # path = "C:\\Users\\yasha\\Github\\ersp22-vigoda\\treegen\\results\\45299_data.nex"
    # print(NexusToRun(path))
    print(generateAndRun())
