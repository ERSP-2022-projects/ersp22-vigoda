import os
import subprocess
import shutil


def generate_mrbayes_script(nexus_filepath, output_filepath=None):
    TEMPLATE_PATH = r"C:\Users\yasha\Github\ersp22-vigoda\treeanalysis\mrbayes_template.nexus"
    temp_path = os.path.normpath(TEMPLATE_PATH)

    nexus_filepath = os.path.normpath(nexus_filepath)
    if output_filepath == None:
        filename = os.path.basename(nexus_filepath)
        output_filepath = os.path.join(os.path.dirname(
            nexus_filepath), f"{filename}.mrbayes")
        # print(f"{(nexus_filepath,output_filepath)}")
    output_filepath = os.path.normpath(output_filepath)
    with open(temp_path, "r") as template_file:
        template = template_file.read()
    script = template.replace("<input_filename>", nexus_filepath)
    log_filename = os.path.basename(nexus_filepath)
    script = script.replace("<log_filename>", log_filename)
    with open(output_filepath, "w") as script_file:
        script_file.write(script)
    return output_filepath


def NexusToRun(nexus_filepath, script_filepath=None):
    if script_filepath == None:
        script_filepath = generate_mrbayes_script(nexus_filepath)

    input_filename = os.path.basename(nexus_filepath)
    mrbayes_filename = os.path.basename(script_filepath)
    # Extract the name before the underscore as the directory name
    dir_name = input_filename.split("_")[0]

    # Create the directory inside the analysis folder
    analysis_dir = os.path.join(
        "C:/Users/yasha/Github/ersp22-vigoda/treeanalysis/analysis", dir_name)
    os.makedirs(analysis_dir, exist_ok=True)

    if not os.path.exists(os.path.join(analysis_dir, input_filename)):
        shutil.move(nexus_filepath, os.path.join(analysis_dir, input_filename))
        nexus_filepath = os.path.join(analysis_dir, input_filename)
    else:
        print(f"{input_filename} has already been moved to {analysis_dir}")

    if not os.path.exists(os.path.join(analysis_dir, mrbayes_filename)):
        shutil.move(script_filepath, os.path.join(
            analysis_dir, mrbayes_filename))
        script_filepath = os.path.join(analysis_dir, mrbayes_filename)

    else:
        print(f"{mrbayes_filename} has already been moved to {analysis_dir}")

    output, error = run_mrbayes(script_filepath)
    return (output, error)


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
    old_directory = os.getcwd()
    # set the working directory to the location of the C++ executable
    os.chdir(treegen_directory)

    # run the C++ executable with the specified arguments
    command = [executable_path, f"species={species}",
               f"seqlen={seqlen}", f"p_mutate={p_mutate}",
               f"mutation_model={mutation_model}"]
    if seed is not None:
        command.append(f"seed={seed}")
    # process = subprocess.run(["echo","works!"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)    # print the output of the C++ executable
    process = subprocess.run(["./treegen.out"], capture_output=True, text=True)

    # change the working directory back to the original location
    os.chdir(old_directory)
    return (process.stdout, process.stderr)


if __name__ == "__main__":
    path = "C:\\Users\\yasha\\Github\\ersp22-vigoda\\treegen\\results\\44350_data.nex"
    run_mrbayes(path)
