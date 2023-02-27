import os
import subprocess
import shutil


def generate_mrbayes_script(nexus_filepath, output_filepath=None, same_directory=True):
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


def run_mrbayes(script_filepath, directory=None):
    if directory == None:
        directory = os.path.dirname(script_filepath)
    old_path = os.getcwd()
    os.chdir(directory)
    command = ["mb", script_filepath]
    process = subprocess.run(command, capture_output=True, text=True)

    os.chdir(old_path)
    return (process.stdout, process.stderr)


def NexusToRun(nexus_filepath, script_filepath=None):
    nexus_filepath = os.path.normpath(nexus_filepath)

    if script_filepath == None:
        script_filepath = generate_mrbayes_script(nexus_filepath)
    script_filepath = os.path.normpath(script_filepath)

    output, error = run_mrbayes(script_filepath)

    input_filename = os.path.basename(nexus_filepath)
    prefix = input_filename.split("_")[0]

    # Create the directory inside the analysis folder
    analysis_dir = os.path.join(
        "C:/Users/yasha/Github/ersp22-vigoda/treeanalysis/analysis", prefix)
    os.makedirs(analysis_dir, exist_ok=True)

    src_dir = os.path.dirname(nexus_filepath)
    dest_dir = analysis_dir

    # Iterate over all files in the source directory
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

    return (output, error)


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
    print(NexusToRun(path))
