import os
import subprocess

def generate_mrbayes_script(input_filename, output_filename):
    with open("mrbayes_template.nexus", "r") as template_file:
        template = template_file.read()
    script = template.replace("<input_filename>", input_filename)
    log_filename = os.path.basename(input_filename)
    script = script.replace("<log_filename>", log_filename)
    with open(output_filename, "w") as script_file:
        script_file.write(script)
        

# def treegen(species = 10,seqlen = 1000,p_mutate = 0.2 ,mutation_model = "jc69",seed=None):
    
#     # set the path to the C++ executable and the command-line arguments
#     executable_path ="C:/Users/yasha/Github/ersp22-vigoda/treegen/treegen.out"
    
#     command = [executable_path, f"species={species}", f"seqlen={seqlen}", 
#         f"p_mutate={p_mutate}", f"mutation_model={mutation_model}"]    
#     if(seed is not None):
#         command.append(f"seed={seed}")
#     # run the C++ executable with the specified arguments
#     #print(command)
#     process = subprocess.run(command, capture_output=True, text=True)
#     #print the output of the C++ executable
#     print(f"{process.stdout=}")
#     print(f"{process.stderr=}")
    
def treegen(species=10, seqlen=1000, p_mutate=0.2, mutation_model="jc69", seed=None):
    # set the path to the C++ executable and the command-line arguments
    executable_path = "./treegen.out"
    result_directory = "C:/Users/yasha/Github/ersp22-vigoda/treegen/result"
    
    # set the working directory to the location of the C++ executable
    os.chdir(os.path.dirname(executable_path))
    
    command = [os.path.basename(executable_path), f"species={species}", f"seqlen={seqlen}", f"p_mutate={p_mutate}", f"mutation_model={mutation_model}"]
    if seed is not None:
        command.append(f"seed={seed}")
    # run the C++ executable with the specified arguments
    
    #print(" ".join(command))
    
    # process = subprocess.run(["echo","works!"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)    # print the output of the C++ executable
    process = subprocess.run(["./treegen.out"], capture_output=True, text=True)

    if process.stdout != "":
        print(f"{process.stdout=}")
    if process.stderr != "":
        print(f"{process.stderr=}")
    
    # change the working directory back to the original location
    os.chdir(os.path.dirname(os.getcwd()))


if __name__ == "__main__":
    treegen()