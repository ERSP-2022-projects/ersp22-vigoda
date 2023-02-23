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
        

def treegen(species = 10,seqlen = 1000,p_mutate = 0.2 ,mutation_model = "jc69",seed=None):
    
    # set the path to the C++ executable and the command-line arguments
    executable_path ="C:/Users/yasha/Github/ersp22-vigoda/treegen/treegen.out"
    
    command = f"{executable_path} {species=} {seqlen=} {p_mutate=} {mutation_model=} "
    if(seed is not None):
        command.append(f"{seed=}")
    # run the C++ executable with the specified arguments
    #print(command)
    process = subprocess.run(command, shell=True, capture_output=True, text=True)
    # print the output of the C++ executable
    