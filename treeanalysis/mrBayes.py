import os

def generate_mrbayes_script(input_filename, output_filename):
    with open("mrbayes_template.nexus", "r") as template_file:
        template = template_file.read()
    script = template.replace("<input_filename>", input_filename)
    with open(output_filename, "w") as script_file:
        script_file.write(script)