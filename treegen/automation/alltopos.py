import subprocess
import pandas as pd
import os
sub_directories = next(os.walk('.'))[1]
print(sub_directories)
data = pd.DataFrame(columns=['Topology', 'Config', 'Seed', 'Swapseed', 'Iterations', 'Distance'])
for sub_directory in sub_directories:
    if (not os.path.exists(sub_directory + '/summary.csv')):
        subprocess.run(['sh', 'allconfig.sh'], cwd=sub_directory)
    new_data = pd.read_csv('./' + sub_directory + '/summary.csv', index_col=0)
    data = pd.concat([data, new_data], ignore_index=True)
data.to_csv("allsummary.csv")
