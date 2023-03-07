import pandas as pd
import os
import subprocess



data_dict = {
    'Topology': [],
    'Config': [],
    'Seed': [], 
    'Swapseed': [], 
    'Iterations': [], 
    'Distance': []
}
topology = os.path.basename(os.getcwd())
for config in ['allLong', 'allShort', 'internalLong', 'internalShort']:
    if (not os.path.exists(config + '/iter_summary.txt')):
        subprocess.run(['sh', 'runbatch.sh'], cwd=config + '/')
    iters = open(file = (config + '/iter_summary.txt'), mode='r')
    seeds = open(file = (config + '/seed_summary.txt'), mode='r')
    distances = open(file = (config + '/tree_distances.txt'), mode='r')
    for iter, seed_pair, distance in zip(iters.readlines(), seeds.readlines(), distances.readlines()):
        data_dict['Topology'].append(topology)
        data_dict['Config'].append(config)
        data_dict['Seed'].append(seed_pair.strip().split(" ")[0])
        data_dict['Swapseed'].append(seed_pair.strip().split(" ")[1])
        data_dict['Iterations'].append(iter.strip())
        data_dict['Distance'].append(distance.strip())
        print(topology, config, seed_pair.strip().split(" ")[0], seed_pair.strip().split(" ")[1], iter.strip(), distance.strip(), sep=",")
    iters.close()
    seeds.close()
    distances.close()
data = pd.DataFrame.from_dict(data_dict)
data.to_csv("summary.csv")