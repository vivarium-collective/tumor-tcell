import argparse
import os
import pickle
from vivarium.core.emitter import (
    data_from_database,
    DatabaseEmitter,
    remove_units, deserialize_value
)
#from vivarium.library.units import units, remove_units


OUT_DIR = 'out/analysis/'


def access():
    # parse
    parser = argparse.ArgumentParser(description='access data from db')
    parser.add_argument(
            'experiment_id',
            type=str,
            default=False)
    args = parser.parse_args()
    experiment_id = args.experiment_id
    # mongo client
    config = {
        'host': '{}:{}'.format('localhost', 27017),
        'database': 'simulations'}
    emitter = DatabaseEmitter(config)
    db = emitter.db
    # access
    data, sim_config = data_from_database(
        experiment_id, db)

    return data, experiment_id, sim_config

def run_analysis(data, experiment_id='tumor_tcell', sim_config={}):
    experiment_out_dir = OUT_DIR + str(experiment_id)
    os.makedirs(experiment_out_dir, exist_ok=True)

    for i in data.keys():
        del data[i]['_id']
        del data[i]['experiment_id']

    data = deserialize_value(data)
    data = remove_units(data)

    # Convert seconds to hours
    times = list(data.keys())
    for time in times:
        hour = time / 3600
        data[hour] = data.pop(time)

    data_export = open(experiment_out_dir+'/data_export.pkl', 'wb')
    pickle.dump(data, data_export)
    data_export.close()

    sim_description = sim_config['description']

    config_export = open(experiment_out_dir + '/config_export.pkl', 'wb')
    pickle.dump(sim_description, config_export)
    config_export.close()

    print('saved '+str(data_export))

if __name__ == '__main__':
    data, experiment_id, sim_config = access()
    run_analysis(data, experiment_id, sim_config)
