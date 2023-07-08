import vivarium
import os
import pickle
from tumor_tcell.library.individual_analysis import individual_analysis
from vivarium.library.units import remove_units
from vivarium.core.emitter import deserialize_value, DatabaseEmitter


def get_data_fromdb(experiment_id='', save_dir=None, analyze=True, bounds=[1200, 1200]):
    db = vivarium.core.emitter.get_experiment_database()
    data = vivarium.core.emitter.data_from_database(experiment_id=experiment_id, client=db)
    # data = deserialize_value(data)
    # emitter = DatabaseEmitter(dict(experiment_id=experiment_id))
    # emitter.client = db
    # data = emitter.get_data()
    data_export = data[0]
    config_export = data[1]

    # Specify the directory path to save the file
    if save_dir is None:
        home_dir = "/home/nolanlab/Documents/GitHub/tumor-tcell/tumor_tcell/experiments/analysis/"
        save_dir = home_dir + experiment_id

        # Create the directory if it doesn't exist
        os.makedirs(save_dir, exist_ok=True)

    # Save the file with the full path
    data_path = os.path.join(save_dir, "data_export.pkl")
    with open(data_path, "wb") as file:
        pickle.dump(data_export, file)

    file_path = os.path.join(save_dir, "config_export.pkl")
    with open(file_path, "wb") as file:
        pickle.dump(config_export, file)

    if analyze==True:
        individual_analysis(analysis_dir=home_dir, experiment_id=experiment_id, bounds=bounds, tcells=True, lymph_nodes=True)

    return data


if __name__ == '__main__':
    data = get_data_fromdb(experiment_id='tumor_tcell_20230626.211422', analyze=True)
    data