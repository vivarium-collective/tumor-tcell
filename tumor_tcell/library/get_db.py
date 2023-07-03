import vivarium
def get_data_fromdb(experiment_id=''):
    db=vivarium.core.emitter.get_experiment_database()
    data = vivarium.core.emitter.data_from_database(experiment_id=experiment_id, client=db)
    return data


if __name__ == '__main__':
    data = get_data_fromdb(experiment_id='tumor_tcell_20230623.211233')
    data