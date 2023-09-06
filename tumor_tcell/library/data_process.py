import pandas as pd

def data_to_dataframes(data, lymph_nodes=False):

    # Create a new dictionary with the 'tumor_environment' data based on restructuring
    new_data = {}

    for key, value in data.items():
        # Extract 'tumor_environment' data and keep the outer structure
        new_data[key] = value['tumor_environment']
    new_data.keys()

    #Convert to initial dataframe
    df_data = pd.DataFrame(new_data)

    #Extract agents from the data
    df_copy = df_data.copy()
    df_agents = df_copy.loc['agents']
    agents_dict = df_agents.to_dict()

    #reformat the dictionary into mulitiindexed dataframe
    reform_agents = {(level1_key, level2_key): values
              for level1_key, level2_dict in agents_dict.items()
              for level2_key, values in level2_dict.items()
              }
    df_agents_exp1 = pd.DataFrame(reform_agents)
    df_agents_multi = df_agents_exp1.T
    names=['time', 'cell']
    df_agents_multi.index.set_names(names, inplace=True)



    ###################################
    #Subset only T cells from all agents
    df_tcell_agents = df_agents_multi.iloc[df_agents_multi.index.get_level_values('cell').str.contains('tcell'), :]

    #Subset categories and recombine in T cells
    df_tcell_trans = df_tcell_agents.T
    tcell_categories = []

    #Extract each feature - boundary, internal, neighbors
    for category in df_tcell_trans.index.values:
      df_boundary_sub = df_tcell_trans.loc[category,:]
      tcell_boundary_dict = df_boundary_sub.to_dict()
      df_boundary_sub2 = pd.DataFrame(tcell_boundary_dict)
      tcell_categories.append(df_boundary_sub2.T)

    #concatenate dataframes
    tcell_data = pd.concat(tcell_categories, axis=1)

    #reformat T cell data for plotting
    tcell_data['IFNg'] = tcell_data['external'].apply(lambda x: x.get('IFNg'))
    tcell_data['transferable_cytotoxic_packets'] = tcell_data['transfer'].apply(lambda x: x.get('cytotoxic_packets'))
    tcell_data['X'] = tcell_data['location'].apply(lambda x: x[0])
    tcell_data['Y'] = tcell_data['location'].apply(lambda x: x[1])

    #Only select columns of interest
    tcell_columns = ['cell_state', 'IFNg', 'transferable_cytotoxic_packets', 'X', 'Y']
    tcell_data_form = tcell_data[tcell_columns]
    tcell_data_form.index.set_names(names, inplace=True)

    # reset index for plotting
    tcell_plot = tcell_data_form.reset_index()


    ########################################3
    # Subset only Tumor cells from all agents
    df_tumor_agents = df_agents_multi.iloc[df_agents_multi.index.get_level_values('cell').str.contains('tumor'), :]

    # Subset categories and recombine in Tumor cells
    df_tumor_trans = df_tumor_agents.T
    tumor_categories = []

    # Extract each feature - boundary, internal, neighbors
    for category in df_tumor_trans.index.values:
        df_boundary_sub = df_tumor_trans.loc[category, :]
        tumor_boundary_dict = df_boundary_sub.to_dict()
        df_boundary_sub2 = pd.DataFrame(tumor_boundary_dict)
        tumor_categories.append(df_boundary_sub2.T)

    # concatenate dataframes
    tumor_data = pd.concat(tumor_categories, axis=1)
    tumor_data;

    # reformat Tumor cell data for plotting
    tumor_data['IFNg'] = tumor_data['external'].apply(lambda x: x.get('IFNg'))
    tumor_data['cytotoxic_packets'] = tumor_data['receive'].apply(lambda x: x.get('cytotoxic_packets'))
    tumor_data['X'] = tumor_data['location'].apply(lambda x: x[0])
    tumor_data['Y'] = tumor_data['location'].apply(lambda x: x[1])
    tumor_data.columns

    # Only select columns of interest
    tumor_columns = ['cell_state', 'IFNg', 'cytotoxic_packets', 'X', 'Y']
    tumor_data_form = tumor_data[tumor_columns]
    tumor_data_form.index.set_names(names, inplace=True)

    # reset index for plotting
    tumor_plot = tumor_data_form.reset_index()


    ################################
    ####Extract death log statistics
    df_death = df_copy.loc['log']
    death_dict = df_death.to_dict()

    # reformat the dictionary into mulitiindexed dataframe
    reform_death = {(level1_key, level2_key): values
                    for level1_key, level2_dict in death_dict.items()
                    for level2_key, values in level2_dict.items()}

    # Make dataframe
    df_death_exp1 = pd.DataFrame(reform_death)
    df_death_multi = df_death_exp1.T
    names = ['time', 'cell']
    df_death_multi.index.set_names(names, inplace=True)
    df_death_multi.columns = ['time', 'death']

    # subset only where death is not equal to false
    df_death_sub = df_death_multi[~(df_death_multi['death']==False)]

    if df_death_sub.empty:
        df_tcell_death = pd.DataFrame({})
        df_tumor_death = pd.DataFrame({})
    else:
        # Only get the final log of the death than contains all the death information
        df_last_death = df_death_sub.loc[df_death_sub.index.levels[0][-1]]
        df_last_death['time'] = df_last_death['time'] / 3600

        # Subset only T cells from all agents
        df_tcell_death = df_last_death.iloc[df_last_death.index.get_level_values('cell').str.contains('tcell'), :]
        df_tumor_death = df_last_death.iloc[df_last_death.index.get_level_values('cell').str.contains('tumor'), :]

        ########################################
        ##Do for T cells
        # sort deaths by time
        df_tcell_death.sort_values(by=['time'], inplace=True)

        # Get different death groupings and count total over time
        death_types_list = list(df_tcell_death['death'].unique())
        for death_type in death_types_list:
            df_tcell_death[death_type] = df_tcell_death['death'].apply(lambda x: 1 if x == death_type else 0)
            df_tcell_death['total_' + str(death_type)] = df_tcell_death[death_type].cumsum()

        # get total death count over time
        total_col_t = [col for col in df_tcell_death.columns if 'total' in col]
        df_tcell_death['total_death'] = df_tcell_death[total_col_t].sum(axis=1)

        ##Do for Tumors
        # sort deaths by time
        df_tumor_death.sort_values(by=['time'], inplace=True)

        # Get different death groupings and count total over time
        death_types_list = list(df_tumor_death['death'].unique())
        for death_type in death_types_list:
            df_tumor_death[death_type] = df_tumor_death['death'].apply(lambda x: 1 if x == death_type else 0)
            df_tumor_death['total_' + str(death_type)] = df_tumor_death[death_type].cumsum()

        # get total death count over time
        total_col = [col for col in df_tumor_death.columns if 'total' in col]
        df_tumor_death['total_death'] = df_tumor_death[total_col].sum(axis=1)

    if lymph_nodes==True:
        ###################################
        # Subset only DC cells from all agents
        df_dendritic_agents = df_agents_multi.iloc[
                              df_agents_multi.index.get_level_values('cell').str.contains('dendritic'), :]

        # Subset categories and recombine in T cells
        df_dendritic_trans = df_dendritic_agents.T
        dendritic_categories = []

        # Extract each feature - boundary, internal, neighbors
        for category in df_dendritic_trans.index.values:
            df_boundary_sub = df_dendritic_trans.loc[category, :]
            dendritic_boundary_dict = df_boundary_sub.to_dict()
            df_boundary_sub2 = pd.DataFrame(dendritic_boundary_dict)
            dendritic_categories.append(df_boundary_sub2.T)

        # concatenate dataframes
        dendritic_data = pd.concat(dendritic_categories, axis=1)

        # reformat T cell data for plotting
        dendritic_data['tumor_debris'] = dendritic_data['external'].apply(lambda x: x.get('tumor_debris'))
        dendritic_data['X'] = dendritic_data['location'].apply(lambda x: x[0])
        dendritic_data['Y'] = dendritic_data['location'].apply(lambda x: x[1])

        # Only select columns of interest
        dendritic_columns = ['cell_state', 'tumor_debris', 'X', 'Y']
        dendritic_data_form = dendritic_data[dendritic_columns]
        dendritic_data_form.index.set_names(names, inplace=True)

        # reset index for plotting
        dendritic_plot = dendritic_data_form.reset_index()

        if df_death_sub.empty:
            df_dendritic_death = pd.DataFrame({})
        else:
            #get dendritic death stats
            df_dendritic_death = df_last_death.iloc[df_last_death.index.get_level_values('cell').str.contains('dendritic'),:]

            ########################################
            ##Do for dendritic cells
            # sort deaths by time
            df_dendritic_death.sort_values(by=['time'], inplace=True)

            # Get different death groupings and count total over time
            death_types_list = list(df_dendritic_death['death'].unique())
            for death_type in death_types_list:
                df_dendritic_death[death_type] = df_dendritic_death['death'].apply(lambda x: 1 if x == death_type else 0)
                df_dendritic_death['total_' + str(death_type)] = df_dendritic_death[death_type].cumsum()

            # get total death count over time
            total_col_t = [col for col in df_dendritic_death.columns if 'total' in col]
            df_dendritic_death['total_death'] = df_dendritic_death[total_col_t].sum(axis=1)

        return df_tumor_death, df_tcell_death, tumor_plot, tcell_plot, df_dendritic_death, dendritic_plot
    else:
        return df_tumor_death, df_tcell_death, tumor_plot, tcell_plot



#Function for processing data when there are no T cells

def control_data_to_dataframes(data):
    # Convert to initial dataframe
    df_data = pd.DataFrame(data)

    # Extract agents from the data
    df_copy = df_data.copy()
    df_agents = df_copy.iloc[0, :]
    agents_dict = df_agents.to_dict()

    # reformat the dictionary into mulitiindexed dataframe
    reform_agents = {(level1_key, level2_key): values
                     for level1_key, level2_dict in agents_dict.items()
                     for level2_key, values in level2_dict.items()
                     }
    df_agents_exp1 = pd.DataFrame(reform_agents)
    df_agents_multi = df_agents_exp1.T
    names = ['time', 'cell']
    df_agents_multi.index.set_names(names, inplace=True)

    ########################################3
    # Subset only Tumor cells from all agents
    df_tumor_agents = df_agents_multi.iloc[df_agents_multi.index.get_level_values('cell').str.contains('tumor'), :]

    # Subset categories and recombine in Tumor cells
    df_tumor_trans = df_tumor_agents.T
    tumor_categories = []

    # Extract each feature - boundary, internal, neighbors
    for category in df_tumor_trans.index.values:
        df_boundary_sub = df_tumor_trans.loc[category, :]
        tumor_boundary_dict = df_boundary_sub.to_dict()
        df_boundary_sub2 = pd.DataFrame(tumor_boundary_dict)
        tumor_categories.append(df_boundary_sub2.T)

    # concatenate dataframes
    tumor_data = pd.concat(tumor_categories, axis=1)
    tumor_data;

    # reformat Tumor cell data for plotting
    tumor_data['IFNg'] = tumor_data['external'].apply(lambda x: x.get('IFNg'))
    tumor_data['cytotoxic_packets'] = tumor_data['receive'].apply(lambda x: x.get('cytotoxic_packets'))
    tumor_data['X'] = tumor_data['location'].apply(lambda x: x[0])
    tumor_data['Y'] = tumor_data['location'].apply(lambda x: x[1])
    tumor_data.columns

    # Only select columns of interest
    tumor_columns = ['cell_state', 'IFNg', 'cytotoxic_packets', 'X', 'Y']
    tumor_data_form = tumor_data[tumor_columns]
    tumor_data_form.index.set_names(names, inplace=True)

    # reset index for plotting
    tumor_plot = tumor_data_form.reset_index()

    ################################
    ####Extract death log statistics
    df_death = df_copy.iloc[1, :]
    death_dict = df_death.to_dict()

    # reformat the dictionary into mulitiindexed dataframe
    reform_death = {(level1_key, level2_key): values
                    for level1_key, level2_dict in death_dict.items()
                    for level2_key, values in level2_dict.items()}

    # Make dataframe
    df_death_exp1 = pd.DataFrame(reform_death)
    df_death_multi = df_death_exp1.T
    names = ['time', 'cell']
    df_death_multi.index.set_names(names, inplace=True)
    df_death_multi.columns = ['time', 'death']

    # subset only where death is not equal to false
    df_death_sub = df_death_multi[~(df_death_multi['death']==False)]

    if df_death_sub.empty:
        return pd.DataFrame({}), tumor_plot

    # Only get the final log of the death than contains all the death information
    df_last_death = df_death_sub.loc[df_death_sub.index.levels[0][-1]]
    df_last_death['time'] = df_last_death['time'] / 3600

    # Subset only T cells from all agents
    df_tumor_death = df_last_death.iloc[df_last_death.index.get_level_values('cell').str.contains('tumor'), :]

    ##Do for Tumors
    # sort deaths by time
    df_tumor_death.sort_values(by=['time'], inplace=True)

    # Get different death groupings and count total over time
    death_types_list = list(df_tumor_death['death'].unique())
    for death_type in death_types_list:
        df_tumor_death[death_type] = df_tumor_death['death'].apply(lambda x: 1 if x == death_type else 0)
        df_tumor_death['total_' + str(death_type)] = df_tumor_death[death_type].cumsum()

    # get total death count over time
    total_col = [col for col in df_tumor_death.columns if 'total' in col]
    df_tumor_death['total_death'] = df_tumor_death[total_col].sum(axis=1)

    return df_tumor_death, tumor_plot