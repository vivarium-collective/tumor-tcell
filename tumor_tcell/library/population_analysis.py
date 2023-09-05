from collections import Counter
from tumor_tcell.library.phylogeny import get_phylogeny
import pandas as pd


def division_analysis(cell_plot):
    #############################3
    # How to extract individual T cells
    df_divide_T = cell_plot.drop_duplicates('cell', keep='first')

    # Get unique agent IDs
    unique_T_cell = list(df_divide_T.cell.unique())

    # run phylogeny function
    phylogeny_T = get_phylogeny(unique_T_cell)

    # get initial ancestors, daughters, and mothers
    daughters_T = list(phylogeny_T.values())
    daughters_T = set([item for sublist in daughters_T for item in sublist])
    descendents_T = list(daughters_T)
    mothers_T = set(list(phylogeny_T.keys()))
    ancestors_T = list(mothers_T - daughters_T)

    # Time for plotting cell divisions
    div_list_T = []
    for cell in descendents_T:
        div = df_divide_T[df_divide_T["cell"] == cell]['time'].min()
        div_list_T.append(div)

    # get unique counts from the list
    div_counts_T = Counter(div_list_T)
    divide_time_T = pd.DataFrame.from_dict(div_counts_T, orient='index').reset_index()

    if not divide_time_T.empty:

        # convert to dataframe
        column_names = ['time', 'counts']
        divide_time_T.columns = column_names

        # divide counts by 2 because each daughter and original cell is counted twice
        divide_time_T['counts'] = divide_time_T['counts'] / 2

        # add 0, 0 initial point
        divide_time_T.loc[-1] = [0, 0]
        divide_time_T.index = divide_time_T.index + 1  # shifting index
        divide_time_T = divide_time_T.sort_values(by='time')

        # accumulate the counts as progresses
        divide_time_T['total_division'] = divide_time_T.counts.cumsum()

    else:
        divide_time_T = pd.DataFrame()

    return divide_time_T