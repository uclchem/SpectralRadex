from radexwrap import *
from pandas import DataFrame, concat
import numpy as np
from functools import partial
import os

_ROOT = os.path.dirname(os.path.abspath(__file__))


def run(parameters, output_file=None):
    """
    Run a single RADEX model

    :param parameters: A dictionary containing the RADEX inputs that the user wishes to set,
        all other parameters will use the default values. See :func:`get_default_parameters`
        for a list of possible parameters.
    :type parameters: dict

    :param output_file: If not ``None``, the RADEX results are stored to this file in csv format/
    :type output_file: str
    """
    columns = ['E_UP (K)', 'freq', 'WAVEL (um)', 'T_ex', 'tau',
               'T_R (K)', 'POP UP', 'POP LOW', 'FLUX (K*km/s)', 'FLUX (erg/cm2/s)']

    if parameters["molfile"][0] != "/":
        parameters["molfile"] = add_data_path(parameters["molfile"])

    nlines, qup, qlow, output = radex(parameters)
    output = DataFrame(columns=columns, data=output[:, :nlines].T)
    output["QN Upper"] = qup.reshape(-1, 6).view('S6')[:nlines]
    output["QN Lower"] = qlow.reshape(-1, 6).view('S6')[:nlines]
    output["Qup"] = output["QN Upper"].map(lambda x: x.decode('UTF-8')).str.strip()
    output["Qlow"] = output["QN Lower"].map(lambda x: x.decode('UTF-8')).str.strip()
    output=output.drop(["QN Upper","QN Lower"],axis=1)
    if output_file is not None:
        output.to_csv(output_file, index=False)
    return output


def format_run_for_grid(line_count, parameters, target_value, columns,grid_variables, grid_parameters):
    """
    Simple function to set up and reformat the output of :func:`run` for :func:`run_grid`
    :meta private:
    """
    for i,variable in enumerate(grid_variables):
        parameters[variable] = grid_parameters[i]
    radex_output = run(parameters)
    transition_value = radex_output.iloc[:line_count][target_value].to_list()
    return DataFrame([[parameters[x] for x in grid_variables] + transition_value], columns=columns)


def run_grid(parameters,
             target_value="FLUX (K*km/s)", freq_range=[0, 3.0e7], pool=None):
    """
    Runs a grid of RADEX models using all combinations of any iterable items in the parameters dictionary whilst keeping other parameters constant. Returns a dataframe of results and can be parallelized with the ``pool`` parameter.

    :param parameters: A dictionary of parameters as provided by :func:`get_default_parameters` or :func:`get_example_grid_parameters`. Parameters should take a single value when they are constant over the grid and contain and interable if they are to be varied.

    :param molfile: Collisional data file for the molecule for which the emission should be calculated.
    :type molfile: str

    :param target_value: RADEX output column to be returned. Select one of 'T_R (K)', 'FLUX (K*km/s)', 'FLUX (erg/cm2/s)'
    :type target_value: str,optional

    :param freq_range: Limit output lines to be those with frequencies in the range (fmin,fax).
    :type freq_range: iterable,float,optional

    :param pool: a Pool object with ``map()``, ``close()`` , and ``join()`` methods such as multiprocessing.Pool or schwimmbad.MPIPool.
         If supplied, the grid will be calculated in parallel. 
    :type pool: Pool, optional
    """

    #cleaning up and checking inputs
    if type(target_value) != str:
        if (target_value != "T_R (K)" or target_value != "FLUX (K*km/s)" or target_value != "FLUX (erg/cm2/s)"):
            print("target_value must be a string and one of the following three options: 'T_R (K)', 'FLUX (K*km/s)', 'FLUX (erg/cm2/s)'")
            exit()
    if parameters["molfile"][0] != "/":
        parameters["molfile"] = add_data_path(parameters["molfile"])

    #get list of varying parameters
    variables=[key for key,value in parameters.items() if is_iter(value) and key!="molfile"]

    parameter_grid = np.array(np.meshgrid(*[parameters[x] for x in variables]))
    parameter_grid = parameter_grid.T.reshape(-1, len(variables))
    parameters=parameters.copy()
    for i,variable in enumerate(variables):
        parameters[variable]=parameter_grid[0,i]
    parameter_combinations = np.delete(parameter_grid, 0, axis=0)
    radex_output = run(parameters)

    dataframe_columns = variables[:]
    line_count = np.shape(radex_output)[0]
    transition_value = []
    for line in range(line_count):
        transition_name = "(" + radex_output.iloc[line]['Qup'] + ")-(" + \
                          radex_output.iloc[line]['Qlow'] + ")[" + \
                          str(radex_output.iloc[line]['freq']) + " GHz]"
        dataframe_columns += [transition_name]
        transition_value += [radex_output.iloc[line][target_value]]
    output = DataFrame(columns=dataframe_columns, data=[[parameters[x] for x in variables] + transition_value])

    if pool is not None:
        func = partial(format_run_for_grid, line_count, parameters, target_value, dataframe_columns,variables)
        pool_results = pool.map(func, parameter_combinations)
        pool.close()
        pool.join()
        pool_results_df = concat(pool_results, axis=0)
        output = output.append(pool_results_df, ignore_index=True)
    else:
        for grid_point in range(len(parameter_combinations)):
            output = output.append(format_run_for_grid(line_count, parameters,
                                            target_value, dataframe_columns,variables,
                                            parameter_combinations[grid_point]), ignore_index=True)

    return output


def get_default_parameters():
    """
    Get the default RADEX parameters as a dictionary, this largely serves as an example for the
    input required for :func:`run`.

    molfile can be a complete path to a collisional datafile in the LAMDA database format or one of the files
    listed by :func:`list_data_files`.

    method is 1 (uniform sphere), 2 (LVG), or 3 (slab)
    """
    parameters = {
        "molfile": "co.dat",
        "tkin": 30.0,
        "tbg": 2.73,
        "cdmol": 1.0e13,
        "h2": 1.0e5,
        "h": 0.0,
        "e-": 0.0,
        "p-h2": 0.0,
        "o-h2": 0.0,
        "h+": 0.0,
        "linewidth": 1.0,
        "fmin": 0.0,
        "fmax": 3.0e7,
        "geometry":1
    }
    return parameters

def get_example_grid_parameters():
    """
    Returns a dictionary of parameters for RADEX with iterables which can be used with :func:`run_grid`.
    """

    parameters = {
        "molfile": "co.dat",
        "tkin": np.linspace(10,300,5),
        "tbg": 2.73,
        "cdmol": np.logspace(14,18,5),
        "h2": np.logspace(4,7,5),
        "h": 0.0,
        "e-": 0.0,
        "p-h2": 0.0,
        "o-h2": 0.0,
        "h+": 0.0,
        "linewidth": 1.0,
        "fmin": 0.0,
        "fmax": 3.0e7,
        "geometry":1
    }
    return parameters


def add_data_path(filename):
    """
    Adds the path to the packaged datafiles to a filename.
    :meta private:
    """
    return os.path.join(_ROOT, "data", filename)


def list_data_files():
    """
    SpectralRadex is packaged with a selection of LAMDA collisional datafiles. 
    This function prints the list of available files. You can provide the full path to another
    file in the parameter dictionary to use one not packaged with SpectralRadex.
    """
    print(os.listdir(os.path.join(_ROOT, "data")))

def is_iter(x):
    return hasattr(x, '__iter__')