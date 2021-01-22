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


def format_run_for_grid(line_count, parameters, target_value, columns, grid_parameters):
    """
    Simple function to set up and reformat the output of :func:`run` for :func:`run_grid`
    """
    parameters["h2"] = grid_parameters[0]
    parameters["tkin"] = grid_parameters[1]
    parameters["cdmol"] = grid_parameters[2]
    radex_output = run(parameters)
    transition_value = radex_output.iloc[:line_count][target_value].to_list()
    return DataFrame([[parameters["h2"], parameters["tkin"], parameters["cdmol"]] + transition_value], columns=columns)


def run_grid(density_values, temperature_values, column_density_values, molfile,
             target_value="FLUX (K*km/s)", freq_range=[0, 3.0e7], pool=None):
    """
    Runs a grid of RADEX models using all combinations of input lists of H2 density, temperature, and the molecular column 
    density. Returns a dataframe of results and can be parallelized with the ``pool`` parameter.

    :param density_values: A list of density values for which RADEX models should be run.
    :type density_values: iterable

    :param temperature_values: A list of temperature values for which RADEX models should be run.
    :type temperature_values: iterable

    :param column_density_values: A list of column density values for which RADEX models should be run.
    :type column_density_values: iterable

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

    if type(target_value) != str:
        if (target_value != "T_R (K)" or target_value != "FLUX (K*km/s)" or target_value != "FLUX (erg/cm2/s)"):
            print("target_value must be a string and one of the following three options: 'T_R (K)', 'FLUX (K*km/s)', 'FLUX (erg/cm2/s)'")
            exit()
    if type(molfile) != str:
        print("molfile must be a string of the moleculer data file to be used.")
        exit()
    parameters = get_default_parameters()
    parameters["fmin"] = freq_range[0]
    parameters["fmax"] = freq_range[1]
    parameters["molfile"] = molfile
    parameter_grid = np.array(np.meshgrid(density_values, temperature_values, column_density_values))
    parameter_combinations = parameter_grid.T.reshape(-1, 3)
    del parameter_grid

    parameters["h2"] = parameter_combinations[0, 0]
    parameters["tkin"] = parameter_combinations[0, 1]
    parameters["cdmol"] = parameter_combinations[0, 2]
    parameter_combinations = np.delete(parameter_combinations, 0, axis=0)
    radex_output = run(parameters)

    dataframe_columns = ["Density", "Temperature", "Column Density"]
    line_count = np.shape(radex_output)[0]
    transition_value = []
    for line in range(line_count):
        transition_name = "(" + radex_output.iloc[line]['Qup'] + ")-(" + \
                          radex_output.iloc[line]['Qlow'] + ")[" + \
                          str(radex_output.iloc[line]['freq']) + " GHz]"
        dataframe_columns += [transition_name]
        transition_value += [radex_output.iloc[line][target_value]]
    output = DataFrame(columns=dataframe_columns, data=[[parameters["h2"], parameters["tkin"],
                                                         parameters["cdmol"]] + transition_value])

    if pool is not None:
        func = partial(format_run_for_grid, line_count, parameters, target_value, dataframe_columns)
        pool_results = pool.map(func, parameter_combinations)
        pool.close()
        pool.join()
        pool_results_df = concat(pool_results, axis=0)
        print(np.shape(pool_results_df))
        output = output.append(pool_results_df, ignore_index=True)
    else:
        for grid_point in range(len(parameter_combinations)):
            output = output.append(format_run_for_grid(line_count, parameters, target_value, dataframe_columns,
                                             parameter_combinations[grid_point]), ignore_index=True)

    return output


def get_default_parameters():
    """
    Get the default RADEX parameters as a dictionary, this largely serves as an example for the
    input required for :func:`run`.
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
        "fmax": 3.0e7
    }
    return parameters


def add_data_path(filename):
    return os.path.join(_ROOT, "data", filename)


def list_data_files():
    """
    Prints the lambda datafiles packaged with spectralradex
    """
    print(os.listdir(os.path.join(_ROOT, "data")))
