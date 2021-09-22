from radexwrap import *
from pandas import DataFrame, concat
import numpy as np
from functools import partial
import os

_ROOT = os.path.dirname(os.path.abspath(__file__))


def run(parameters, output_file=None):
    """
    Run a single RADEX model using a dictionary to set parameters.

    :param parameters: A dictionary containing the RADEX inputs that the user wishes to set,
        all other parameters will use the default values. See :func:`get_default_parameters`
        for a list of possible parameters and :func:`run_params` for descriptions.
    :type parameters: dict

    :param output_file: If not ``None``, the RADEX results are stored to this file in csv format/
    :type output_file: str
    """
    columns = ['E_UP (K)', 'freq', 'WAVEL (um)', 'T_ex', 'tau',
               'T_R (K)', 'POP UP', 'POP LOW', 'FLUX (K*km/s)', 'FLUX (erg/cm2/s)']



    parameters["molfile"] = add_data_path(parameters["molfile"])
    success,nlines, qup, qlow, output = from_dict(parameters)
    if success==1:
        output = DataFrame(columns=columns, data=output[:, :nlines].T)
        output["QN Upper"] = qup.reshape(-1, 6).view('S6')[:nlines]
        output["QN Lower"] = qlow.reshape(-1, 6).view('S6')[:nlines]
        output["Qup"] = output["QN Upper"].map(lambda x: x.decode('UTF-8')).str.strip()
        output["Qlow"] = output["QN Lower"].map(lambda x: x.decode('UTF-8')).str.strip()
        output=output.drop(["QN Upper","QN Lower"],axis=1)
        output=output[output["freq"]>parameters["fmin"]]
        output=output[output["freq"]<parameters["fmax"]]

        if output_file is not None:
            output.to_csv(output_file, index=False)
        return output
    else:
        print("RADEX Failed, check RADEX error messages\nYour parameters were:\n")
        print(parameters)
        return None

def run_params(molfile,tkin,cdmol,nh=0.0,nh2=0.0,op_ratio=3.0,ne=0.0,nhe=0.0,nhx=0.0,
        linewidth=1.0,fmin=0.0,fmax=500.0,tbg=2.73,geometry=1, output_file=None):
    """
    Run a single RADEX model from individual parameters

    :param molfile: Either the full path, starting from . or / to a datafile in the Lamda database
                    format or the filename of a datafile from `list_data_files()`.
    :type molfile: str

    :param tkin: Temperature of the Gas in Kelvin
    :type molfile: float

    :param cdmol: Column density of the emitting species in cm :math:`^{-2}`
    :type molfile: float

    :param nh: Number density of H atoms
    :type nh: float, optional

    :param nh2: Total number density of H2 molecules, set this to o-H2 + p-H2 if using ortho and para H2 as collisional partners.
    :type nh2: float, optional
    
    :param op_ratio: Ortho to para ratio for H2. Defaults to statistical limit of 3 and used to set o-H2 and p-H2 densities from nh2.
    :type op_ratio: float, optional

    :param ne: Number density of electron.
    :type ne: float, optional
    
    :param nhe: Number density of He atoms.
    :type nhe: float, optional 
    
    :param nhx: Number density of H+ ions.
    :type nh: float, optional

    :param linewidth: FWHM of the line in km s :math:`^{-1}`.
    :type linewidth: float, optional

    :param fmin: Minimum frequency below which a line is not included in the results.
    :type fmin: float, optional
    
    :param fmax: Maximum frequency above which a line is not included in the results.
    :type fmax: float, optional

    :param tbg: Background temperature, defaults to CMB temperature 2.73 K.
    :type tbg: float, optional  

    :param geometry: Choice of geometry of emitting object. 1 for sphere, 2 for LVG, 3 for slab.
    :type geometry: int, optional
    """
    columns = ['E_UP (K)', 'freq', 'WAVEL (um)', 'T_ex', 'tau',
               'T_R (K)', 'POP UP', 'POP LOW', 'FLUX (K*km/s)', 'FLUX (erg/cm2/s)']


    molfile = add_data_path(molfile)
    ortho=op_ratio/(op_ratio+1.0)
    para=1.0-ortho
    densities=[nh2,nh2*ortho,nh2*para,ne,nh,nhe,nhx]
    success,nlines, qup, qlow, output = from_params(molfile,tkin,tbg,cdmol,densities,
                                                    linewidth,fmin,fmax,geometry)
    if success==1:
        output = DataFrame(columns=columns, data=output[:, :nlines].T)
        output["QN Upper"] = qup.reshape(-1, 6).view('S6')[:nlines]
        output["QN Lower"] = qlow.reshape(-1, 6).view('S6')[:nlines]
        output["Qup"] = output["QN Upper"].map(lambda x: x.decode('UTF-8')).str.strip()
        output["Qlow"] = output["QN Lower"].map(lambda x: x.decode('UTF-8')).str.strip()
        output=output.drop(["QN Upper","QN Lower"],axis=1)
        output=output[output["freq"]>fmin]
        output=output[output["freq"]<fmax]
        if output_file is not None:
            output.to_csv(output_file, index=False)
        return output
    else:
        print("RADEX Failed, check RADEX error messages\nYour parameters were:\n")
        print(parameters)
        return None

def run_grid(parameters,
             target_value="FLUX (K*km/s)", pool=None):
    """
    Runs a grid of RADEX models using all combinations of any iterable items in the parameters dictionary whilst keeping other parameters constant. Returns a dataframe of results and can be parallelized with the ``pool`` parameter.

    :param parameters: A dictionary of parameters as provided by :func:`get_default_parameters` or :func:`get_example_grid_parameters`. Parameters should take a single value when they are constant over the grid and contain and interable if they are to be varied.

    :param molfile: Collisional data file for the molecule for which the emission should be calculated.
    :type molfile: str

    :param target_value: RADEX output column to be returned. Select one of 'T_R (K)', 'FLUX (K*km/s)', 'FLUX (erg/cm2/s)'
    :type target_value: str,optional

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
        "fmax": 1000.0,
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
        "fmax": 800.0,
        "geometry":1
    }
    return parameters

def get_transition_table(molfile):
    """
    Reads a collisional data file and returns a pandas DataFrame for the molecule with one row per transition containing the Einstein coefficients, upper level energy and frequency.

    :param molfile: Either the full path to a collisional datafile or the filename of one supplied with SpectralRadex
    :type molfile: str

    """
    molfile=add_data_path(molfile)
    with open(molfile) as f:
        f.readline()
        molecule=f.readline()
        f.readline()
        mass=f.readline()
        f.readline()
        levels=int(f.readline())
        f.readline()

        for i in range(levels):
            f.readline()
        f.readline()
        n_transitions=int(f.readline())
        f.readline()
        line_df=DataFrame(columns=["Upper level","Lower level","Aij","Frequency","E_u"])
        for i in range(n_transitions):
            line_df.loc[len(line_df)]=f.readline().split()[1:6]
    line_df[["Aij","Frequency","E_u"]]=line_df[["Aij","Frequency","E_u"]].astype(float)
    return line_df


def add_data_path(filename):
    #Adds the path to the packaged datafiles to a filename.
    if filename[0] not in ["~",".","/"]:
        return os.path.join(_ROOT, "data", filename)
    else:
        return filename


def list_data_files():
    """
    SpectralRadex is packaged with a selection of LAMDA collisional datafiles. 
    This function prints the list of available files. You can provide the full path to another
    file in the parameter dictionary to use one not packaged with SpectralRadex.
    """
    print(os.listdir(os.path.join(_ROOT, "data")))

def is_iter(x):
    return hasattr(x, '__iter__')

def format_run_for_grid(line_count, parameters, target_value, columns,grid_variables, grid_parameters):
    #Simple function to set up and reformat the output of :func:`run` for :func:`run_grid`
    for i,variable in enumerate(grid_variables):
        parameters[variable] = grid_parameters[i]
    radex_output = run(parameters)
    if radex_output is not None:
        transition_value = radex_output.iloc[:line_count][target_value].to_list()
    else:
        transition_value=[np.nan]*line_count
    return DataFrame([[parameters[x] for x in grid_variables] + transition_value], columns=columns)