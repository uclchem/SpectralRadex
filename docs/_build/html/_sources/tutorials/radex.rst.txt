**Note:** This tutorial was generated from an IPython notebook that can be
downloaded `here <https://github.com/uclchem/SpectralRadex/tree/master/examples>`_.

.. _radex:

Radex
=====

.. code:: python

    from spectralradex import radex
    from multiprocessing import Pool
    import numpy as np
    import time


The simplest use case for SpectralRadex is to be a simple python wrapper
for RADEX. This allows large grids of RADEX models or complex parameter
inference procedures to be run in an environment suited to those tasks.

If one wishes to run radex, we simply need a dictionary of the
parameters RADEX expects. An example can be obtained using the
``get_default_parameters()`` function like so

.. code:: python

    params = radex.get_default_parameters()
    print("{")
    for key,value in params.items():
        print(f"\t{key} : {value}")
    print("}")


.. parsed-literal::

    {
    	molfile : co.dat
    	tkin : 30.0
    	tbg : 2.73
    	cdmol : 10000000000000.0
    	h2 : 100000.0
    	h : 0.0
    	e- : 0.0
    	p-h2 : 0.0
    	o-h2 : 0.0
    	h+ : 0.0
    	linewidth : 1.0
    	fmin : 0.0
    	fmax : 30000000.0
    	geometry : 1
    }


and then we pass that to the ``run()`` function.

.. code:: python

    output = radex.run(params)
    output.head()




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>E_UP (K)</th>
          <th>freq</th>
          <th>WAVEL (um)</th>
          <th>T_ex</th>
          <th>tau</th>
          <th>T_R (K)</th>
          <th>POP UP</th>
          <th>POP LOW</th>
          <th>FLUX (K*km/s)</th>
          <th>FLUX (erg/cm2/s)</th>
          <th>Qup</th>
          <th>Qlow</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>5.53</td>
          <td>115.271202</td>
          <td>2600.757633</td>
          <td>31.666252</td>
          <td>0.000223</td>
          <td>0.006275</td>
          <td>0.246666</td>
          <td>0.097917</td>
          <td>0.006680</td>
          <td>1.317591e-10</td>
          <td>1</td>
          <td>0</td>
        </tr>
        <tr>
          <th>1</th>
          <td>16.60</td>
          <td>230.538000</td>
          <td>1300.403656</td>
          <td>29.262261</td>
          <td>0.000735</td>
          <td>0.017551</td>
          <td>0.281677</td>
          <td>0.246666</td>
          <td>0.018683</td>
          <td>2.947981e-09</td>
          <td>2</td>
          <td>1</td>
        </tr>
        <tr>
          <th>2</th>
          <td>33.19</td>
          <td>345.795990</td>
          <td>866.963374</td>
          <td>26.640080</td>
          <td>0.001112</td>
          <td>0.021294</td>
          <td>0.211510</td>
          <td>0.281677</td>
          <td>0.022667</td>
          <td>1.207049e-08</td>
          <td>3</td>
          <td>2</td>
        </tr>
        <tr>
          <th>3</th>
          <td>55.32</td>
          <td>461.040768</td>
          <td>650.251515</td>
          <td>24.363876</td>
          <td>0.001022</td>
          <td>0.015261</td>
          <td>0.109663</td>
          <td>0.211510</td>
          <td>0.016246</td>
          <td>2.050309e-08</td>
          <td>4</td>
          <td>3</td>
        </tr>
        <tr>
          <th>4</th>
          <td>82.97</td>
          <td>576.267931</td>
          <td>520.231028</td>
          <td>22.798547</td>
          <td>0.000605</td>
          <td>0.007078</td>
          <td>0.039845</td>
          <td>0.109663</td>
          <td>0.007535</td>
          <td>1.856956e-08</td>
          <td>5</td>
          <td>4</td>
        </tr>
      </tbody>
    </table>
    </div>



Parameter Grids
---------------

It is more likely that one will want to run the code over many
combinations of input parameters. This can be achieved via the
``run_grid()`` function. This function also takes a parameter dictionary
of the same format as ``run()``. However, variables which are too be
varied over the grid should be supplied as iterables.

Furthermore, to keep things simple, the desired RADEXtakes iterables for
the three variables (density, temperature and column density) as well as
fixed values for the other RADEX parameters. It then produces the RADEX
output for all combinations of the three iterables.

We’ll use an example grid which can be acquired using the
``get_example_grid_parameters()`` function.

.. code:: python

    parameters=radex.get_example_grid_parameters()
    parameters




.. parsed-literal::

    {'molfile': 'co.dat',
     'tkin': array([ 10. ,  82.5, 155. , 227.5, 300. ]),
     'tbg': 2.73,
     'cdmol': array([1.e+14, 1.e+15, 1.e+16, 1.e+17, 1.e+18]),
     'h2': array([   10000.        ,    56234.13251903,   316227.76601684,
             1778279.41003892, 10000000.        ]),
     'h': 0.0,
     'e-': 0.0,
     'p-h2': 0.0,
     'o-h2': 0.0,
     'h+': 0.0,
     'linewidth': 1.0,
     'fmin': 0.0,
     'fmax': 30000000.0,
     'geometry': 1}



.. code:: python

    tic = time.perf_counter()
    
    grid_df = radex.run_grid(parameters,target_value="T_R (K)")
    toc = time.perf_counter()
    print(f"run_grid took {toc-tic:0.4f} seconds without a pool")


.. parsed-literal::

    run_grid took 3.0571 seconds without a pool


.. code:: python

    grid_df.iloc[:,0:6].head()




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>tkin</th>
          <th>cdmol</th>
          <th>h2</th>
          <th>(1)-(0)[115.2712018 GHz]</th>
          <th>(2)-(1)[230.538 GHz]</th>
          <th>(3)-(2)[345.7959899 GHz]</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>10.0</td>
          <td>1.000000e+14</td>
          <td>10000.0</td>
          <td>0.114622</td>
          <td>0.108152</td>
          <td>0.022018</td>
        </tr>
        <tr>
          <th>1</th>
          <td>10.0</td>
          <td>1.000000e+15</td>
          <td>10000.0</td>
          <td>1.048925</td>
          <td>0.958338</td>
          <td>0.215099</td>
        </tr>
        <tr>
          <th>2</th>
          <td>10.0</td>
          <td>1.000000e+16</td>
          <td>10000.0</td>
          <td>5.189712</td>
          <td>4.045272</td>
          <td>1.567682</td>
        </tr>
        <tr>
          <th>3</th>
          <td>10.0</td>
          <td>1.000000e+17</td>
          <td>10000.0</td>
          <td>6.561081</td>
          <td>5.156221</td>
          <td>3.411413</td>
        </tr>
        <tr>
          <th>4</th>
          <td>10.0</td>
          <td>1.000000e+18</td>
          <td>10000.0</td>
          <td>6.639451</td>
          <td>5.259944</td>
          <td>3.822848</td>
        </tr>
      </tbody>
    </table>
    </div>



Parallelization
~~~~~~~~~~~~~~~

In order to be as flexible as possible, SpectralRadex has no built in
multiprocessing. However, the ``run_grid()`` function does take the
optional parameter ``pool`` which should be an object with ``map()``,
``join()``, and ``close()`` methods that allow functions to be evaluated
in parallel. For example, the python standard
`multiprocessing.pool <https://docs.python.org/3.6/library/multiprocessing.html>`__
obect or Schwimmbad’s
`MPIPool <https://schwimmbad.readthedocs.io/en/latest/examples/#using-mpipool>`__.

If such an object is supplied, the grid will be evaluated in parallel.
Note the time in the example below compared to the grid above.

.. code:: python

    tic = time.perf_counter()
    pool=Pool(8)
    grid_df = radex.run_grid(parameters,target_value="T_R (K)",pool=pool)
    toc = time.perf_counter()
    print(f"run_grid took {toc-tic:0.4f} seconds with a pool of 8 workers")
    grid_df.iloc[:,0:6].head()


.. parsed-literal::

    run_grid took 0.7459 seconds with a pool of 8 workers




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>tkin</th>
          <th>cdmol</th>
          <th>h2</th>
          <th>(1)-(0)[115.2712018 GHz]</th>
          <th>(2)-(1)[230.538 GHz]</th>
          <th>(3)-(2)[345.7959899 GHz]</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>10.0</td>
          <td>1.000000e+14</td>
          <td>10000.0</td>
          <td>0.114622</td>
          <td>0.108152</td>
          <td>0.022018</td>
        </tr>
        <tr>
          <th>1</th>
          <td>10.0</td>
          <td>1.000000e+15</td>
          <td>10000.0</td>
          <td>1.048925</td>
          <td>0.958338</td>
          <td>0.215099</td>
        </tr>
        <tr>
          <th>2</th>
          <td>10.0</td>
          <td>1.000000e+16</td>
          <td>10000.0</td>
          <td>5.189712</td>
          <td>4.045272</td>
          <td>1.567682</td>
        </tr>
        <tr>
          <th>3</th>
          <td>10.0</td>
          <td>1.000000e+17</td>
          <td>10000.0</td>
          <td>6.561081</td>
          <td>5.156221</td>
          <td>3.411413</td>
        </tr>
        <tr>
          <th>4</th>
          <td>10.0</td>
          <td>1.000000e+18</td>
          <td>10000.0</td>
          <td>6.639451</td>
          <td>5.259944</td>
          <td>3.822848</td>
        </tr>
      </tbody>
    </table>
    </div>


