from spectralradex import radex

params={'molfile': 'SO-pH2.dat', 'tkin': 240.23288848051274, 'tbg': 2.73, 'cdmol': 2.7117385194476626e+19, 'h2': 1477400.1189838066, 'h': 0.0, 'e-': 0.0, 'p-h2': 369350.02974595164, 'o-h2': 1108050.0892378548, 'h+': 0.0, 'linewidth': 125.49076959987242, 'fmin': 0.0, 'fmax': 30000000.0}
result=radex.run(params)
print("hi")
print(result)