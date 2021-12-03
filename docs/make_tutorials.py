import subprocess
import glob
import os

# Convert the tutorials
for fn in glob.glob("../examples/*.ipynb"):
    name = os.path.splitext(os.path.split(fn)[1])[0]
    outfn = os.path.join("tutorials", name + ".rst")
    print("Building {0}...".format(name))
    subprocess.check_call(
        "jupyter nbconvert --template _templates/tutorial_rst.tpl --to rst "
        + fn
        + " --output-dir tutorials",
        shell=True,)