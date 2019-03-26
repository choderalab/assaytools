import pymc
from pymc import MCMC
import assaytools
from pymc.Matplot import plot
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

def entry_point():

    import sys

    #This allows us to import local inputs.py
    sys.path.append('./')

    import argparse

    # Define argparse stuff

    parser = argparse.ArgumentParser(description="""Plot the trace from the pickled object written by quickmodel """)
    parser.add_argument("--input", help="The pickle object you'd like to plot",required=True)
    parser.add_argument("--name", help="Name and extension of file that will be saved",required=True)
    args = parser.parse_args()

    # Define inputs
    # If an inputs###.py is defined, it is run and used. An inputs.json is also created.
    inputs_file = args.input
    filename = args.name

    db = pymc.database.pickle.load(inputs_file)
    ntraces = len(db.trace_names[0])

    sns.set_style('white')
    sns.set_palette('muted', 12)
    fig = plt.figure(figsize=(6, 4 * ntraces))
    for i, trace_name in enumerate(db.trace_names[0]):
        plt.subplot(ntraces, 1, i + 1)
        plt.plot(range(len(db.trace(trace_name)[:])), db.trace(trace_name)[:])
        plt.title(trace_name)
    fig.savefig(filename, dpi=500, bbox_inches='tight')


if __name__ == '__main__':
    entry_point()
