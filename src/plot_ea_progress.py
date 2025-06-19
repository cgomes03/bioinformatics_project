import os
import sys
import argparse
import logging

from EA_rMTA.plotting.plotter import ResultPlotter

#-----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description='Generate EA result plots using ResultPlotter'
    )
    parser.add_argument(
        'run_folder',
        help='Path to EA run folder containing top.csv and final.csv'
    )
    parser.add_argument(
        '--top_n', '-n',
        type=int,
        default=None,
        help='If set, only use the top N candidates per generation'
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s'
    )

    run_folder = args.run_folder
    top_n = args.top_n

    # validate presence of top.csv
    top_csv = os.path.join(run_folder, 'top.csv')
    if not os.path.isfile(top_csv):
        print(f"ERROR: top.csv missing in {run_folder}", file=sys.stderr)
        sys.exit(1)

    # create plots
    plotter = ResultPlotter(run_folder)
    plotter.plot_all(top_n=top_n, save=True)
    print('Plots created: convergence_plot.png, boxplot.png')

if __name__ == '__main__':
    main()
