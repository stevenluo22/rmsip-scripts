import numpy as np
import argparse
import importlib.util
import matplotlib.pyplot as plt
import csv

def readcsv(csvfile):
    rmsip_values = {}
    with open(csvfile, mode='r') as file:
        reader = csv.reader(file)
        for row in reader:
            row_key = int(row[0])  # First column as key
            rmsip_values[row_key] = {i: float(value) for i, value in enumerate(row[1:], start=1)}
    return rmsip_values

def run(args):
    plot_values = {}
    if args.csvinput == "Use pythonlist input":
        spec = importlib.util.spec_from_file_location("csvlist", args.csvsinput)
        CSVlist = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(CSVlist)
        csvlist = CSVlist.csvList()
        
        for csvname, csv in csvlist.items():
            data = readcsv(csv)
            for run_name, run_data in data.items():
                plot_values[f"{csvname}_{run_name}"] = run_data
    else: 
        plot_values[args.name] = readcsv(args.csvinput)
    

    time = 0
    times = []
    for i in range(len(run_data)):
        time += args.time/len(run_data)
        times.append(time)

    plt.figure(figsize=(10, 6))

    # Plot each run
    r = 0
    for name, rmsip_data in plot_values.items():
        r += 1
        plt.plot(times, rmsip_data, label=f"{name}", marker='o')

    # Customize the plot
    plt.title(args.figureTitle)
    plt.xlabel('Time')
    plt.ylabel('RMSIP')
    plt.xticks(rotation=90)  # Rotate x-axis labels for better readability
    plt.legend(title='Runs', bbox_to_anchor=(1, 1), loc='upper left')
    plt.tight_layout()

    plt.savefig(args.outputplot, dpi=300, bbox_inches='tight')

    if args.addPrint:
        print(f"Plot saved as {args.outputplot}")    

def main(args=None):
    parser = argparse.ArgumentParser(
        description="Generate or modify a plot of RMSIPs")
    parser.add_argument("-y", "--csvinput", help="input csv file", default="Use pythonlist input", type=str)
    parser.add_argument("-n", "--name", help="Name of the input simulation", default="", type=str)
    parser.add_argument("-s", "--csvsinput", help="input csv file in python input", type=str)
    parser.add_argument("-t", "--time", help="How long was the simulation?", default=400, type=float)
    parser.add_argument("-o", "--outputplot", help="Name of output plot", default="RMSIPs.jpg", type=str)
    parser.add_argument("-f", "--figureTitle", help="output figure title", default="RMSIP plot", type=str)
    parser.add_argument("--addPrint", help="Add printing of the process", action="store_true", default=False)

    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)


    run(args)

if __name__=="__main__":
    main()