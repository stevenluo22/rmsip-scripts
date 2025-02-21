import MDAnalysis as mda
import MDAnalysis.analysis.pca as pca
from MDAnalysis.analysis import align
import numpy as np
import argparse

def run(args):
    # Select backbone atoms
    selection = "name CA"

    n_components = args.components  # Use top 10 eigenvectors
    frame_step = args.interval  # Analyze every 20 frames

    rmsip_values = []

    # Load the trajectory
    native = args.pdb
    movie = args.dcd
    u = mda.Universe(native, movie)  # Use "native.pdb" if needed
    
    atoms = u.select_atoms(selection)
    
    rmsip_data = []
        
    # Split the trajectory into segments
    for ts in range(2, len(u.trajectory), frame_step):
        start = 2
        half = int(min((start + ts + frame_step)/2, len(u.trajectory)))
        end = int(min(ts + frame_step, len(u.trajectory)))
    
        # Collect atomic positions for this segment
        first_half_coordinates = [atoms.positions.copy() for ts in u.trajectory[start:half]]
        second_half_coordinates = [atoms.positions.copy() for ts in u.trajectory[half:end]]
        
        # Convert to NumPy array
        first_half_coordinates = np.array(first_half_coordinates)
        second_half_coordinates = np.array(second_half_coordinates)
    
        # Create a new Universe for this segment
        first_half_u = mda.Merge(atoms)  # Clone atom selection
        first_half_u.load_new(first_half_coordinates, format="memory")  # Load coordinates
        
        # Run PCA on this segment
        first_half_pca_analysis = pca.PCA(first_half_u, select=selection).run()
    
        # Store top 10 eigenvectors (each column is an eigenvector)
        first_half_eigenvectors = first_half_pca_analysis.p_components[:, :n_components]
        
        second_half_u = mda.Merge(atoms)
        second_half_u.load_new(second_half_coordinates, format="memory")

        second_half_pca_analysis = pca.PCA(second_half_u, select=selection).run()
        
        second_half_eigenvectors = second_half_pca_analysis.p_components[:, :n_components]

        rmsip = np.sqrt(np.sum(np.dot(first_half_eigenvectors.T, second_half_eigenvectors) ** 2) / n_components)
        if args.addPrint:
            print(f"frame start analysis: {start}, halfway: {half}, end: {end}, RMSIP: {rmsip}")
        rmsip_data.append(rmsip)
    rmsip_values.append(rmsip_data)

    import matplotlib.pyplot as plt

    time = 0
    times = []
    for i in range(1, int(len(u.trajectory)/frame_step)+1):
        time += args.time/frame_step
        times.append(time)

    plt.figure(figsize=(10, 6))

    # Plot each run
    r = 0
    for rmsip_data in rmsip_values:
        r += 1
        plt.plot(times, rmsip_data, label=f"run {r}", marker='o')

    # Customize the plot
    plt.title(args.figureTitle)
    plt.xlabel('Time')
    plt.ylabel('RMSIP')
    plt.xticks(rotation=90)  # Rotate x-axis labels for better readability
    plt.legend(title='Runs', bbox_to_anchor=(1, 1), loc='upper left')
    plt.tight_layout()

    plt.savefig(args.outputplot, dpi=300, bbox_inches='tight')
    # Show the plot

def main(args=None):
    parser = argparse.ArgumentParser(
        description="Calculating Cross-Q/Mutual-Q of pdb files")
    parser.add_argument("-p", "--pdb", help="Where is reference pdb?", type=str)
    parser.add_argument("-d", "--dcd", help="Where is trajectory dcd?", type=str)
    parser.add_argument("-c", "--components", help="How many top PCA components to use?", default=10, type=int)
    parser.add_argument("-s", "--interval", help="analysis interval in frames", default=20, type=int)
    parser.add_argument("-t", "--time", help="How long was the simulation?", default=400, type=float)
    parser.add_argument("-v", "--outputCSV", help="Name of csv output data file", default="RMSIP.csv", type=str)
    parser.add_argument("-o", "--outputplot", help="Name of output plot", default="RMSIP.jpg", type=str)
    parser.add_argument("-f", "--figureTitle", help="output figure title", default="RMSIP plot", type=str)
    parser.add_argument("--addPrint", help="Add printing of the process", action="store_true", default=False)

    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)


    run(args)

if __name__=="__main__":
    main()