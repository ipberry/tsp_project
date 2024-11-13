"""
CSE 6140 Fall 2024
TSP Project: Felicity Nielson, Isabel Berry, Austin Chemelli, Prasad Shetye

Solves the Traveling Salesman Problem using three different algorithm options: Exact (Brute-Force), Approximate (Deterministic 2-Approximation Algorithm), or Local Search (a stochastic algorithm)
"""

import argparse
import pandas as pd

def get_arguments():
    """
    this function will create and parse the command line arguments 
    Instance, algorithm, and time are required, but seed is optional and defaulted to None since it is only used in LS
    Options for algorithm can be brutre force, local search, or approximate

    Parameters
    ----------
    None

    Returns
    -------
    parse.parse_args
        inputted arguments parsed and accessible by saying args.inst etc
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("-inst", type = str, required=True, help = 'filename of the dataset')
    parser.add_argument("-alg", type = str, choices = ['BF', 'LS', 'Approx'], required =True, help = 'search method (BF, Approx, LS)')
    parser.add_argument("-time", type = int, required = True, help = 'cutoff time in seconds')
    parser.add_argument("-seed", default = None, type = int, help = 'random seed for LS alg')
    return parser.parse_args()

def read_inputfile(filename):
    """
    this function will read the input file and create a list of the nodes and coordinates in the input file

    Parameters
    ----------
    filename: str
        filename of the input instance (without the extension)

    Returns
    -------
    cities: list
        Dimension num_nodes x num_coords+1, a list of a list for each node where the first entry is the node number and the rest of the entrys are the coordinates
    """
    
    with open(filename) as inp_file:
        inp_lines = inp_file.readlines()
    start = False
    cities = []
    for num, line in enumerate(inp_lines):
        if 'EOF' in line:
            break
        if start:
            cities.append(line.split())
        if 'NODE_COORD' in line:
            start = True
    return cities

def write_output(instance, method, cutoff, quality, tour_ordered_list, seed = None):
    """
    writes the output file for a given set of parameters with a tour and cost already computed

    Parameters
    ----------
    instance: str
        filename of the input instance (without the extension)
    method: str
        BF, LS, or Approx for the type of algorithm used to find the tour
    cutoff: int
        cutoff time for the algorithm
    quality: float
        quality of the tour found
    tour_ordered_list: list
        ordered list of the tour found

    Returns
    -------
    None
    """
    file_name = f'{instance}_{method}_{cutoff}'
    if seed is not None:
        file_name = file_name + '_seed'
    file_name = file_name + '.sol'

    with open(file_name, 'w') as output:
        output.write(f'{quality}\n')
        for num, city in enumerate(tour_ordered_list):
            if num ==0:
                output.write(f'{city}')
            else:
                output.write(f', {city}')

def brute_force(cutoff, city_coords):
    return quality, tour_ordered_list 

def approximate_mst(cutoff, city_coords):
    return quality, tour_ordered_list 

def local_search(cutoff, city_coords, seed):
    return quality, tour_ordered_list 

def main():
    """
    Solves the Traveling Salesman Problem using three different algorithm options: Exact (Brute-Force), Approximate (Deterministic 2-Approximation Algorithm), or Local Search (a stochastic algorithm)
    
    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    args = get_arguments()
    city_coords = read_inputfile(f'{args.inst}.tsp')
    if args.alg == 'BF':
        quality, tour_ordered_list = brute_force(args.time, city_coords)
    elif args.alg == 'LS':
        quality, tour_ordered_list = local_search(args.time, city_coords, args.seed)
    else:
        quality, tour_ordered_list = approximate_mst(args.time, city_coords)
    write_output(args.inst, args.alg, args.time, quality, tour_ordered_list, args.seed)

if __name__ == "__main__":
    main()
