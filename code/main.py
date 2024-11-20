"""
CSE 6140 Fall 2024
TSP Project: Felicity Nielson, Isabel Berry, Austin Chemelli, Prasad Shetye

Solves the Traveling Salesman Problem using three different algorithm options: Exact (Brute-Force), Approximate (Deterministic 2-Approximation Algorithm), or Local Search (a stochastic algorithm)
"""
import math
import time
import argparse
import pandas as pd
import random
import numpy as np
import itertools

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
        file_name = file_name + f'_{seed}'
    file_name = file_name + '.sol'

    with open(file_name, 'w') as output:
        output.write(f'{quality}\n')
        for num, city in enumerate(tour_ordered_list):
            if num ==0:
                output.write(f'{city}')
            else:
                output.write(f', {city}')


def brute_force(cutoff, city_coords):
    """
    Solve the TSP using brute force with a time cutoff.

    Parameters
    ----------
    cutoff: int
        Time cutoff in seconds.
    city_coords: list
        List of city coordinates.

    Returns
    -------
    tuple
        (best_distance, best_tour) where best_distance is the shortest distance found,
        and best_tour is the corresponding tour as an ordered list of city indices.
    """
    base = time.time()
    ret = [math.inf, None]

    def helper(current, remaining):
        if time.time() < cutoff + base:
            if not remaining:
                distance = 0
                for i in range(len(current)):
                    curr1, curr2 = current[i], current[(i+1) % len(current)]
                    distance += math.sqrt((float(curr1[1]) - float(curr2[1])) **2 + (float(curr1[2]) - float(curr2[2]))**2)
                if distance < ret[0]:
                    ret[0] = distance
                    ret[1] = [city[0] for city in current]
                return
            else:
                for ind, city in enumerate(remaining):
                    current.append(city)
                    helper(current, remaining[:ind] + remaining[ind+1:])
                    current.pop()
        else:
            return

    helper([], city_coords)

    return tuple(ret)

def approximate_mst(cutoff, city_coords):
    return 

def local_search(cutoff, city_coords, seed):
    """
    Solves the Traveling Salesman Problem for a given list of city coordinates through a simulated annealing algorithm 
    
    Parameters
    ----------
    cutoff: int
        cutoff time for algorithm to exit after
    city_coords: list 
        list of lists of city coordinates and city numbers

    Returns
    -------
    quality: float
        quality of the solution that was found
    tour_ordered_list: list
        list of the order to visit the nodes for the found solution to TSP
    """

    kT = 500000
    coolingFraction = 0.98 
    steps_to_lower = 10000
    iters = 0
    boltzmann = np.inf
    old_S = list(range(len(city_coords)))
    distances = create_distance_matrix(city_coords)
    best_S = old_S
    best_cost = cost(distances, best_S)

    start = time.time()
    while(float((time.time()-start)) < cutoff):
        index1, index2 = random.sample(range(len(old_S)), 2)
        new_S = old_S.copy()
        new_S[index1] = old_S[index2]
        new_S[index2] = old_S[index1]
        cost_new_S = cost(distances, new_S)
        cost_old_S = cost(distances, old_S)

        if cost_new_S <= cost_old_S:
            old_S = new_S
        else:
            deltaE = cost_new_S - cost_old_S
            #generate random number
            random_float = random.random()
            boltzmann = np.exp(-deltaE/kT)
            if random_float < np.exp(-deltaE/(kT)):
                old_S = new_S
        if cost_old_S < best_cost:
            best_cost = cost_old_S
            best_S = old_S
        if iters%steps_to_lower == 1:
            kT = kT * coolingFraction
        if boltzmann < 0.0001:
            kT = 500000
        iters+=1
    return best_cost, best_S 

def calculate_distance(city1, city2):
    """
    Calculates distance between two cities
    
    Parameters
    ----------
    city1: list
        list of city 1 with first element as city number, second element as X coord, and third element as Y coord
    city2: list
        list of city 2 with first element as city number, second element as X coord, and third element as Y coord

    Returns
    -------
    distance: float
        distance between the two provided cities
    """
    x1, y1 = float(city1[1]), float(city1[2])  # Extracting X, Y coordinates of node1
    x2, y2 = float(city2[1]), float(city2[2])  # Extracting X, Y coordinates of node2
    return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)

def create_distance_matrix(city_coords):
    """
    creates a distance matrix containing the distances between all citiies provided in the city_coords list

    Parameters
    ----------
    city_coords: list
        list of all cities where each city contains a list of the index, x coord, and y coord

    Returns
    -------
    distance_matrix: list of lists
        distance matrix between all cities
    """
    num_cities = len(city_coords)
    distance_matrix = [[0] * num_cities for city in range(num_cities)]  # Initialize a num_nodes x num_nodes matrix

    for i in range(num_cities):
        for j in range(i, num_cities):  # Only compute upper triangle (distance matrix is symmetric)
            dist = calculate_distance(city_coords[i], city_coords[j])
            distance_matrix[i][j] = dist
            distance_matrix[j][i] = dist  # Since distance is symmetric (dist(i,j) == dist(j,i))

    return distance_matrix

def cost(distance_matrix, tour_ordered_list):
    """
    Computes the cost of a possible tour route between cities

    Parameters
    ----------
    tour_ordered_list: list
        a list specifying the order to visit the cities in on the tour
    distance_matrix: list of lists
        distance matrix between all cities

    Returns
    -------
    total_distance: float
        cost of the tour route
    """
    total_distance = 0
    num_cities = len(tour_ordered_list)

    for i in range(num_cities):
        # Get the current node and the next node in the tour
        current_city_index = tour_ordered_list[i] - 1  # Nodes are 1-indexed, so adjust to 0-indexed
        next_city_index = tour_ordered_list[(i + 1) % num_cities] - 1  # Use modulo to loop back to the start

        # Calculate distance between current node and next node
        total_distance += distance_matrix[current_city_index][next_city_index] 

    return total_distance


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
