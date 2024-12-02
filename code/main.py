"""
CSE 6140 Fall 2024
TSP Project: Felicity Nielson, Isabel Berry, Austin Chemelli, Prasad Shetye

Solves the Traveling Salesman Problem using three different algorithm options:
Exact (Brute-Force), 
Approximate (Deterministic 2-Approximation Algorithm), 
or Local Search (a stochastic algorithm)
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
    Parses command line arguments 
    Instance, algorithm, and time are required, but seed is optional and defaulted to None since it is only used in LS
    Options for algorithm are brute force, approximate, or local search

    example command:
        python main.py -inst Cincinnati -alg LS -time 20

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
    parser.add_argument("-alg", type = str, choices = ['BF', 'LS', 'Approx'], required=True, help='search method (BF, Approx, LS)')
    parser.add_argument("-time", type = int, required = True, help = 'cutoff time in seconds')
    parser.add_argument("-seed", default = None, type = int, help = 'random seed for LS alg')
    return parser.parse_args()

def read_inputfile(filename):
    """
    Reads the input file and creates a list of the nodes and coordinates in the input file

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
    Writes the output file for a given set of parameters with a tour and cost already computed

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
    #Base time for cutoff condition
    base = time.time()
    #Stores the return values [shortest distance found thus far, the order of the shortest tour thus far]
    ret = [math.inf, None]
    #Helper to explore all permutations
    def helper(current, remaining):
        if time.time() < cutoff + base:
            #If all cities explored calculate tour distance
            if not remaining:
                distance = 0
                for i in range(len(current)):
                    curr1, curr2 = current[i], current[(i+1) % len(current)]
                    distance += math.sqrt((float(curr1[1]) - float(curr2[1])) **2 + (float(curr1[2]) - float(curr2[2]))**2)
                # If the current tour is less than the shortest tour we are storing, update the shortest tour we store
                if distance < ret[0]:
                    ret[0] = distance
                    ret[1] = [city[0] for city in current]
                return
            else:
                #Explore the unexplored cities
                for ind, city in enumerate(remaining):
                    current.append(city)
                    # Recursively call helper with the updated path and remaining cities
                    helper(current, remaining[:ind] + remaining[ind+1:])
                    #Backtrack be removing the last city we added 2 lines above
                    current.pop()
        else:
            #Cutoff time exceeded
            return
    helper([], city_coords)

    return tuple(ret)

def approximate_mst(cutoff, city_coords):
    ########################### Define the helper methods ##################################
    # MST-Prim
    def mst_prim(distance_matrix, num_vertices, r):
        # initialize data structures
        visited = [False] * num_vertices
        pq = [i for i in range(num_vertices)]
        key = [math.inf] * num_vertices
        key[r] = 0
        parent = [-1] * num_vertices

        # pop one vertex off at a time from pq until all vertices are considered
        while pq:
            # extract the vertex u with the minimum key value and remove from pq
            u = min(pq, key=lambda x: key[x])
            pq.remove(u)
            visited[u] = True

            # consider each vertex v adjacent to u
            for v in range(num_vertices):
                if not visited[v]:
                    # if this edge is shorter than the current shortest edge to v
                    if distance_matrix[u][v] < key[v]:
                        # attach to the MST via parent
                        parent[v] = u
                        # update the key value of v
                        key[v] = distance_matrix[u][v]

        # create an adjecency list representation for the MST
        mst = [[] for _ in range(num_vertices)]
        for v in range(1, num_vertices):
            u = parent[v]
            mst[v].append(u)
            mst[u].append(v)

        return mst
    
    # preorder walk
    def preorder(curr, visited, walk):
        walk.append(curr + 1)
        visited[curr] = True
        for adj in mst[curr]:
            if not visited[adj]:
                preorder(adj, visited, walk)

    ########################### Generate the tour ##################################

    start_time = time.time()

    # get the distance matrix using the city_coords
    distance_matrix = create_distance_matrix(city_coords)
    num_vertices = len(city_coords)

    # choose a vertex r as the root of the MST
    r = 0

    # get the MST using Prim's algorithm
    mst = mst_prim(distance_matrix, num_vertices, r)

    # preorder walk of MST and make a list h of vertices according to first visited
    visited = [False] * num_vertices
    h = []
    preorder(0, visited, h)

    end_time = time.time()

    # get the quality and runtime of the tour
    quality = cost(distance_matrix, h)
    runtime = end_time - start_time

    return quality, runtime, h

def local_search(cutoff, city_coords, seed, kT=500000, coolingFraction=0.98, steps_to_lower=10000):
    """
    Solves the Traveling Salesman Problem for a given list of city coordinates through a simulated annealing algorithm 
    
    Parameters
    ----------
    cutoff: int
        cutoff time for algorithm to exit after
    city_coords: list 
        list of lists of city coordinates and city numbers
    seed: int 
        random seed
    kT: int
        boltzmann's costant * temperature
    coolingFraction: float
        the percentage of kT the system cools at each step. Used to modify kT at each step
    steps_to_lower: int
        number of steps that will be taken to lower kT to 0 

    Returns
    -------
    quality: float
        quality of the solution that was found
    tour_ordered_list: list
        list of the order to visit the nodes for the found solution to TSP
    """
    start = time.time()
    iters = 0
    boltzmann = np.inf
    old_S = list(range(len(city_coords)))
    distances = create_distance_matrix(city_coords)
    best_S = old_S
    best_cost = cost(distances, best_S)
    random.seed(seed)
    kT_mod = kT

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
            boltzmann = np.exp(-deltaE/kT_mod)
            if random_float < np.exp(-deltaE/(kT_mod)):
                old_S = new_S
        if cost_old_S < best_cost:
            best_cost = cost_old_S
            best_S = old_S
        if iters%steps_to_lower == 1:
            kT_mod = kT_mod * coolingFraction
        if boltzmann < 0.0001: # this number was estimated through trial and error as a good cutoff
            kT_mod = kT
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
    Solves the Traveling Salesman Problem using three different algorithm options: 
    Exact (Brute-Force), 
    Approximate (Deterministic 2-Approximation Algorithm), 
    or Local Search (a stochastic algorithm)
    
    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    args = get_arguments()
    city_coords = read_inputfile(f'../data/{args.inst}.tsp')
    if args.alg == 'BF':
        quality, tour_ordered_list = brute_force(args.time, city_coords)
    elif args.alg == 'LS':
        quality, tour_ordered_list = local_search(args.time, city_coords, args.seed)
    else:
        quality, runtime, tour_ordered_list = approximate_mst(args.time, city_coords)

    write_output(args.inst, args.alg, args.time, quality, tour_ordered_list, args.seed)
    print(args.inst, args.time, quality)

def test_LS():
    """
    simple logic to automate the testing. this is just for me and felicities sake, we will manipulate the df further to make it up to quality of the deliverable, but i figured i would push it so that y'all could do something similar if you want. but do be aware I have not tested this yet
    """
    cities = ['Atlanta', 'Berlin', 'Boston', 'Champaign', 'Cincinnati', 'Denver', 'NYC', 
              'Philadelphia', 'Roanoke', 'SanFransisco', 'Toronto', 'UKansasState', 'UMissouri']
    times = [15,30,60,300]
    num_iter = 10    

    
    df = pd.DataFrame(index=cities, columns=times)

    for city in cities:
        city_coords = read_inputfile(f'../data/{city}.tsp')
        for time in times:
            quality_sum = 0
            for i in range(0,num_iter):
                seed = random.randint(0,1000)
                quality, tour_ordered_list = local_search(time, city_coords, seed)
                quality_sum += float(quality)
                write_output(city, 'LS', time, quality, tour_ordered_list, seed)
            quality_ave = quality_sum/num_iter
            df.loc[city, time] = quality_ave
            df.to_csv('ls_results.csv')
    df['full tour'] = 'yes'
    df.to_csv('ls_results.csv')


def hp_tune_LS():
    """
    Applies grid search to hyper parameters for LS. 
    """
    # cities = ['Atlanta', 'Berlin', 'Boston', 'Champaign', 'Cincinnati', 'Denver', 'NYC', 
    #           'Philadelphia', 'Roanoke', 'SanFransisco', 'Toronto', 'UKansasState', 'UMissouri']
    cities = ['Toronto']
    times = [30]
    num_iter = 1

    grid = [[400000, 500000, 600000],[0.9, 0.95, 0.98],[8000, 10000, 12000]] # columns: kT, coolingFraction, steps_to_lower (starting: 500000, 0.98, 10000)
    
    df = pd.DataFrame(index=cities, columns=times)

    for kT in grid[0]:
        for coolingFraction in grid[1]:
            for steps_to_lower in grid[2]:
                for city in cities:
                    city_coords = read_inputfile(f'{city}.tsp')
                    for time in times:
                        quality_sum = 0
                        for i in range(0,num_iter):
                            #seed = random.randint(0,1000)
                            seed = 0
                            quality, tour_ordered_list = local_search(time, city_coords, seed, kT, coolingFraction, steps_to_lower)
                            quality_sum += float(quality)
                            end_tag = f'{seed}_kT{kT}_CF{coolingFraction}_StL{steps_to_lower}'
                            write_output(city, 'LS', time, quality, tour_ordered_list, end_tag)
                        quality_ave = quality_sum/num_iter
                        df.loc[city, time] = quality_ave
                df['full tour'] = 'yes'
                df.to_csv('ls_results.csv')

def test_Approx():
    cities = ['Atlanta', 'Berlin', 'Boston', 'Champaign', 'Cincinnati', 'Denver', 'NYC', 
              'Philadelphia', 'Roanoke', 'SanFrancisco', 'Toronto', 'UKansasState', 'UMissouri']
    times = [5,15,30,60,120,180,300]
    
    # create a df where rows are the cities and the column are runtime and quality
    df = pd.DataFrame(index=cities, columns=['time', 'quality', 'runtime', 'full tour'])

    for city in cities:
        city_coords = read_inputfile(f'..\\data\\{city}.tsp')
        for time in times:
            quality, runtime, tour_ordered_list = approximate_mst(time, city_coords)
            write_output(city, 'Approx', time, quality, tour_ordered_list)
            df.loc[city, 'time'] = runtime
            df.loc[city, 'quality'] = quality
            df.loc[city, 'runtime'] = runtime
            if runtime < time:
                df.loc[city, 'full tour'] = 'yes'
            else:
                df.loc[city, 'full tour'] = 'no'
            
    df.to_csv('approx_results.csv')

def test_BF():
    cities = ['Atlanta', 'Berlin', 'Boston', 'Champaign', 'Cincinnati', 'Denver', 'NYC', 
              'Philadelphia', 'Roanoke', 'SanFrancisco', 'Toronto', 'UKansasState', 'UMissouri']
    cut_offs = [5, 15, 30, 60, 120, 180, 300]
    
    df = pd.DataFrame(index=cities, columns=['time', 'quality', 'runtime', 'full tour'])

    for city in cities:
        city_coords = read_inputfile(f'../data/{city}.tsp')
        for cut_off in cut_offs:
            start_time = time.time()
            quality, tour_ordered_list = brute_force(cut_off, city_coords)
            end_time = time.time()
            runtime = end_time - start_time
            write_output(city, 'BF', cut_off, quality, tour_ordered_list)

            df.loc[city, 'time'] = cut_off
            df.loc[city, 'quality'] = quality
            df.loc[city, 'runtime'] = runtime
            df.loc[city, 'full tour'] = 'yes' if runtime < cut_off else 'no'    
    df.to_csv('bf_results.csv')


if __name__ == "__main__":
    main()
