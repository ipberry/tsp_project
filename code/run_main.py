import os
import subprocess

# Define the datasets and times
instances = [
    "Atlanta", "Berlin", "Boston", "Champaign", "Cincinnati",
    "Denver", "NYC", "Philadelphia", "Roanoke", "SanFrancisco",
    "Toronto", "UKansasState", "UMissouri"
]

times = [15, 30, 60, 300]

# Path to the main.py file
main_py = "main.py"

# Iterate through instances and times
for inst in instances:
    for time in times:
        command = [
            "python", main_py,
            "-inst", inst,
            "-alg", "BF",
            "-time", str(time)
        ]
        print(f"Running: {' '.join(command)}")
        try:
            # Run the command
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error while executing: {' '.join(command)}")
            print(e)
