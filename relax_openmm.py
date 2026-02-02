from openmm.app import *
from openmm import *
from openmm.unit import *
import numpy as np
import sys

gro_file = "last_frame_pull1.gro"  # Replace with your GRO file
top_file = "topol.top"  # Replace with your topology file


# Define DistanceReporter Class
class DistanceReporter:
    def __init__(self, file_name, reportInterval, group1, group2, initial_distance):
        """
        A custom reporter to log the distance and strain between two groups.

        Parameters:
        - file_name (str): The name of the output file.
        - reportInterval (int): Number of steps between reports.
        - group1 (list): Atom indices for group 1.
        - group2 (list): Atom indices for group 2.
        - initial_distance (float): The initial distance between the two groups in nanometers.
        """
        self.file_name = file_name
        self.reportInterval = reportInterval
        self.group1 = group1
        self.group2 = group2
        self.initial_distance = initial_distance

        # Initialize the output file
        with open(self.file_name, "w") as f:
            f.write("Step, Distance (nm), Strain (%)\n")

    def describeNextReport(self, simulation):
        """Describe the next report this reporter will generate."""
        return (True, True, False, False, self.reportInterval)

    def report(self, simulation, state):
        """Log the distance and strain between two groups at the current step."""
        # Get current positions
        positions = state.getPositions(asNumpy=True)

        # Compute centroids
        centroid1 = np.mean([positions[i] for i in self.group1], axis=0)
        centroid2 = np.mean([positions[i] for i in self.group2], axis=0)

        # Compute current distance
        current_distance = np.linalg.norm(centroid1 - centroid2) * nanometer  # Add units

        # Calculate strain as a percentage
        strain_percent = ((current_distance - self.initial_distance) / self.initial_distance) * 100

        # Log the step, distance, and strain
        with open(self.file_name, "a") as f:
            f.write(f"{simulation.currentStep}, {current_distance / nanometer:.4f}, {strain_percent:.2f}\n")


group1 = [15071, 15072, 15073, 15074, 15075, 15076, 15077, 15078, 15079, 15080, 15081, 15082, 15083, 15084, 15085, 15086, 15087, 15088, 15089, 15090]
group2 = [60893, 60894, 60895, 60896, 60897, 60898, 60899, 60900]




# Load the GROMACS files
gro = GromacsGroFile(gro_file)
top = GromacsTopFile(top_file, periodicBoxVectors=gro.getPeriodicBoxVectors(), includeDir='/home/liuchw/Softwares/gromacs-2021.3/share/top')

# Calculate the initial distance between the two groups
positions = gro.positions  # Assuming 'gro' is the GromacsGroFile object

def calculate_centroid(positions, group_indices):
    return np.mean([positions[i] for i in group_indices], axis=0)

# Ensure positions are numeric values
centroid1 = calculate_centroid(positions, group1)
centroid2 = calculate_centroid(positions, group2)

centroid1 = np.array([c.value_in_unit(nanometer) for c in centroid1])
centroid2 = np.array([c.value_in_unit(nanometer) for c in centroid2])

initial_distance = np.linalg.norm(centroid1 - centroid2) * nanometer


# Load input files
system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)

# Create simulation
integrator = LangevinIntegrator(300*kelvin, 1.0/picoseconds, 2.0*femtoseconds)
integrator.setConstraintTolerance(0.00001)
platform = Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'single'}
simulation = Simulation(top.topology, system, integrator, platform, properties)
simulation.context.setPositions(gro.positions)
simulation.minimizeEnergy()

max_steps = 500000  # Example total steps
check_interval = 1000

# Add reporters
simulation.loadCheckpoint('relax1_7.chk')
simulation.reporters.append(DCDReporter('relax1_8.dcd', 1000))
simulation.reporters.append(StateDataReporter('relax1_8.log', 1000, step=True, potentialEnergy=True, temperature=True,  volume=True, progress=True, remainingTime=True, speed=True, totalSteps=max_steps))
distance_file = "distance_relax1_8.log"
distance_reporter = DistanceReporter(distance_file, 1000, group1, group2, initial_distance)
simulation.reporters.append(distance_reporter)
simulation.reporters.append(CheckpointReporter('relax1_8.chk', 1000))

print('Running dynamics')
simulation.step(max_steps)

