import matplotlib
matplotlib.use('Agg')
import numpy as np
import sys
import os
# import geode as distance_module
import cartesianDist as distance_module
from jaccard_plot import *
from itertools import product
# from shapely.geometry import LineString, Point
import matplotlib.pyplot as plt
import train_geo
import seaborn

# Try to remove this line eventually!!!
np.seterr(all='ignore')

# Method Definition
# =================

def get_migration_matrix(best_distM):
	migration_matrix = np.zeros((NB_DEMES, NB_DEMES))
	migration_matrix[best_distM <= MAXIMUM_DISTANCE_MIGRATION] = MIGRATION_PARAMETER
	return migration_matrix

def v_color():
    return '#ffcc33'


def plot_point(ax, p):
    x, y = p
    ax.plot(x, y, 'o', color='#999999', zorder=1)


def plot_line(ax, p1, p2):
    x, y = [p1[0], p2[0]], [p1[1], p2[1]]
    ax.plot(x, y, color=v_color(), alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)


def draw_grid(ax, pos, migration_matrix):
	points = [ plot_point(ax, p)  for p in pos ]
	for i in range(len(migration_matrix)):
		for j in range(i, len(migration_matrix)):
			if migration_matrix[i,j] > 0 :
				plot_line(ax, grid[i], grid[j])


def arg_sampling_number_per_demes():
	output_list = [TOTAL_SAMPLE/NB_DEMES] * NB_DEMES
	return '-I %d ' % len(output_list) +  ' '.join(map(str,output_list))


def arg_migration_matrix(M):
	output_list = []
	for i in range(len(M)):
		for j in range(len(M)):
			output_list.append(str(M[i,j]) if i != j else 'x' )

	return '-ma ' + ' '.join(map(str,output_list))


def arg_create():
	args = []
	N_0 = 10**(5)
	NB_BASE_PAIRS = 2e6
	args.append('macs %d %d' % (TOTAL_SAMPLE, NB_BASE_PAIRS))  # first args
	args.append('-t %s' % ((4*N_0)/100 / (1.4 * (10)**8)))  # mutation rate theta = 4N per site per generation
	# args.append('-t %f' % THETA)  # mutation rate theta
	# args.append('-p 9')  # The number of decimals after the dot for precision
	# args.append('-L') # Should give some sort of time.
	args.append('-T')
	# args.append('-r %f %d' % (RECOMBINATION_PARAMETER, NUMBER_OF_SITES ) )  # recombination
	args.append('-r %s' % ((4*N_0)/100 / float(1*(10)**7)) )  # recombination
	args.append('-h 1e2')
	args.append(arg_sampling_number_per_demes())
	args.append(arg_migration_matrix(migration_matrix))
	N_0_per_demes = 1./NB_DEMES
	for i in range(NB_DEMES):
		args.append('-n %d %f' % (i + 1, N_0_per_demes))
	# args.append('-R example_input/hotspot.txt')
	# args.append('-F example_input/ascertainment.txt 0')
	args.append('2>trees.txt | msformatter > haplotypes.txt')

	f = open('c.macs', 'w')
	for l in args:
		f.write('%s ' % l)
	f.close()
	f = open('markers.txt', 'w')
	markers = reduce(lambda x, y: x+y, [[i] * (TOTAL_SAMPLE/NB_DEMES) for i in range(NB_DEMES)] )
	for m in markers:
		f.write('%d\n' % m )
	f.close()

# Program setup
# =============
searchsphere_no = 17
MAXIMUM_DISTANCE_MIGRATION = 1.1 # in km
# searchspace_euclidean, searchspace_radians, searchspace_degree = train_geo.load_and_filter_searchspace(BOUNDS["SCCS"], searchsphere_no, prefix_path='../myscripts/gen_searchspheres/searchspheres/')

grid = np.array([(i, j) for i in range(-7, 7) for j in range(-5, 5) ])

# NB_DEMES = len(searchspace_radians)
NB_DEMES = len(grid)
SAMPLE_PER_DEME = 1
TOTAL_SAMPLE = SAMPLE_PER_DEME * NB_DEMES
DIPLOID_POPULATION_SIZE = 10000  # AKA N_0
NUMBER_OF_SITES = 4000  # The number of places recombination can happen
# NEUTRAL_MUTATION_RATE = 3 *  (1. / NUMBER_OF_SITES)  # Per site
BETWEEN_SITES_RECOMBINATION_PROBABILITY = 10**(-6)  # 1. / NUMBER_OF_SITES
RECOMBINATION_PARAMETER = 4 * DIPLOID_POPULATION_SIZE * (BETWEEN_SITES_RECOMBINATION_PROBABILITY * NUMBER_OF_SITES)  # 4N_0r, where r is the probability of cross-over per generation between the ends of the locus being simulated.
# THETA = 4 * DIPLOID_POPULATION_SIZE * NEUTRAL_MUTATION_RATE  # 4*N_0*mu


# reorder grid according to longitudes
# sorter = np.argsort(searchspace_radians[:,0])
# searchspace_euclidean, searchspace_radians, searchspace_degree = searchspace_euclidean[sorter], searchspace_radians[sorter], searchspace_degree[sorter]

indiv_position = []
# for p in searchspace_radians:
# 	for i in range(TOTAL_SAMPLE / len(searchspace_radians)):
# 		indiv_position.append(p)

for p in grid:
	for i in range(TOTAL_SAMPLE / NB_DEMES ):
		indiv_position.append(p)

np.savetxt('macs_sim/indiv_position.txt', indiv_position, delimiter='\t')


# Train3
# TRAIN_PATH = '/Volumes/gravel/barakatt_projects/dist/src/simulations/trains/train0.txt'
# TRAIN_PATH = '../dist/src/simulations/trains/train1.txt'
# mytrain = geode.line(*np.loadtxt(TRAIN_PATH), inputformat='degrees')
mytrain = distance_module.line(0, 0, 10, 0)


m_param = 0.000025  #0.05  # m is the fraction of each subpopulation made up of new migrants each generation.
mu_param = 1. / NUMBER_OF_SITES  # neutral mutation rate for the entire genome. I think 1 would mean 1 mutation per generation on average
THETA = 4 * DIPLOID_POPULATION_SIZE * mu_param

# distM = distance_module.get_dist(*searchspace_radians.T)
distM = distance_module.get_dist(*zip(*grid))

best_distM = distance_module.computeDistances(mytrain, distM, *np.array(grid).T)
MIGRATION_PARAMETER = 4 * DIPLOID_POPULATION_SIZE * m_param

# Some book keeping things

# ref_matrix = np.arange(NB_DEMES).reshape((5,5))
# d = {}
# for i,p in enumerate(product(range(np.array(5)), range(np.array(5)))):
# 	d[i] = p

# for i,j in product(range(np.array(NB_DEMES)), range(np.array(NB_DEMES))):
# 	coor1 = d[i]
# 	coor2 = d[j]
# 	if np.absolute(coor1[0] - coor2[0]) + np.absolute(coor1[1] - coor2[1]) == 1:
# 		migration_matrix[j,i] = 1.0
# 		print j, i, coor1, coor2

migration_matrix = get_migration_matrix(best_distM)
fig = plt.figure(1)
ax = fig.add_subplot(111)
fit = lambda x, y : (x - 0.05 * (y - x), y + 0.05 * (y - x))
ax.set_xlim(fit(np.min(np.array(grid)[:,0]), np.max(np.array(grid)[:,0])))
ax.set_ylim(fit(np.min(np.array(grid)[:,1]), np.max(np.array(grid)[:,1])))
draw_grid(ax, grid, migration_matrix)
plt.savefig('grid_look.png')
arg_create()
# fcommand = open('run_c.sh','w')
# fcommand.write('ms %s 1 -f c.ms' % TOTAL_SAMPLE)
# fcommand.close()
