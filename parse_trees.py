import dendropy as dy
from time import time
import numpy as np
import matplotlib.pylab as plt

fin = open('only_trees.txt', 'r')
fpositions = open('trees.txt', 'r')

while fpositions.readline().split(':')[0] != "Time for build prior tree":
	continue

recombination_position = []
for l in fpositions:
	if l[0:2] == 'Tr':
		recombination_position.append(float(l.split(',')[1].split(':')[1]))

recombination_position = np.array(recombination_position)

query_list = np.linspace(0,np.max(recombination_position),400)
# print 'query_list:',  query_list

indices_trees_of_interest = np.unique(np.searchsorted(recombination_position, query_list))

# while fin.readline() != "//\n":
# 	continue

class mynode:
	def __init__(self, t):
		d = 0
		triple = []
		open_pos = []
		for i in range(len(t)):
			if t[i] == '(':
				open_pos.append(i)
				d += 1
			elif t[i] == ')':
				triple.append((d,open_pos.pop(),i))
				d -= 1
		print triple



# def tree2dict(tree_as_text):


trees = []
c = 0
trees_of_interest = np.array(fin.readlines())[indices_trees_of_interest]
print len(trees_of_interest)
for l in trees_of_interest:
	# t = fin.readline().split(']')[1][0:-1]
	# t = fin.readline()

	trees.append(dy.Tree.get_from_string(l, schema="newick"))
	# print len(trees[-1].taxon_set.labels())
	c += 1
	if c > 600:
		break

NUMBER_LEAFS = len(trees[0].taxon_set.labels())


# Precomputing time to tip for each node of the tree for each trees.
t0 = time()
for t in trees:
	t.calc_node_ages(check_prec=1e-03)

print time() - t0

# print tree1.taxon_set.labels()
t0 = time()

ys = np.zeros((len(trees), NUMBER_LEAFS, NUMBER_LEAFS))

for t in range(len(trees)):
	for i in range(NUMBER_LEAFS):
		for j in range(i + 1, NUMBER_LEAFS):
			ys[t,i,j] = trees[t].mrca(taxon_labels=[str(j),str(i)]).age

print time() - t0
f = plt.figure()
ax = f.add_subplot(111)
ax.set_yscale('log')

ax.plot(recombination_position[indices_trees_of_interest], ys, '-')
plt.show()
# print [tree1.mrca(taxon_labels=['1',str(i)]).edge_length for i in range(100)]
# print tree1.mrca(taxon_labels=tree1.taxon_set.labels()).edge_length


