import dendropy as dy
import numpy as np

fin = open('haplotypes.txt', 'r')


while fin.readline() != "//\n":
	continue



t = fin.readline().split(']')[1][0:-2]
tree1 = dy.Tree.get_from_string(t, schema="newick")

print tree1
