from __future__ import division
import sys
import networkx as nx
import glob
import operator
import math



in_folder = "/Users/jessicaxinhu/Documents/noa_bc_pilot/diffCoExpress_analysis/diffCoEx/with_vest/diffmodules/diffNetw_inweb"
out_folder = "/Users/jessicaxinhu/Documents/noa_bc_pilot/diffCoExpress_analysis/diffCoEx/with_vest/diffmodules/diffNetw_inweb/network_properties/pos_neg_edges"


G = glob.glob(in_folder + "/*_inweb.tsv")
cutoff_degree = 5
highDegree_nodes = set()


for f in G:
	netw = open(f)
	g = nx.Graph()
	for line in netw:
		tmp = line.strip().split("\t")
		node1 = tmp[0]
		node2 = tmp[1]
		
		g.add_edge(node1, node2, diff=tmp[3])
	
	# print module names as check
	mod_name = f.split("/")[-1].split("_")[0]
	print mod_name
	



        # degree
	#fh1 = open((out_folder + "/" + mod_name + '_degree.txt'),'w')
	deg = nx.degree(g)
	sorted_deg = sorted(deg.items(), key=operator.itemgetter(1), reverse=True)

	#for item in sorted_deg:
	#        fh1.write("%s\t%s\n" % (item[0], item[1]))
	#        if item[1] >= cutoff_degree:
	#                highDegree_nodes.add(item[0])


	# rewiring score
	node_info = {}
	node_info2 = {}
	for n in g.nodes():
		n_edges = g.edges(n)	# all edges for a spec node
		count_pos_edges = 0
		count_neg_edges = 0
		for e in n_edges:
			if g.edge[e[0]][e[1]]['diff'] == '+-':
				count_neg_edges += 1
			elif g.edge[e[0]][e[1]]['diff'] == '-+':
				count_pos_edges += 1

		rewir_score = math.log((count_pos_edges + 1)/(count_neg_edges + 1))	# natural log
		node_info[n] = (deg[n], rewir_score)
		

		node_info2[n] = (count_pos_edges, count_neg_edges) 


	# print pos and neg edges
	#fh1 = open((out_folder + "/" + mod_name + '_NumberPosNegEdges.txt'),'w')
	#for key in node_info2.keys():
	#	fh1.write("%s\t%s\t%s\n" % (key, str(node_info2[key][0]), str(node_info2[key][1])))



	# sort and print node info (degree and rewiring ratio)
	sorted_info = sorted(node_info.items(), key=operator.itemgetter(1), reverse=True)
	fh1 = open((out_folder + "/" + mod_name + '_nodeInfo.txt'),'w')
	for item in sorted_info:
		fh1.write("%s\t%s\t%s\n" % (item[0], str(item[1][0]), str(item[1][1])))



	





#fh2 = open((out_folder + '/high_degree_nodes.txt'),'w')
#for n in highDegree_nodes:
#	fh2.write("%s\n" % (n))
	






