import sys
import math
import numpy as np
import igraph
import networkx as nx
import networkx.algorithms.link_analysis as la


from helper import Helper

class Quark(Helper):

    def __init__(self, ref, mutation_rate = 10**-6):
        self.expected_mutations = math.ceil(len(ref) * mutation_rate * 2)
        self.simulations = 10000000
        self.p_restart = 0.9

        self.freq = {}
        self.total = 0
    

    def distance_matrix(self, hap_list):

        self.samples = len(hap_list)
        self.personalize = {}
        
        x = []
        for i in range(self.samples):
            y = []
            for j in range(self.samples):
                y.append(0)
            x.append(y)

        self.matrix = np.array(x)
        for n1 in range(self.samples):
            k1, v1 = hap_list[n1]
            for n2 in range(n1):
                k2, v2 = hap_list[n2]

                if n1 != n2:
                    diff = len(k1) + len(k2)
                    for var in k1:
                        if var in k2:
                            diff -= 2

                    self.matrix[n1][n2] = diff
                    self.matrix[n2][n1] = diff

            self.freq[n1] = v1
            self.total += v1

        for n1 in range(self.samples):
            self.personalize[n1] = float(self.freq[n1]) / float(self.total)
    

    def graph_it(self):
        G = nx.Graph()
        G.add_nodes_from(range(len(self.matrix)))

        edges = []
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[i][:i])):
                if i != j:
                    if self.matrix[i][j] <= int(self.expected_mutations):
                        edges.append((i, j))
        G.add_edges_from(edges)
        self.graph = G

    
    def graph_it_igraph(self):
        G = igraph.Graph()
        G.add_vertices(len(self.matrix))

        edges = []
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[i][:i])):
                if i != j:
                    if self.matrix[i][j] <= int(self.expected_mutations):
                        edges.append((i, j))
        G.add_edges(edges)
        self.graph = G
    
    
    def remove_singletons(self):
        pass


    def rank_it(self, out):

        f_out = open(out, "w")

        rank = la.pagerank(self.graph, personalization=self.personalize,
                           max_iter=self.simulations, alpha=self.p_restart)

        size_of_rank = len(rank)

        f_out.write("\t".join(["id", "pagerank", "odds", "freq-abs", "freq-rel"]) + "\n")

        for k, v in sorted(rank.items(), key=lambda q: q[1],
                           reverse=True):
            f_out.write("%s:\t%s\t%s\t%s\t%s\n" % (k, v,
                                                   float(v) * float(size_of_rank),
                                                   self.freq[k],
                                                   self.freq[k]/float(self.total)))

        f_out.close()

    
    def rank_it_igraph(self, out):

        f_out = open(out, "w")
        
        rank = self.graph.pagerank(niter=self.simulations, damping=self.p_restart, directed=False)

        size_of_rank = len(rank)

        f_out.write("\t".join(["id", "pagerank", "odds", "freq-abs", "freq-rel"]) + "\n")
        
        for k, v in sorted(zip(range(self.graph.vcount()), rank), key=lambda q: q[1],
                           reverse=True):
            f_out.write("%s:\t%s\t%s\t%s\t%s\n" % (k, v,
                                                   float(v) * float(size_of_rank),
                                                   self.freq[k],
                                                   self.freq[k]/float(self.total)))

        f_out.close()

    
    def predict_it(self, odds):
        rank = self.graph.pagerank(niter=self.simulations, damping=self.p_restart)

        sum_of_rank = len(rank)

        prediction = []
        new_total = 0

        for k, v in sorted(rank.items(), key=lambda q: q[1], reverse = True):
            if (float(v) * float(sum_of_rank)) > odds:
                prediction.append([k, self.freq[k],
                                   self.freq[k]/float(self.total),
                                   v])

                new_total += self.freq[k]

        # normalize
        for each in prediction:
            each[2] = float(each[1])/float(new_total)

        return prediction
