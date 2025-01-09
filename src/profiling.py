import pickle
import networkx as nx 

with open('ITS_graphs.pkl.gz', 'rb') as file: 
	data = pickle.load(file)


# generate a list of reaction graphs
# find a way to recover the reaction id from the reaction centre 
	# consider adding key to data dict with rc object 
	# pickle for quick 



# retrieve reaction centre 
# in next iteration consider making it independent of reaction centre method.
# i.e. to include WP6 neighbourhood function 

# cluster method should receive the list of dicts and return a copy of this 
# good so that we can still relate the id of the reaction to the reaction centre and the cluster 
# to which it belongs