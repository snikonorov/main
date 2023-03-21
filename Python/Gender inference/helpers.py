""" helper functions """

import numpy
from numpy.random import permutation

# splits arrays `x` and `y` into a pair of (training, testing) arrays each;
# returns a pair of dictionaries of the form {'train': [...], 'test': [...]}
#
def dataset_split(x, y, train_fraction = 0.9, shuffle = True):

	X = {'train': None, 'test': None}
	Y = {'train': None, 'test': None}

	idxs = numpy.array(range(len(x)))

	if shuffle:
		idxs = permutation(idxs)

	idx = round(len(x)*train_fraction)

	X['train'], X['test'] = x[idxs][:idx], x[idxs][idx:]
	Y['train'], Y['test'] = y[idxs][:idx], y[idxs][idx:]

	return (X, Y)

# 'unzip' an array (convert from array of pairs to a pair of arrays)
#
def unzip(a):

	return tuple \
	(
		map(numpy.array, [q for q in zip(*a)])
	)

# returns the 'inverse' of a dict `d`,
# i.e. a dict which is a mapping from `d.values()` to `d.keys()`
#
# if there are multiple keys corresponding to the same `value`
# the resulting dict will map the `value` to a tuple of those keys
#
def inv_dict(d):

	res = dict()

	for key in d:
	
		value = d[key]

		if value not in res:
			res[value] = key
		else:
			res[value] = (*res[value], key)

	return res
