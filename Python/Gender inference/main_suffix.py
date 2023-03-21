import os
import numpy

from helpers import *
from model_suffix import Model

#-------------------------------------------------------------------------

# path to the current directory
base = os.path.dirname(os.path.realpath(__file__)) + os.sep

# path to the 'name_gender' table
data_path = base + "name_gender.csv"

X = None		# for names
Y = None		# for the corresponding labels

# fraction of the dataset that will be used for training
#
train_fraction = 0.9

"""
experiments showed that using only about 10% of the provided dataset for training
allows to get a smaller model which is comparable to models which use more training data
(in terms of accuracy)
"""

#-------------------------------------------------------------------------

""" data acquisition and pre-processing """

with open(data_path, 'r', encoding = "utf-8") as f:

	lines = f.readlines()

	# convert lines to an array of pairs [name, gender]
	data = list(map(lambda s: s.strip(' \r\n').split(','), lines[1:]))

	# creating training and testing datasets
	#
	X, Y = dataset_split(*unzip(data), train_fraction)

#-------------------------------------------------------------------------

"""  model fitting """

model = Model()

# hyperparameters
#
#    obtained experimentally via hyperparameter tuning
#    for a range of different `train_fraction` values
#
max_suffix_length = 6
threshold = 0.1

# hyperparameter tuning---------------------------------------------
#
if False:

	max_suffix_lengths = range(1, 7)
	thresholds = [0.01, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9]

	accuracies = numpy.zeros((len(max_suffix_lengths), len(thresholds)))
	model_sizes = numpy.zeros((len(max_suffix_lengths), len(thresholds)), dtype = int)

	for idx_0, max_suffix_length in enumerate(max_suffix_lengths):
		for idx_1, threshold in enumerate(thresholds):
		
			model.fit(X['train'], Y['train'], max_suffix_length, threshold)
			
			y_predicted = numpy.array(model.predict(X['test']))

			q = ~(y_predicted == Y['test'])
			error_rate = 100*q.astype(int).sum()/len(q)

			accuracies[idx_0, idx_1] = numpy.round(100 - error_rate, 1)
			model_sizes[idx_0, idx_1] = len(model.suffixes)

	# locate all maxima
	idxs = numpy.argwhere(accuracies == accuracies.max())

	# maxima scores
	# (the goal is to minimize the `suffix_length` and maximize the `threshold`)
	#
	scores = (idxs[:, 0] - idxs[:, 1])

	# locate all minima
	score_idxs = numpy.argwhere(scores == scores.min()).transpose()[0]

	# selecting a pair of indices which minimizes both the score and the model size
	#
	q = idxs[score_idxs].transpose()
	min_idx = numpy.argmin(model_sizes[q[0], q[1]]) 
	p = idxs[score_idxs][min_idx]

	# adjusting hyperparameters
	#
	max_suffix_length = max_suffix_lengths[p[0]]
	threshold = thresholds[p[1]]

	print()
	print(accuracies)

	print()
	print('train_fraction:', train_fraction)
	print()
	print('Optimal model values:')
	print('     max_suffix_length:',  max_suffix_length)
	print('             threshold:',  threshold)
	print('              accuracy: ', accuracies[p[0], p[1]], '%', sep='')
	print('            model size:',  model_sizes[p[0], p[1]])

#-------------------------------------------------------------------

# actual fitting
#
model.fit(X['train'], Y['train'], max_suffix_length, threshold)

# the trained model can be saved and loaded:
#
#model.save(base + 'model.dict')
#model.load(base + 'model.dict')

#-------------------------------------------------------------------------

""" model testing """

y_predicted = numpy.array(model.predict(X['test']))

q = ~(y_predicted == Y['test'])
error_rate = 100*q.astype(int).sum()/len(q)

# percentage of names which were failed to be classified
#
unclassified_rate = 100*(y_predicted == None).astype(int).sum()/len(q)

print()
print('Overall average precision: {0:.2f}%'.format(100 - error_rate))
print('Average unclassified rate: {0:.2f}%'.format(unclassified_rate))
print()
print('Model size: ', len(model.suffixes))

#-------------------------------------------------------------------------

""" [optional] model visualization """

# prints all the suffixes with a corresponding label and score
#
if False:

	print()

	for suffix in model.suffixes:

		# only display short suffixes
		#
		if len(suffix) < 4:
			print(suffix, ': ', model.suffixes[suffix], sep='')
