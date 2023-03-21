import os
import numpy

from helpers import *
from model_cnn import Model

#-------------------------------------------------------------------------

# suppress tensorflow warnings
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

# path to the current directory
base = os.path.dirname(os.path.realpath(__file__)) + os.sep

""" path to the 'name_gender' table """
data_path = base + "name_gender.csv"

#-------------------------------------------------------------------------

alphabet_size = 26			# assuming English alphabet by default
							# (may be adjusted during pre-processing)

max_word_length = 18		# (may be adjusted during data pre-processing)

train_fraction = 0.9

# data acquisition and pre-processing-------------------------------------

X = None		# for names
Y = None		# for the corresponding labels

with open(data_path, 'r', encoding = "utf-8") as f:

	lines = f.readlines()

	#----------------------------------------------

	def split_and_normalize(s):

		p = s.strip(' \r\n').split(',')
		p[0] = Model.normalize(p[0])

		return p

	#----------------------------------------------

	# convert lines to an array of pairs [name, gender] and normalize the names
	#
	data = list(map(split_and_normalize, lines[1:]))

	# [optional] adjusting the alphabet size
	#
	#     (should be used if non-ASCII characters may be present in data)
	#
	if True:

		# convert all names to bytearrays and shift all the codes by 'A'
		# (so that [A..Z] is converted to [0..25])
		#
		A_code = bytearray('A', "utf-8")[0]
		q = list(map(lambda p: numpy.array(bytearray(p[0], "utf-8")) - A_code, data))

		alphabet = set(numpy.hstack(q))
		alphabet_size = max(alphabet_size, max(alphabet) + 1)

	# adjusting `max_word_length`

	data_max_word_length = max(map(lambda p: len(p[0]), data))
	max_word_length = max(max_word_length, data_max_word_length)

	# [optional] printing the dataset stats
	#
	if True:
	
		print()
		print('Dataset stats:')
		print('   Name count      =', len(data))
		print('   Max word length =', data_max_word_length)
		print('   Alphabet size   =', alphabet_size)
		print()

	# creating training and testing datasets
	#
	X, Y = dataset_split(*unzip(data), train_fraction)

# model fitting-----------------------------------------------------------

# hyperparameters
#
params = \
{
	'n_filters': 128,
	'kernel_sizes': [(17, 9), (7, 7)],		# these values were obtained via hyperparameter tuning
	'dense_dim': 32
}

model = Model(alphabet_size, max_word_length, **params)

# hyperparameter tuning-----------------------------------------------
#
if False:

	# tuning `kernel_sizes`
	#
	p0 = [(11, 3), (13, 5), (15, 7), (17, 9), (19, 11)]
	p1 = [(3, 3), (5, 5), (7, 7)]

	accuracies = numpy.zeros((len(p0), len(p1)))

	for idx0 in range(len(p0)):
		for idx1 in range(len(p1)):

			params['kernel_sizes'] = [p0[idx0], p1[idx1]]

			model = Model(alphabet_size, max_word_length, **params)

			model.fit(X['train'], Y['train'])

			y_predicted = model.predict(X['test'])

			q = ~(y_predicted == Y['test'])
			error_rate = 100*q.astype(int).sum()/len(q)

			accuracies[idx0, idx1] = round(100 - error_rate, 1)

	# locate all maxima (and take the first one in case there are multiple)
	idx = numpy.argwhere(accuracies == accuracies.max())[0]

	# adjusting hyperparameters
	#
	params['kernel_sizes'] = [p0[idx[0]], p1[idx[1]]]

	print()
	print('Optimal model values:')
	print('          kernel_sizes:',  params['kernel_sizes'])

#---------------------------------------------------------------------

print()
model.summary()
print()

model.fit(X['train'], Y['train'])

# the trained model can be saved and loaded:
#
#model.save(base + 'model.tf')
#model.load(base + 'model.tf')

# model testing-----------------------------------------------------------

y_predicted = model.predict(X['test'])

q = ~(y_predicted == Y['test'])
error_rate = 100*q.astype(int).sum()/len(q)

print()
print('Overall average precision: {0:.2f}%\n'.format(100 - error_rate))
