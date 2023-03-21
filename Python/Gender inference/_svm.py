import os

import numpy
from sklearn import svm, metrics
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

#import tensorflow.keras as kr

#print(kr.models.Sequential())

# path to the current directory
base = os.path.dirname(os.path.realpath(__file__)) + os.sep

""" path to the 'name_gender' table """
data_path = base + "name_gender.csv"

#-------------------------------------------------------------------------

alphabet_size = 26			# assuming English alphabet by default
							# (may be adjusted during pre-processing)

max_word_length = 18

# input features dimensionality
#
#   for names containing only ASCII characters is equal to
#   `alphabet_size*max_word_length`
#
#   the maximum between the initial value of `input_dim`
#   and the one determined from the datatset will be used
#
input_dim = alphabet_size*max_word_length

# converts a `word` (represented by an array of integers)
# into a flattened array of 1-hot vectors
#
def array2vec(word, alphabet_size, max_word_length):

	n = len(word)

	if n > max_word_length:

		print('The word is too long. Truncating')
		n = max_word_length

	v = numpy.zeros(alphabet_size*max_word_length)

	for idx in range(n):

		v[alphabet_size*idx + word[idx]] = 1

	return v

# data acquisition and pre-processing-------------------------------------

lines = None

X = {'train': None, 'test': None}
Y = {'train': None, 'test': None}

with open(data_path, 'r', encoding = "utf-8") as f:

	lines = f.readlines()

if lines is not None:

	#----------------------------------------------------

	def sanitize_and_split(s):

		p = s.upper().strip(' \r\n').split(',')
		p[0] = ''.join((ch for ch in p[0] if ch.isalpha()))

		return p

	#----------------------------------------------------

	# convert lines to an array of pairs [name, gender]
	data = map(sanitize_and_split, lines[1:])

	# 'unzip' data (convert from array of pairs to a pair of arrays)
	data = [a for a in zip(*data)]

	# convert to a bytearray and shift all the codes by 'A'
	# (so that [A..Z] is converted to [0..25])
	#
	A_code = bytearray('A', "utf-8")[0]
	data[0] = list(map(lambda s: numpy.array(bytearray(s, "utf-8")) - A_code, data[0]))

	# [optional] adjusting the alphabet size
	#
	#     (should be used if non-ASCII characters may be present in data)
	#
	if True:
		alphabet = set(numpy.hstack(data[0]))
		alphabet_size = max(alphabet_size, max(alphabet) + 1)

	# adjusting `input_dim`

	data_max_word_length = max(map(len, data[0]))
	max_word_length = max(max_word_length, data_max_word_length)

	input_dim = max(input_dim, alphabet_size*max_word_length)

	# converting arrays of integers to flattened 1-hot representation
	data[0] =  list\
	(
		map(lambda word: array2vec(word, alphabet_size, max_word_length), data[0])
	)

	# [optional] printing the dataset stats
	#
	if True:
	
		print()
		print('Dataset stats:')
		print('   Name count      =', len(data[0]))
		print('   Max word length =', data_max_word_length)
		print('   Alphabet size   =', alphabet_size)

	# creating training and testing datasets

	x = numpy.array(data[0])		# input features
	y = numpy.array(data[1])		# labels

	#----------------------------------------------------

	# splits an array `a` into a pair of (training, testing) arrays
	#
	def dataset_split(a, train_fraction = 0.9):

		idx = round(len(a)*train_fraction)
		return (a[:idx], a[idx:])

	#----------------------------------------------------

	X['train'], X['test'] = dataset_split(x)
	Y['train'], Y['test'] = dataset_split(y)

else:
	print("Error reading data from file")
	exit(0)

# model fitting-----------------------------------------------------------

def Model(C, max_iter = -1):

	return make_pipeline \
	(
		StandardScaler(),
		#svm.SVC(kernel = 'sigmoid', gamma = gamma, C = C, coef0 = r, max_iter = max_iter)
		#svm.SVC(kernel = 'rbf', gamma = gamma, C = C, max_iter = max_iter)
		svm.LinearSVC(dual = False, C = C)
	)

# hyperparameters
#
gamma = 'auto'
C = 1.0

# hyperparameter tuning---------------------------------------------------
#
if False:

	gammas = [1.0]
	Cs = [1.0, 10., 100., 1000., 10000.]

	error_rates = numpy.zeros((len(gammas), len(Cs)))

	for idx_0, gamma in enumerate(gammas):
		for idx_1, C in enumerate(Cs):

			model = Model(C = C, max_iter = 20)
			model.fit(X['train'], Y['train'])

			y_predicted = model.predict(X['test'])

			q = ~(y_predicted == Y['test'])
			error_rate = 100*q.astype(int).sum()/len(q)

			error_rates[idx_0, idx_1] = error_rate

	print()
	print(error_rates)

# model fitting-----------------------------------------------------------

model = Model(C)

model.fit(X['train'], Y['train'])

#print(model['linearsvc'].coef_.shape)

# model testing-----------------------------------------------------------

y_predicted = model.predict(X['test'])

q = ~(y_predicted == Y['test'])
error_rate = 100*q.astype(int).sum()/len(q)

print()
print('Overall average precision: {0:.2f}%\n'.format(100 - error_rate))
#print(metrics.classification_report(Y['test'], y_predicted))

disp = metrics.ConfusionMatrixDisplay.from_predictions(Y['test'], y_predicted)
print(f"Confusion matrix:\n{disp.confusion_matrix}")
