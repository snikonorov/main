import os

import tensorflow.keras as keras
from helpers import *

#------------------------------------------------------------------------

""" CNN-based name classification model """

class Model:

	"""

	The main idea of this model is to represent each word as a matrix (essentially, a stack of one-hot vectors
	each of which encodes a character) and then use 2D convolutions to both learn 'contextual' information (in 
	terms of character neighborhood) and to extract more meaningful and compact representation of the initial 
	one-hot character encoding.

	"""

	# converts a `word` to a 'standard form':
	#
	#   - converts to upper case
	#   - removes all non-alphabetical characters
	#
	@staticmethod
	def normalize(word):

		return ''.join((ch for ch in word.upper() if ch.isalpha()))

	# converts a string `s` into an array of 1-hot vectors
	# (transforms it into a bytearray and shifts the codes first)
	#
	#    [return].shape == (alphabet_size, max_word_length)
	#
	@staticmethod
	def string2vec(s, alphabet_size, max_word_length):

		A_code = 65
		word = numpy.array(bytearray(s, "utf-8")) - A_code

		n = len(word)

		if n > max_word_length:

			print('The string is too long. Truncating.')
			n = max_word_length

		v = numpy.zeros((alphabet_size, max_word_length))

		for idx in range(n):

			v[word[idx], idx] = 1

		return v

	#-----------------------------------------------------------------

	def __init__ \
	(
		self, alphabet_size = None, max_word_length = None, *,
		n_filters = 64,
		kernel_sizes = [(9, 7), (9, 7)],
		dense_dim = 100
	):

		if alphabet_size is not None and max_word_length is not None:

			self.m = keras.models.Sequential \
			([
				keras.Input(shape = (alphabet_size, max_word_length, 1)),
				keras.layers.Conv2D(n_filters, kernel_sizes[0], activation = 'relu'),
				keras.layers.Conv2D(n_filters, kernel_sizes[1], activation = 'relu'),
				keras.layers.Flatten(),
				keras.layers.Dense(dense_dim, activation = 'relu'),
				keras.layers.Dense(1, activation = 'sigmoid')
			])

			self.m.build((alphabet_size, max_word_length, 1))

			self.m.compile \
			(
				loss = "binary_crossentropy",
				optimizer = keras.optimizers.Adam(learning_rate = 1e-3),
				metrics = ["accuracy"]
			)

		else:
			self.m = None

		self.alphabet_size = alphabet_size
		self.max_word_length = max_word_length

		self.number_to_label = None
		self.label_to_number = None

	# convert `words` into input features
	#
	def words_to_features(self, words):

		return numpy.array \
		(
			[
				self.string2vec(self.normalize(word), self.alphabet_size, self.max_word_length)
				for word in words
			]
		)

	def fit(self, words, labels, epochs = 10, batch_size = 64, validation_split = 0.1):

		# converting words and labels into matrices and numbers

		x_train = self.words_to_features(words)

		self.number_to_label = dict(enumerate(set(labels)))
		self.label_to_number = inv_dict(self.number_to_label)

		y_train = numpy.array \
		(
			[self.label_to_number[label] for label in labels]
		)
		
		# x_train.shape == (N, alphabet_size, max_word_length)
		# y_train.shape == (N, )

		self.m.fit \
		(
			x_train[:, :, :, None], y_train[:, None],
			batch_size = batch_size,
			epochs = epochs,
			validation_split = validation_split
		)

	def predict(self, words):

		x = self.words_to_features(words)

		# x.shape        == (N, alphabet_size, max_word_length)
		# [return].shape == (N, )

		y = numpy.round(self.m.predict(x[:, :, :, None])[:, 0]).astype(int)

		return numpy.array \
		(
			[self.number_to_label[a] for a in y]
		)

	def summary(self):

		print('CNN-based classifier')
		print()
		print \
		(
			'input.shape: (N: <any>, alphabet_size: {}, max_word_length: {})\n'
			.format(self.alphabet_size, self.max_word_length)
		)

		if self.m is not None:
			self.m.summary()

	def save(self, path):

		self.m.save(path)

		params = \
		{
			'alphabet_size':    self.alphabet_size,
			'max_word_length':  self.max_word_length,
			'number_to_label':  self.number_to_label,
			'label_to_number':  self.label_to_number
		}

		with open(path + os.sep + '_params', 'w') as f:
			f.write(str(params))

	def load(self, path):

		self.m = keras.models.load_model(path)

		with open(path + os.sep + '_params', 'r') as f:

			params = eval(f.read())

			for attr in params:
				setattr(self, attr, params[attr])
