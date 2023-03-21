""" Suffix-based name classification model """

class Model:

	"""

	One of the main problem with word classification is the fact that words should be
	embedded (implicitly or explicitly) into some vector space beforehand. Or, equivalently,
	a meaningful distance function between words should be defined, which has proven to be
	quite a challenge. word2vec may serve as one of the examples of obtaining such a function,
	the core idea behind the model being "two words are the more similar the more they are used
	in the same context". This idea is rather intuitive and, given a large enough text corpus,
	yields remarkable results in terms of obtaining a meaningful representation of words in a 
	metric space.

	Name, however, is a rather special case of a word. In terms of 'contextual similarity' all
	the names are used in very similar contexts (especially in English which doesn't have gender-specific
	verb and adjective endings). One may argue though that names originated from meaningful
	words, however there seems to be no obvious correlation between a words meaning and its gender.
	Moreover, the same concept may have different genders in different languages (for example, the
	German word 'Mädchen' ('young women') has a neutral gender while virtually the same word in Spanish
	'niña' has feminine gender). Thus, it appears that the gender of each given name is essentially arbitrary
	and the only way to categorize names is to use a dictionary.

	Still, sometimes it seems like there is a certain pattern when it comes to somewhat 'common' names:
	even if you have never heard a name before you can sometimes 'feel' what gender is it likely to be.
	For example, I personally have never heard neither 'Aadhya' nor 'Aadhvik', but I feel like the first
	should be female and the second one should be male. On the contrary, it feels to me that 'Aadhyan'
	should be female, while it is apparently male. So, this apparent 'feel' is very elusive, but upon
	further investigation I came up with a heuristic which allowed me to formalize it.

	The main idea behind the name classification model is to use suffix matching, i.e. to try to estimate
	the gender by only looking at 1, 2, 3 etc. last letters of the name. For example, all the names 
	'Derielle', 'Elizabelle', 'Elle', 'Ercelle' and 'Gabrielle' are female, and their common feature
	is the suffix 'elle'. There are male names with the same suffix however ('Danzelle', 'Darrelle', 'Jamelle')
	but the ratio of the number of female to male names with this suffix is about 63 to 7, which means that about
	90% of the names ending with 'elle' are female. To improve the classification accuracy one may look at longer
	suffixes, though at the cost of the model size.

	The general algorithm for training the model is as follows: starting with suffixes of length 1 up to some
	maximal value check how well does each suffix categorizes the names. The categorization score used is the
	so called 'decisiveness', which is defined as the difference between the number of male and female names
	with the suffix in question divided by the total number of the names with this suffix. The idea is that when
	the counts for each gender are close that means that the suffix does not help to make a decision about the
	gender, so this suffix should not be considered during the prediction phase. After training the model contains
	a map (dict) between suffixes and the corresponding label and 'decisiveness'.

	The algorithm for prediction the gender using the trained model is as follows: starting with longer suffixes
	check whether the analyzed name has this suffix and in case it does append the score array with its label and
	score. The idea is to go from more specific suffixes to less and less specific, ending with one-letter suffixes
	(which still provide some information about the name in case all the other suffixes have failed). If none of
	the suffixes retained by the model are found in the name the predicted label is `None`.

	"""

	# converts a `word` to a 'standard form':
	#
	#   - converts to upper case
	#   - removes all non-alphabetical characters
	#
	@staticmethod
	def normalize(word):

		return ''.join((ch for ch in word.upper() if ch.isalpha()))

	#-----------------------------------------------------------------

	def __init__(self):

		self.suffixes = {}
		self.max_suffix_length = 0

	# threshold: suffix 'decisiveness' threshold
	#
	def fit(self, words, labels, max_suffix_length = 4, threshold = 0.9):

		self.suffixes.clear()
		self.max_suffix_length = max_suffix_length

		label_list = list(set(labels))
		initial_dict = lambda: dict(zip(label_list, [0]*len(label_list)))

		for suffix_length in range(1, max_suffix_length+1):
		
			# calculating suffix support
			#
			for idx in range(len(words)):

				word = self.normalize(words[idx])

				if suffix_length <= len(word):

					suffix = word[-suffix_length:]
					label = labels[idx]

					if suffix not in self.suffixes:
						self.suffixes[suffix] = initial_dict()

					self.suffixes[suffix][label] += 1

			# cutting the least 'decisive' suffixes

			suffixes = [key for key in self.suffixes]

			for suffix in suffixes:

				if len(suffix) == suffix_length:

					support = self.suffixes[suffix]

					vals = sorted(support.values())
					min_diff = min([q - p for p, q in zip(vals, vals[1:])])

					# relative 'decisiveness'
					#
					#  the bigger the value the more effective the current suffix is for classification
					#
					d = min_diff/sum(vals)

					if d < threshold:

						self.suffixes.pop(suffix)
					else:

						# extract the pair (label, support) with the maximal support
						m = max(support.items(), key = lambda p: p[1])

						self.suffixes[suffix] = (m[0], d)		# (label, decisiveness)

	def predict(self, words):

		labels = [None]*len(words)

		for idx in range(len(words)):

			word = self.normalize(words[idx])

			max_suffix_length = min(len(word), self.max_suffix_length)
			cumulative_scores = {}

			for suffix_length in range(max_suffix_length, 0, -1):

				suffix = word[-suffix_length:]

				if suffix in self.suffixes:

					s_label, s_score = self.suffixes[suffix][0:2]

					if s_label not in cumulative_scores:
						cumulative_scores[s_label] = 0

					cumulative_scores[s_label] += s_score

			if cumulative_scores:

				# after all the suffixes have been tested `cumulative_scores` contains
				# cumulative scores for each of the possible labels (except when this score is 0);
				#
				# the final label is the one that corresponds to the maximum cumulative score:
				#
				labels[idx] = max(cumulative_scores.items(), key = lambda p: p[1])[0]

		return labels

	def save(self, path):
	
		with open(path, 'w', encoding = 'utf-8') as f:

			f.write \
			(
				str \
				(
					{
						'suffixes': self.suffixes,
						'max_suffix_length': self.max_suffix_length
					}
				)
			)

	def load(self, path):

		with open(path, 'r', encoding = 'utf-8') as f:

			d = eval(f.read())

			self.suffixes = d['suffixes']
			self.max_suffix_length = d['max_suffix_length']
