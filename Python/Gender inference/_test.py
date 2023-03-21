import os

import numpy

import model_suffix
import model_cnn

#--------------------------------------------------------------------

# suppress tensorflow warnings
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

base = os.path.dirname(os.path.realpath(__file__)) + os.sep

s_model = model_suffix.Model()
s_model.load(base + 'model.dict')

cnn_model = model_cnn.Model()
cnn_model.load(base + 'model.tf')

#--------------------------------------------------------------------

# none of the names below can be found in 'name_gender.csv'
#
names = ['Qventin', 'Lollita', 'Sero', 'Bertold', 'Geh', 'Lloid', 'Qvar', 'Llavar', 'Hurrem', 'Gulfem']

s_labels   = s_model.predict(names)
cnn_labels = cnn_model.predict(names)

print()
print('Names outside of the provided dataset and/or nonsense names:\n')
print('suffix-based model:\n', numpy.array(list(zip(names, s_labels))))
print()
print('CNN-based model:\n', numpy.array(list(zip(names, cnn_labels))))
print()
