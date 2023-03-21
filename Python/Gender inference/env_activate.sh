#!/bin/sh

# paths (please change accordingly)

# path to the virtual environment containing tensorflow
VENV_PATH=

# path to CUDA Toolkit
CUDA_PATH=

# path to cuDNN
cuDNN_PATH=

#----------------------------------------------------------------------

# activate the virtual environment containing tensorflow
$VENV_PATH/activate

# change path to the current directory
cd $PWD

# adding paths to CUDA toolkit and cuDNN
PATH=$CUDA_PATH/bin:$PATH
PATH=$CUDA_PATH/extras/CUPTI/lib64:$PATH
PATH=$CUDA_PATH/include:$PATH
PATH=$cuDNN_PATH/bin:$PATH