@echo off

:: paths (please change accordingly)

:: path to the virtual environment containing tensorflow
set VENV_PATH=C:\Files\Python\_venv\TF\Scripts

:: path to CUDA Toolkit
set CUDA_PATH=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.2

:: path to cuDNN
set cuDNN_PATH=C:\Files\cuDNN\CUDA 11.2

::----------------------------------------------------------------------

:: activate the virtual environment containing tensorflow
%VENV_PATH%/activate

:: change path to the current directory
cd %~dp0

:: adding paths to CUDA toolkit and cuDNN
set PATH=%CUDA_PATH%/bin;%PATH%
set PATH=%CUDA_PATH%/extras/CUPTI/lib64;%PATH%
set PATH=%CUDA_PATH%/include;%PATH%
set PATH=%cuDNN_PATH%/bin;%PATH%