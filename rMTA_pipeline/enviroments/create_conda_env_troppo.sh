#!/bin/bash

# Name of the Conda environment
ENV_NAME="troppo_env"

# Python version
PYTHON_VERSION="3.7"

# List of dependencies to install via Conda (available in conda-forge)
CONDA_DEPENDENCIES=(
cobra
#cmake
#swig
#swiglpk
)

# List of dependencies to install via pip (not available in Conda)
PIP_DEPENDENCIES=(

troppo
psutil
#python-libsbml
#cobra

)

# Function to check if a Conda environment exists
function conda_env_exists {
    conda env list | grep -q "^${ENV_NAME}[[:space:]]"
}

# Check if the environment exists
if conda_env_exists; then
    echo "Environment '${ENV_NAME}' exists. Deleting..."
    conda env remove -n "${ENV_NAME}" -y
    if [ $? -ne 0 ]; then
        echo "Failed to remove environment '${ENV_NAME}'. Exiting."
        exit 1
    fi
else
    echo "Environment '${ENV_NAME}' does not exist. Proceeding to create it."
fi

# Create the new environment with specified Python version and dependencies via conda-forge
echo "Creating environment '${ENV_NAME}' with Python ${PYTHON_VERSION} and specified dependencies..."
conda create -y -n "${ENV_NAME}" -c conda-forge python="${PYTHON_VERSION}" "${CONDA_DEPENDENCIES[@]}"
if [ $? -eq 0 ]; then
    echo "Conda environment '${ENV_NAME}' created successfully."
else
    echo "Failed to create environment '${ENV_NAME}'. Exiting."
    exit 1
fi

# Activate the environment
echo "Activating environment '${ENV_NAME}'..."
# Initialize Conda for the script
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"
if [ $? -ne 0 ]; then
    echo "Failed to activate environment '${ENV_NAME}'. Exiting."
    exit 1
fi
echo "Environment '${ENV_NAME}' is now active."

# Install pip dependencies
if [ ${#PIP_DEPENDENCIES[@]} -gt 0 ]; then
  echo "Installing pip dependencies..."
  pip install "${PIP_DEPENDENCIES[@]}"
  #pip install --no-cache-dir --force-reinstall "${PIP_DEPENDENCIES[@]}"

  if [ $? -eq 0 ]; then
      echo "Pip dependencies installed successfully."
  else
      echo "Failed to install pip dependencies. Exiting."
      exit 1
  fi

else
    echo "No pip dependencies to install."

fi

# Install your package in editable mode
#echo "Installing your package in editable mode..."
#pip install -e .
#if [ $? -eq 0 ]; then
#    echo "Package installed successfully in environment '${ENV_NAME}'."
#else
#    echo "Failed to install package in environment '${ENV_NAME}'. Exiting."
#    exit 1
#fi

# Inform the user
echo "Conda environment '${ENV_NAME}' is set up and your package is installed."
