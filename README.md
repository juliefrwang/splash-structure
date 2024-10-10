# Splash-structure

## 1. Setting up the Python Environment

First, ensure you have Python `>= 3.9` installed. Then, follow these steps:

```bash
# Clone the repository
git clone https://github.com/juliefrwang/splash-structure.git
cd splash-structure

# Create a virtual environment (you can replace 'my_env' with any name)
python3 -m venv my_env

# Activate the virtual environment
# For Linux/Mac:
source my_env/bin/activate

# Install the required dependencies
pip install -r requirements.txt
```
## 2. Installing Julia programming language

Make sure Julia is installed on your system. You can download and install Julia from the official website: https://julialang.org/downloads/.

After installing, verify your Julia installation:
```bash
julia --version
```
To install the required Julia packages for this project, run:
```bash
julia -e 'using Pkg; Pkg.add("DataFrames"); Pkg.add("CSV"); Pkg.add("Combinatorics")'
```
__Note__: This method for installing Julia packages is tentative and will be updated in the future.

## 3. Running SPLASH-strucure on test data

File `test_runs/Test.sh` contains the commands to run SPLASH-structure on the test data, under both target mode and compactor mode. To run the test data, execute the following command:
```bash
# Make sure to navigate to the test_runs directory
cd test_runs

# Run the test script 
./Test.sh
```
