# Get the directory of the script
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
MAIN_FILE=../../build/main_json #TODO: make this more robust

# cd to the script directory
cd $SCRIPT_DIR

# Create the config files
mkdir -p inputs
rm -f inputs/*
python3 make_configs.py

# Remove the old output files
mkdir -p outputs
rm -f outputs/*

# Log some information about the run environment
echo "Running on $(hostname)"
echo "Running on $(date)"

# Run all the config files
for f in inputs/*.json; do
    $MAIN_FILE $f ${f/inputs/outputs}
done

# Make the plots
python3 make_plots.py