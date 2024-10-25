#!/bin/bash  

## parallel_config.bash
#
# A script to configure a FireFOAM case to run in parallel
#
# No modifications to the script should be needed


# Function to display usage information  
usage() {  
    echo "Recursively modify the 'numberOfSubdomains' and 'SBATCH' lines in all files in the present directory."
    echo "Usage: $0 --num-nodes=<num_nodes> --cores-per-node=<cores_per_node>"  
    echo "Options:"  
    echo "  --num-nodes=<num_nodes>       Number of nodes to use"  
    echo "  --cores-per-node=<cores_per_node> Number of cores per node"  
    echo "  --help                        Display this help message and exit"  
    exit 1  
}  
  
# Initialize variables  
num_nodes=""  
cores_per_node=""  
  
# Parse command line arguments  
while [ "$#" -gt 0 ]; do  
    case "$1" in  
        --num-nodes=*)  
            num_nodes="${1#*=}"  
            shift  
            ;;  
        --cores-per-node=*)  
            cores_per_node="${1#*=}"  
            shift  
            ;;  
        --help)  
            usage  
            ;;  
        *)  
            echo "Unknown parameter: $1"  
            usage  
            ;;  
    esac  
done  
  
# Check if the required arguments are provided  
if [ -z "$num_nodes" ] || [ -z "$cores_per_node" ]; then  
    echo "Error: Both --num-nodes and --cores-per-node must be provided."  
    usage  
fi  
  
# Calculate the total number of tasks (subdomains)  
total_tasks=$((cores_per_node * num_nodes))  
  
# Find and modify the relevant lines in all files containing SBATCH commands  
echo "Updating SBATCH lines in all relevant files..."  
find . -type f -print0 | xargs -0 grep -l "^#SBATCH" | while read -r file; do  
    sed -i "s/^#SBATCH --nodes=.*/#SBATCH --nodes=${num_nodes}/" "$file"  
    sed -i "s/^#SBATCH --ntasks=.*/#SBATCH --ntasks=${total_tasks}/" "$file"  
    sed -i "s/^#SBATCH --ntasks-per-node=.*/#SBATCH --ntasks-per-node=${cores_per_node}/" "$file"  
    echo "Modified: $file"  
done  
  
echo "SBATCH lines have been updated in relevant files."  
  
# Update numberOfSubdomains in input files  
echo "Updating numberOfSubdomains to ${total_tasks}"  
find . -type f -exec sh -c '  
    new_value="$1"  
    shift  
    for file in "$@"; do  
        if grep -q "numberOfSubdomains [0-9]\+" "$file"; then  
            sed -i "s/numberOfSubdomains [0-9]\+/numberOfSubdomains $new_value/" "$file"  
            echo "Modified: $file"  
        fi  
    done  
' sh "${total_tasks}" {} +  
  
echo "numberOfSubdomains has been updated in input files."  

