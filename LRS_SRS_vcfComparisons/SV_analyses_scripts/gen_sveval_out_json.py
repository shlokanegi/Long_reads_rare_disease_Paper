import os
import json
import sys
import pandas as pd

# Function to update output JSON file
def update_output(output_file, stats):
    if os.path.exists(output_file):
        with open(output_file, 'r') as file:
            output_data = json.load(file)
    else:
        output_data = {}

    output_data.update(stats)

    with open(output_file, 'w') as file:
        json.dump(output_data, file, indent=4)
        print("done processing output file")

# Function to process sample directories
def process_sample_files(directory, output_file, algorithm):
    for root, dirs, files in os.walk(directory):
        for file in files:
            sveval_out_file = os.path.join(directory, file)  #E.g. /private/groups/migalab/shnegi/raredis/SV_analysis/broad_SVs/sveval/sveval_out/M11AO_sveval_out.tsv
            if os.path.exists(sveval_out_file):
                df = pd.read_csv(sveval_out_file, sep='\t', header=0)
                print(df)
                # extract values
                tp = int(df.loc[df.type=="Total", "TP"].iloc[0])
                fp = int(df.loc[df.type=="Total", "FP"].iloc[0])
                fn = int(df.loc[df.type=="Total", "FN"].iloc[0])
                tp_base = int(df.loc[df.type=="Total", "TP.baseline"].iloc[0])
                sample = file.split('_')[0]
                lrs_only = f'{algorithm}-only'
                # Construct the stats dictionary
                stats = {
                    sample: {
                        "Intersection":tp,
                        lrs_only:fn,
                        "SRS-only":fp
                    }
                }
                print(stats)
                update_output(output_file, stats)

if __name__ == "__main__":
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_directory_with_sveval_outputs> <output_json_filename> <LRS_algorithm>")
        sys.exit(1)

    input_directory = sys.argv[1]
    output_file = sys.argv[2]
    algorithm = sys.argv[3]

    process_sample_files(input_directory, output_file, algorithm)