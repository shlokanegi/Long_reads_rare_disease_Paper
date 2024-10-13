import os
import json
import sys

# Function to update output JSON file
def update_output(output_json_filename, stats):
    if os.path.exists(output_json_filename):
        with open(output_json_filename, 'r') as file:
            output_data = json.load(file)
    else:
        output_data = {}

    output_data.update(stats)

    with open(output_json_filename, 'w') as file:
        json.dump(output_data, file, indent=4)

# Function to process sample directories
def process_sample_directories(intersection_tsv, lrsonly_tsv, srsonly_tsv, output_json_filename):
    stats = {}
    with open(intersection_tsv) as fi, open(lrsonly_tsv) as fl, open(srsonly_tsv) as fs:
        for line_i, line_l, line_s in zip(fi, fl, fs):
            sample, i_snps, i_indels = line_i.strip().split()[0], line_i.strip().split()[1], line_i.strip().split()[2]
            sample, l_snps, l_indels = line_l.strip().split()[0], line_l.strip().split()[1], line_l.strip().split()[2]
            sample, s_snps, s_indels = line_s.strip().split()[0], line_s.strip().split()[1], line_s.strip().split()[2]

            if sample not in stats:
                stats[sample] = {
                    "SNP": {
                        "Intersection": i_snps,
                        "LRS-only": l_snps,
                        "SRS-only": s_snps
                    },
                    "INDEL": {
                        "Intersection": i_indels,
                        "LRS-only": l_indels,
                        "SRS-only": s_indels
                    }
                }
            else:
                # Update existing stats if the sample already exists
                stats[sample]["SNP"]["Intersection"] = i_snps
                stats[sample]["SNP"]["LRS-only"] = l_snps
                stats[sample]["SNP"]["SRS-only"] = s_snps
                stats[sample]["INDEL"]["Intersection"] = i_indels
                stats[sample]["INDEL"]["LRS-only"] = l_indels
                stats[sample]["INDEL"]["SRS-only"] = s_indels

    update_output(output_json_filename, stats)


if __name__ == "__main__":
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 5:
        print("Usage: python3 generate-intersection-json.py <intersection_tsv> <lrsonly_tsv> <srsonly_tsv> <output_json_filename>")
        sys.exit(1)

    intersection_tsv = sys.argv[1]
    lrsonly_tsv = sys.argv[2]
    srsonly_tsv = sys.argv[3]
    output_json_filename = sys.argv[4]

    intersection_tsv, lrsonly_tsv, srsonly_tsv

    process_sample_directories(intersection_tsv, lrsonly_tsv, srsonly_tsv, output_json_filename)