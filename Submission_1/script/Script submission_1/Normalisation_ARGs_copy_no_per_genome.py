import pandas as pd
import glob
import os
import re

# 1. Metadata for Total Reads sample ID with respective total reads mapped
total_reads_meta = {
    '34': 22294672, '37': 23245086, '38': 22734902, '40': 26687198, '41': 23072510,
    '42': 24207120, '44': 24724472, '45': 26918498, '46': 23765148, '47': 25845720,
    '48': 24391854, '49': 26294204, '51': 26810354, '52': 27661712, '53': 31761902,
    '54': 24776086, '55': 22507018, '56': 26452462, '57': 23918566, '58': 27064976,
    '59': 21468314, '60': 25338294, '61': 28272880, '62': 25410692, '63': 23060554,
    '64': 25679342, '65': 27890904, '66': 22558358, '67': 23540280, '68': 23621048,
    '70': 22238620, '71': 23928220, '72': 30386846, '74': 26159118, '75': 26794054,
    '77': 20371482, '78': 25395848, '80': 23045398, '81': 26007938, '82': 24570802,
    '84': 22432514, '85': 22938332, '86': 29284616, '88': 22328310, '89': 25704536,
    '90': 25001378, '91': 21837176, '94': 23696378, '97': 30853668, '98': 25154708,
    '99': 39184056, '100': 27281504, '101': 27252898, '103': 32362714, '33': 24027406,
    '35': 23421868, '36': 21046116, '39': 22100794, '43': 22760350, '69': 23091460,
    '73': 24820118, '76': 21925516, '79': 23819324, '83': 22365186, '87': 18711424,
    '92': 24096624, '93': 27716560, '96': 38407150, '102': 24664514, '50': 25530932,
    '32': 20042588, '95': 19381202,
}

# 2. Universal Parameters
AVG_GENOME_SIZE = 4.64e6
READ_LENGTH = 150


def run_normalization(input_dir, output_file):
    all_summaries = []

    # 3. Use glob to find files in given folder
    file_pattern = os.path.join(input_dir, "*_output_allele_map.csv")
    files = glob.glob(file_pattern)

    print(f"Found {len(files)} files to process...")

    for filepath in files:
        filename = os.path.basename(filepath)

        # Regex to extract ID from file name
        match = re.search(r'Sample-(\d+)', filename)
        if not match:
            continue
        sample_id = match.group(1)

        if sample_id not in total_reads_meta:
            print(f"Skipping Sample {sample_id}: Total reads not in metadata.")
            continue

        # Formula logic
        total_reads = total_reads_meta[sample_id]
        B = total_reads * READ_LENGTH
        genome_depth = B / AVG_GENOME_SIZE

        # Load file
        try:
            df = pd.read_csv(filepath)
            # Copy per genome = Depth (Reads/Length) / Genome Depth
            df['Normalized_Value'] = df['Depth'] / genome_depth

            # Aggregate by Family
            summary = df.groupby('AMR Gene Family')['Normalized_Value'].sum().reset_index()
            summary['Sample ID'] = sample_id
            all_summaries.append(summary)
            print(f"Successfully processed Sample {sample_id}")
        except Exception as e:
            print(f"Error processing {filename}: {e}")

    if not all_summaries:
        print("Error: No data was processed. Check your paths and filenames.")
        return

    # 4. Create the final matrix
    combined_df = pd.concat(all_summaries)
    final_matrix = combined_df.pivot(index='Sample ID', columns='AMR Gene Family', values='Normalized_Value')

    # Fill NaN with 0 and clean up
    final_matrix = final_matrix.fillna(0).reset_index()

    # Save the file
    final_matrix.to_csv(output_file, index=False)
    print(f"\nDONE! Final file saved at: {output_file}")


#EXECUTION STARTS HERE
# Define paths
input_folder = "/path/to/folder/"
output_path = "/path/to/folder/"

# CALL THE FUNCTION
run_normalization(input_folder, output_path)
