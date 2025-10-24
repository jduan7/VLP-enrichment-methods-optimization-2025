import os
import sys

# Function to extract the required part of the filename
def extract_name(file_name):
    # Split the filename at underscores
    parts = file_name.split('_')
    if len(parts) >= 3:
        # Join the parts before the second-to-last underscore
        return '_'.join(parts[:-2])
    return None

# Function to find all R1_001.fastq.gz files and process their names
def find_fastq_files(directory):
    extracted_names = []
    # Walk through the directory
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith("R1_001.fastq.gz"):
                name_part = extract_name(file)
                if name_part:
                    extracted_names.append(name_part)
    return extracted_names

# Main function to execute the script
def main():
    if len(sys.argv) != 2:
        print("Usage: python3 extract_fastq_names.py <directory_path>")
        sys.exit(1)

    directory = sys.argv[1]

    if not os.path.isdir(directory):
        print("Error: Not a valid directory.")
        sys.exit(1)

    extracted_names = find_fastq_files(directory)
    
    # Print the extracted names, one per line
    for name in extracted_names:
        print(name)

if __name__ == "__main__":
    main()

