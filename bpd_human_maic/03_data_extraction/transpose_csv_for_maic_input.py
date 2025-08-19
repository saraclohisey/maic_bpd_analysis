import os

import csv

import csv

def transpose_csv_to_tsv(csv_file, tsv_file, placeholder='-'):
    # Read the CSV file
    with open(csv_file, 'r', newline='') as infile:
        reader = list(csv.reader(infile))
        
        # Transpose the data
        transposed_data = zip(*reader)
        
        # Write the transposed data to the TSV file
        with open(tsv_file, 'w', newline='') as outfile:
            tsv_writer = csv.writer(outfile, delimiter='\t')
            
            for row in transposed_data:
                # Remove empty cells from the row
                cleaned_row = [cell for cell in row if cell]
                
                # Write only non-empty rows
                if cleaned_row:
                    tsv_writer.writerow(cleaned_row)

# Example usage
csv_file = '240923_human_maic_input_v1.csv'
tsv_file = '240923_human_maic_input_v1T.csv'
transpose_csv_to_tsv(csv_file, tsv_file)#, placeholder='-')  # You can set your own placeholder here

'''
# Example usage
csv_file = '240923_human_maic_input_v1.csv'
tsv_file = 'output.tsv'
transpose_csv_to_tsv(csv_file, tsv_file)

'''




'''
def read_all_files_in_directory(directory):
    
    # Iterate through all files in the given directory
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        print (file_path)
        
        # Check if it's a file
        #if os.path.isfile(file_path):
        #    with open(file_path, 'r', encoding='utf-8') as f:
        #        file_contents[filename] = f.read()
    
    #return file_contents



directory = '/Users/sclohise/Desktop/03_data_extraction/'

print (read_all_files_in_directory(directory))
'''