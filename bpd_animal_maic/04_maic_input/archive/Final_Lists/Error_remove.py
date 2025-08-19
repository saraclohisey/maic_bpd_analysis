import os
import csv

# Folder containing your CSV files
folder_path = '/Users/chris/Dropbox/bpd_animal_maic/04_maic_input/Lists - Gene + HOS Erros Removed'

# Define the terms to look for in the second column
terms_to_remove = ['not found', 'error - retries exceeded']

# Loop through all files in the folder
for filename in os.listdir(folder_path):
    if filename.endswith('.csv'):
        file_path = os.path.join(folder_path, filename)

        # Read the CSV file
        with open(file_path, 'r', newline='', encoding='utf-8') as csvfile:
            reader = csv.reader(csvfile)
            rows = list(reader)

        # Filter out rows where column 2 contains the unwanted terms
        cleaned_rows = [row for row in rows if not any(term in row[1].lower() for term in terms_to_remove)]

        # Write the cleaned data back to the same file
        with open(file_path, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(cleaned_rows)

print("Processing complete. Rows with 'not found' or 'error - retries exceeded' in column 2 have been removed.")
