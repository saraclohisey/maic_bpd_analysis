import pandas as pd

# Load the input CSV file
input_file = 'Al_mudares_37368978_male_genes_combine.csv'  # Replace filename
df = pd.read_csv(input_file)

# Assuming the columns in the CSV are:
# 'Mouse Gene Ranked' - the correct order of mouse genes
# 'Mouse Gene Unranked' - the order after renaming
# 'Human Ortholog Gene Symbol' - human ortholog genes corresponding to column 2

# Create a dictionary that maps 'Mouse Gene Unranked' to 'Human Ortholog Gene Symbol'
mapping = dict(zip(df['Mouse Gene Unranked'], df['Human Ortholog Gene Symbol']))

# Reorder the human ortholog symbols according to 'Mouse Gene Ranked'
df['Correct Human Ortholog'] = df['Mouse Gene Ranked'].map(mapping)

# Create the output dataframe with the desired columns
output_df = df[['Mouse Gene Ranked', 'Correct Human Ortholog']]

# Rename the columns to match the desired output
output_df.columns = ['Mouse Genes', 'Human Ortholog Symbol']

# Save the reordered data to a new CSV file
output_file = 'Al_mudares_37368978_male_genes_reorder.csv'
output_df.to_csv(output_file, index=False)

print(f"Reordered gene list saved to '{output_file}'.")
