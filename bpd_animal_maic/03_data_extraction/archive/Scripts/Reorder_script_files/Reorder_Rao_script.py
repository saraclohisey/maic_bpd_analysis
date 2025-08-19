import pandas as pd

# Load the CSV file (replace 'your_file.csv' with the path to your actual file)
df = pd.read_csv('Reorder_Rao.csv')

# Initialize a list to store the ordered rows for the new DataFrame
ordered_data = []

# Loop through each gene in "Mouse Genes Ranked" and find the corresponding human ortholog
for gene in df['Mouse Genes Ranked']:
    # Check if the gene exists in the 'Mouse Genes Unranked' column
    match = df[df['Mouse Genes Unranked'] == gene]
    
    if not match.empty:
        # If a match is found, get the corresponding human ortholog symbol
        human_symbol = match['Human Ortholog Symbol'].values[0]
    else:
        # If no match is found, leave human symbol as None
        human_symbol = None
    
    # Append the ordered data
    ordered_data.append({'Mouse Genes': gene, 'Human Ortholog Symbol': human_symbol})

# Create a DataFrame from the ordered data
df_ordered = pd.DataFrame(ordered_data)

# Save the ordered data to a new CSV file (replace 'ordered_file.csv' with your desired output file name)
df_ordered.to_csv('ordered_file.csv', index=False)

# Corrected print statement
print("Ordered file created and saved as 'ordered_file.csv'")
