import pandas as pd

# Load the CSV file
input_file = 'Rao_37310763_ranked.csv'  # Change this to your input file name
output_file = 'Rao_37310763_ranked_duplicates_removed.csv'  # Output file for unique genes

try:
    # Read the CSV file
    df = pd.read_csv(input_file)

    # Check if 'Mouse Genes' column exists
    if 'Mouse Genes' not in df.columns:
        raise ValueError("The column 'Mouse Genes' does not exist in the file.")

    # Remove duplicates and keep the first occurrence
    df_unique = df['Mouse Genes'].drop_duplicates()

    # Save the unique genes to a new CSV file
    df_unique.to_csv(output_file, index=False)

    print(f"Duplicates removed. Unique genes saved to '{output_file}'.")

except FileNotFoundError:
    print(f"Error: The file '{input_file}' was not found.")
except Exception as e:
    print(f"An error occurred: {e}")
