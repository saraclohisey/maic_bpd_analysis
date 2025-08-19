import pandas as pd

# Load the CSV file
input_file = 'Rao_37310763_ranked_HOS.csv'  # Change this to your input file name
output_file = 'Rao_37310763_ranked_HOS_duplicates_removed.csv'  # Output file for unique genes

try:
    # Read the CSV file
    df = pd.read_csv(input_file)

    # Check if required columns exist
    required_columns = ['Mouse Genes Unranked', 'Human Ortholog symbol']
    if not all(col in df.columns for col in required_columns):
        raise ValueError("The required columns are not present in the file.")

    # Remove duplicates based on 'Mouse Genes Unranked' while keeping the first occurrence
    df_unique = df.drop_duplicates(subset='Mouse Genes Unranked', keep='first')

    # Save the unique rows to a new CSV file
    df_unique.to_csv(output_file, index=False)

    print(f"Duplicates removed. Unique genes saved to '{output_file}'.")

except FileNotFoundError:
    print(f"Error: The file '{input_file}' was not found.")
except Exception as e:
    print(f"An error occurred: {e}")
