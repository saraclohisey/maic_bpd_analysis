import pandas as pd

# Load the CSV file (replace 'your_file.csv' with the path to your actual file)
df = pd.read_csv('Rao_unrankedHOS_duplicates.csv')

# Ensure column names are correct by stripping any extra whitespace
df.columns = df.columns.str.strip()

# Separate rows with and without an 'Identity' value
df_with_identity = df.dropna(subset=['Identity'])
df_without_identity = df[df['Identity'].isna()]

# For rows with an 'Identity' value, sort by 'Identity' and remove duplicates by keeping the highest 'Identity'
df_with_identity_cleaned = df_with_identity.sort_values(by='Identity', ascending=False).drop_duplicates(subset='Mouse Genes', keep='first')

# For rows without an 'Identity' value, keep only the first occurrence of each duplicate
df_without_identity_cleaned = df_without_identity.drop_duplicates(subset='Mouse Genes', keep='first')

# Combine the two DataFrames
df_final = pd.concat([df_with_identity_cleaned, df_without_identity_cleaned])

# Save the cleaned data to a new CSV file (replace 'cleaned_file.csv' with your desired output file name)
df_final.to_csv('cleaned_file.csv', index=False)

print("Duplicates removed, and data saved to 'cleaned_file.csv'")
