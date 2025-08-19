import csv
import sys

def csv_to_txt(csv_file, txt_file):
    with open(csv_file, mode='r', newline='') as infile, open(txt_file, mode='w') as outfile:
        reader = csv.reader(infile)
        for row in reader:
            # Join the row elements with a tab or space, change as needed
            outfile.write('\t'.join(row) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python csv_to_txt.py input.csv output.txt")
    else:
        csv_to_txt(sys.argv[1], sys.argv[2])
