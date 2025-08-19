import os

def read_all_files_in_directory(directory):
    file_contents = {}
    
    # Iterate through all files in the given directory
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        
        # Check if it's a file
        if os.path.isfile(file_path):
            with open(file_path, 'r', encoding='utf-8') as f:
                file_contents[filename] = f.read()
    
    return file_contents