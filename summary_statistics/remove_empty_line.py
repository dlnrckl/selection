import os
import sys

def remove_empty_line(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Find the line with 'positions' and remove the empty line after it
    for i, line in enumerate(lines):
        if 'positions' in line:
            if lines[i + 1].strip() == "":
                del lines[i + 1]  # Remove the empty line
            break
    
    # Write the updated lines back to the file
    with open(file_path, 'w') as file:
        file.writelines(lines)

def process_all_ms_files(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".ms"):
            file_path = os.path.join(directory, filename)
            print(f"Processing: {file_path}")
            remove_empty_line(file_path)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <directory_path>")
    else:
        directory = sys.argv[1]
        process_all_ms_files(directory)
