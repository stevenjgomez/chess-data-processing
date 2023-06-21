import os
import shutil


def copy_file_to_cwd(file_name):
    # Get the path of the currently executing module file
    module_path = os.path.abspath(__file__)

    # Determine the installation directory by going up one level
    installation_dir = os.path.dirname(module_path)

    # Specify the source file path
    source_file = os.path.join(installation_dir, file_name)

    # Specify the destination file path in the current working directory
    destination_file = os.path.join(os.getcwd(), file_name)

    # Copy the file from the installation directory to the current working directory
    shutil.copyfile(source_file, destination_file)

    print(f"{file_name} copied successfully.")

def checkccode():
    # Get the path to the user's home directory
    home_dir = os.path.expanduser("~")

    # Create the full path for the "chess_hkl_scripts" folder within the home directory
    folder_path = os.path.join(home_dir, "chess_hkl_scripts")

    # Create the folder if it doesn't already exist
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        print(f"The 'chess_hkl_scripts' folder has been created in {folder_path}.")
    else:
        print(f"hkl scripts directory: {folder_path}.")

    os.chdir(folder_path)

    # Get the current working directory
    cwd = os.getcwd()

    copy_file_to_cwd("hkl.c")
    copy_file_to_cwd("../README.md")

    # Define the file names to check
    file_names = ["hkl.o", "libhkl.dll"]

    # Check if the files exist in the cwd
    for file_name in file_names:
        file_path = os.path.join(cwd, file_name)
        if os.path.exists(file_path):
            print(f"{file_name} found.")
        else:
            print(f"{file_name} has not been compiled.")
            print("Please reference README for instructions on how to compile hkl scripts from c code.")
            print(f"c code can be found in {cwd}")
            break
