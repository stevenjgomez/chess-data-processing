import os
from Pil6M_HKLConv_3D import Pil6M_HKLConv_3D
from auto_ormfinder import auto_ormfinder


def orient_index(H, K, L, data_dir, orm_temp, index_temps=None, stacks_list=None):
    if stacks_list is None:
        stacks_list = [1, 2, 3]
    if index_temps is None:
        # Identify all temperature folders
        folders = []  # Empty list to store folder names
        for item in os.listdir(data_dir):
            try:
                folders.append(int(item))  # If folder can be converted to int, add it
            except ValueError:
                pass  # Otherwise don't add it
        folders.sort()  # Sort from low to high T
        folders = [str(folder) for folder in folders]  # Convert to strings
    else:
        folders = index_temps

    print('\n\nSelected folders: ' + str(folders))

    ormdir = data_dir + orm_temp + '/'

    # Check if orientation matrix exists already
    file_path = os.path.join(ormdir, "ormatrix_auto.nxs")
    if os.path.isfile(file_path):
        print(f"\nOrientation matrix found at {file_path}")
        val = input("Resolve orientation matrix? (y/n): ")
    else:
        val = 'y'

    # Solve for orientation matrix
    if val == 'y':
        print(f"\nNew orientation matrix will be saved to {file_path}")
        auto_ormfinder(ormdir)

    # Check if orientation matrix exists again before performing indexing
    if os.path.isfile(file_path):
        # Index specified folders
        for folder in folders:
            workingdir = os.path.join(data_dir, folder)
            Pil6M_HKLConv_3D(H, K, L, workingdir, ormdir, stacks_list)
    else:
        raise FileNotFoundError(f"Orientation matrix not found in {file_path}.")
