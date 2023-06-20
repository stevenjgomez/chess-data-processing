import os
import sys


#########################################################
######################### INPUTS ########################

# Directory with temperature folders
data_dir = r'/nfs/chess/id4baux/2023-2/pokharel-3700-a/ACS-5-11-8-2/ACS-5-11-8-2/'

#########################################################
#########################################################


folders=[] # Empty list to store folder names
for item in os.listdir(data_dir):
  try:
    folders.append(int(item)) # If folder can be converted to int, add it
  except:
    pass # Otherwise don't add it
folders.sort() # Sort from low to high T
folders = [str(i) for i in folders] # Convert to strings


print('\n\nAvailable folders: ')
[print("["+str(i)+"] "+str(name)) for i,name in enumerate(folders)]

orm_selection_str = input("Enter an index to use for orientation matrix: ")
orm_folder = folders[int(orm_selection_str)]
working_dir_orm = data_dir + orm_folder + '/'

selection_str = input("Enter a comma-separated list of indices to select for indexing: ")
indices = [int(i) for i in selection_str.split(',')]
folders = [folders[i] for i in indices]


print('\n\nSelected folders: ' + str(folders))



val = input("Do you want to solve new orientation matrix using\n" + working_dir_orm + "? (y/n): ")
if val=='y':
   os.system("python auto_ormfinder.py "+working_dir_orm)

if True:
  for folder in folders:
    os.system("python Pil6M_HKLConv_3D.py "+ data_dir + folder + '/' + " " + working_dir_orm)
else:
  pass