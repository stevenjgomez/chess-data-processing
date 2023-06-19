import os

for sample in []
  #################################################################
  # INPUTS
  folder_name = 'YO-45-2-7'
  sample_name = 'YO-45-2-7'
  specfile='/nfs/chess/id4b/2023-2/pokharel-3700-a/'+folder_name
  calibfile="/nfs/chess/id4baux/2023-2/pokharel-3700-a/calibration/LaB6_26keV.poni"
  maskfile="/nfs/chess/id4baux/2023-2/pokharel-3700-a/calibration/mask.edf"
  raw_dir = '/nfs/chess/id4b/2023-2/pokharel-3700-a/raw6M/'+folder_name+'/'+sample_name+'/'
  out_dir = '/nfs/chess/id4baux/2023-2/pokharel-3700-a/'+folder_name+'/'+sample_name+'/'
  #################################################################




  folders=[] # Empty list to store folder names
  for item in os.listdir(raw_dir):
    try:
      folders.append(int(item)) # If folder can be converted to int, add it
    except:
      pass # Otherwise don't add it
  folders.sort() # Sort from low to high T
  folders = [str(i) for i in folders] # Convert to strings


  print('\nRaw directory: ' + raw_dir)
  print('\nAvailable folders: ')
  [print("["+str(i)+"] "+str(name)) for i,name in enumerate(folders)]

  selection_str = input("Enter a comma-separated list of indices to select for stacking: ")
  indices = [int(i) for i in selection_str.split(',')]
  folders = [folders[i] for i in indices]

  print('\n\nSelected folders: ' + str(folders))




  val = input('Continue? (y/n): ')
  if val=='y':
    for folder in folders:
      scan_folders=[] # Empty list to store folder names
      for scan_folder in os.listdir(raw_dir+folder+'/'):
        scan_folders.append(scan_folder)
      scan_folders.sort()


      print('\nScan directory: ' + raw_dir + folder)
      print('\nAvailable folders: ')
      [print("["+str(i)+"] "+str(name)) for i,name in enumerate(scan_folders)]

      selection_str = input("Enter a comma-separated list of indices to select for stacking: ")
      indices = [int(i) for i in selection_str.split(',')]
      scan_folders = [scan_folders[i] for i in indices]

      print('\n\nSelected folders: ' + str(scan_folders))

      for scan_folder in scan_folders:
        scan_num=scan_folder[-3:]
        temperature=folder
        raw_path = raw_dir+folder+"/"+scan_folder+"/"
        out_path = out_dir+folder+"/"
        os.system("python SGA_stack_em_all_2022_batch.py"+" "+specfile+" "+calibfile+" "+maskfile+" "+sample_name+" "+scan_num+" "+temperature+" "+raw_path+" "+out_path)
  else:
    pass
