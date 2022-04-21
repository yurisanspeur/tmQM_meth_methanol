import shutil

# read file
file1 = open('./halides.txt', 'r')
Lines = file1.readlines()

# count = 0
xyz_file = []
for line in Lines:
    # count += 1
    line_data = line.split()
    if line_data[0] == 'monodentates':
        xyz_file.append(line_data[2])

# copy candidate xyz files
src_path = "./individual_xyz_data/"
dst_path = "./filtered_xyz_data/"
for file in xyz_file:
    shutil.copy(src_path+file, dst_path)