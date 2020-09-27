# Using readlines() 
file1 = open('hard.dat', 'r') 
Lines = file1.readlines() 
  
count = 0
# Strips the newline character 
for line in Lines: 
    print("Line{}: {}".format(count, len(line)))
    count += 1