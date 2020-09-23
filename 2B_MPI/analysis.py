# Using readlines() 
file1 = open('mseq-big13-example.dat', 'r') 
Lines = file1.readlines() 

# writing to file 
file2 = open('mseq-small13-example.dat', 'w') 



count = 0
# Strips the newline character 
for line in Lines: 
    print("Line{}: {}".format(count, len(line)))
    if count < 3:
        file2.write(line)
    else:
        until = (len(line)-1)//10
        print(until)
        file2.write(line[0:until]) 
        file2.write("\n")
    count += 1
file2.close()