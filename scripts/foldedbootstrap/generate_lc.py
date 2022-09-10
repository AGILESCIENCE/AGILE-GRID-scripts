import os,sys,math,random
import numpy as np

def pop_random(lst):
    idx = random.randrange(0, len(lst))
    return lst.pop(idx)

mode = sys.argv[1]
input_path = sys.argv[2]


# start from eight BIN files
if mode == "1":

    os.chdir(input_path)
    
    os.system("rm merged_bins.txt")
    os.system("rm ordered_merged_bins.txt")

    os.system("cat BIN* > merged_bins.txt")
    os.system("sort merged_bins.txt > ordered_merged_bins.txt")

    file_path = "ordered_merged_bins.txt"

# start from an already ordered and merged BIN file
elif mode == "2":
    file_path = input_path

starting_data = []

# apri il file
count = 0
for line in open(file_path):

    starting_data.append(line)
    count = count + 1

os.system("rm BIN1; rm BIN2; rm BIN3; rm BIN4; rm BIN5; rm BIN6; rm BIN7; rm BIN8")
os.system("rm random_BIN1.txt; rm random_BIN2.txt; rm random_BIN3.txt; rm random_BIN4.txt; rm random_BIN5.txt; rm random_BIN6.txt; rm random_BIN7.txt; rm random_BIN8.txt")

print("Total bins= "+str(count))

bin1 = open("random_BIN1.txt", "a")
bin2 = open("random_BIN2.txt", "a")
bin3 = open("random_BIN3.txt", "a")
bin4 = open("random_BIN4.txt", "a")
bin5 = open("random_BIN5.txt", "a")
bin6 = open("random_BIN6.txt", "a")
bin7 = open("random_BIN7.txt", "a")
bin8 = open("random_BIN8.txt", "a")

bin_length = math.floor(count/8)

print("File length = "+str(bin_length))

index_used = []

for index in range(0,bin_length):
    #print(index)    

    bin1.write(pop_random(starting_data))
    bin2.write(pop_random(starting_data))
    bin3.write(pop_random(starting_data))
    bin4.write(pop_random(starting_data))
    bin5.write(pop_random(starting_data))
    bin6.write(pop_random(starting_data))
    bin7.write(pop_random(starting_data))
    bin8.write(pop_random(starting_data))

print("remaining data= "+str(len(starting_data)))

bin1.close()
bin2.close()
bin3.close()
bin4.close()
bin5.close()
bin6.close()
bin7.close()
bin8.close()

os.system("sort random_BIN1.txt > BIN1")
os.system("sort random_BIN2.txt > BIN2")
os.system("sort random_BIN3.txt > BIN3")
os.system("sort random_BIN4.txt > BIN4")
os.system("sort random_BIN5.txt > BIN5")
os.system("sort random_BIN6.txt > BIN6")
os.system("sort random_BIN7.txt > BIN7")
os.system("sort random_BIN8.txt > BIN8")
os.system("rm random*.txt")

