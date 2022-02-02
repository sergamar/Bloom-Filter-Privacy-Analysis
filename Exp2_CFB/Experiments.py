import random
import sys
import getopt
from CountingBloomFilter import CountingBloomFilter
from GenericHashFunctionsSHA512 import GenericHashFunctionsSHA512
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time

# Testing parameters
filter_size = 65536
# filter_size = 16384
k = 3
n = 14870
# n = 3718
max_val = 100000
falses = n * 4
trials = 1

# Function to generate the random set of elements.
# Current version uses strings
# num is the number of elements to generate
# if a CountingBloomFilter bf is passed, elements are stored there
# if a list ds is passed, elements are stored there
# exclude set are a elements that should not be selected
# maxVal is the maximum int value
def generate_random_elements(num, bf=None, ds=None, max_val=1000000000, exclude=None):
    # Elements are added to the set to check for repetitions
    if exclude is None:
        exclude = set()
    s = set()
    s.clear()

    # Keeps data of stored elements
    stored = 0
    # Generate elements until the value "stored" is reached
    while stored < num:
        # Generate integers between 1 and max_val
        entry = random.randint(1, max_val)
        # if entry was already selected or is in the exclude set,
        # go to next iteration
        if entry in s or entry in exclude:
            continue
        # When an BloomFilter is received
        if bf is not None:
            # Add the entry to the filter
            bf.add(entry)
        # When a list is received
        if ds is not None:
            # Add the element to the list
            ds.append(entry)
        # Another element has been stored
        stored = stored + 1
        # Add it to the set so they are not repeated
        s.add(entry)
    return

# Function to generate a random set of false positives
# num is the number of false positives to generate
# bf is the CountingBloomFilter in which to test the elements to find false positives
# maxVal is the maximum int value
# exclude set are a elements that should not be selected (basically true positives)
# returns a list with the generated false positives
def generate_random_fp(num, bf, max_val=1000000000, exclude=None):
    # Elements are added to the set to check for repetitions
    if exclude is None:
        exclude = set()
    s = set()
    s.clear()

    # Keeps data of stored elements
    stored = 0
    false_positives = []
    # Generate elements until the value "stored" is reached
    while stored < num:
        # Generate integers between 1 and max_val
        entry = random.randint(1, max_val)
        # if entry was already selected or is in the exclude set,
        # go to next iteration
        if entry in s or entry in exclude or not bf.check(entry):
            continue
        # When an CountingBloomFilter is received
        false_positives.append(entry)
        # Another element has been stored
        stored = stored + 1
        # Add it to the set so they are not repeated
        s.add(entry)
    return false_positives

# Function to find all elements from the universe that returns a positive from CBF
# bf is the Counting Bloom Filter
# max_val is the maximum integer value. Universe will include elements from 1 to max_val
def find_p(bf, max_val):
    # Create the list P of (true and false) positive elements
    p = list()
    # Check all elements of the universe, from 1 to max_val
    for i in range(max_val+1):
        # If one of the positions is 0, then it is a negative
        # Otherwise, add it to P
        if bf.check(i, 1):
            p.append(i)

    return p

# Function to find all elements from a set that returns a positive from CBF
# bf is the Counting Bloom Filter
# set is the set of elements to be tested against the filter
def find_p_set(bf, set):
    # Create the list P of (true and false) positive elements
    p = list()
    # Check all elements of the set
    for element in set:
        # If one of the positions is 0, then it is a negative
        # Otherwise, add it to P
        if bf.check(element, 1):
            p.append(element)

    return p

# Function that decides whether an posible true positive (by removing it and doing
# some checkings with the rest of the positives) is really a true positive or unknown
# Returns True when it is a tp, False if unknown
# element is the element to be tested
# bf is the Counting Bloom Filter
# positives_set is the set of all positives accepted by the filter
def test_element(element, bf, positives_set, removals):
    # We remove the element to be tested
    bf.remove(element)
    # And obtain the difference between the original positives and the new ones
    new_positives = set(find_p_set(bf, list(positives_set)))
    diff = positives_set - new_positives
    # If the only difference is the element itself, it is a true positive
    if len(diff) == 1:
        removals.add(element)
        return True
    # Otherwise, we see what happens if we add all the elements that disappeared
    # from the filter not taking the element itself into account
    elif len(diff) > 1:
        diff = diff - {element}
        for e in list(diff):
            bf.add(e)
        # If x isn't in the filter, then it is a true positive
        if not bf.check(element):
            for e in list(diff):
                bf.remove(e)
                removals.add(e)
            removals.add(element)
            return True
        # If it is in the filter, we can't say it is a TP for sure
        for e in list(diff):
            bf.remove(e)
        bf.add(element)
        return False
    # Otherwise, we can't decide whether it is a true positive or a false one
    else:
        bf.add(element)
        return False

x_axis = []
y_axis = []

dots = [x*(n//10) for x in range(0,51)]

for fals in dots:
    avg = 0
    print("FPs:", fals)
    for _ in range(trials):
        # Generate a standard bloom filter with the testing parameters
        bf = CountingBloomFilter(filter_size, k)

        # Fill the filter with random elements
        true_positives = []
        generate_random_elements(n, bf, true_positives, max_val)
        # print(true_positives)

        # Generate a certain number of false positives
        false_positives = generate_random_fp(fals, bf, max_val, true_positives)

        print("Start:", time.ctime(time.time()))
        # Run the algorithm
        all_positives = true_positives + false_positives
        all_positives_set = set(all_positives)
        found_tps = []
        new_tp_found = True
        rounds = 0
        while new_tp_found:
            new_tp_found = False
            removals = set()
            for pos in list(all_positives_set):
                # If we have already removed a tp or fp, we don't take it into account
                if pos in removals:
                    continue
                if test_element(pos, bf, all_positives_set, removals):
                    found_tps.append(pos)
                    new_tp_found = True
            all_positives_set = all_positives_set - removals
            rounds += 1
        print("Rounds:", rounds)
        avg += len(found_tps)/len(true_positives)
        print("Finish:", time.ctime(time.time()))

    x_axis.append(fals/n)
    y_axis.append(avg/trials*100)


print(x_axis)
print(y_axis)

plt.xlabel('Ratio False positives/True positives')
plt.ylabel('% of successfully obtained elements')
plt.title('m='+str(filter_size)+', '+'k='+str(k)+', '+'n='+str(n))
plt.yticks([0,20,40,60,80,100])
plt.ylim(bottom=0)
plt.xlim(left=0)
plt.xlim(right=4.25)
plt.grid()
plt.plot(x_axis, y_axis)
plt.savefig('initial_test_65536.png')





