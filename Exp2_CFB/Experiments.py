import random
import sys
import getopt
from CountingBloomFilter import CountingBloomFilter
from GenericHashFunctionsSHA512 import GenericHashFunctionsSHA512
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
import argparse

# Testing parameters
# filter_size = 1024
# k = 4
# n = 175
# trials = 100
max_val = 1000000000

parser = argparse.ArgumentParser()
parser.add_argument("-m", dest="m", type=int, help="Filter size (default 1024)", default=1024)
parser.add_argument("-n", dest="n", type=int, help="Number of true positives (default 256)", default=256)
parser.add_argument("-t", dest="t", type=int, help="Number of iterations (default 100)", default=100)
parser.add_argument("-k", dest="k", type=int, help="Number of hashes (default 3)", default=3)
args = parser.parse_args()
filter_size = args.m
n = args.n
trials = args.t
k = args.k

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

# Function that clears the element from its positions in T and also clears all the related false positives
# elements is the T array
# positives is the list of elements to be removed
# count_cbf is the list of counters from the CBF
# count is the list of counters from T
# hashf is the hash function used in the CBF
# k is the number of positions
# is_positive indicates if it is a real positive (true) or a false positive (false)
def clear_positions(elements, positives, count_cbf, count, hashf, k, is_positive):
    # Additional elements to be removed
    additional = list()
    # and iterate over them
    num = len(positives)

    # Traverse the positives list
    for i in range(num):
        # get next element to be removed
        next_positive = positives[i]
        # for the k hash functions
        for j in range(k):
            # Get the position mapped for the element and the jth hash function
            jpos = hashf.getbit_idx(next_positive, j)
            # Element might have been removed in a different level of recursion
            if elements[jpos].count(next_positive) == 0:
                break
            # Remove the element from the position
            elements[jpos].remove(next_positive)
            # Reduce the T counter for that position
            count[jpos] -= 1
            # Reduce the CBF counter only when it is a real positive
            if is_positive:
                count_cbf[jpos] -= 1
            # If no more elements are mapped to this position in the CBF
            # we can remove all the pending elements from T and they are false positives
            if count_cbf[jpos] == 0 and count[jpos] != 0:
                # Add those elements to additional list
                additional.extend(elements[jpos])

    # Recursive call to remove the false positive elements
    if len(additional) > 0:
        # Pass False as last parameter as they are false positives
        clear_positions(elements, additional, count_cbf, count, hashf, k, False)

    return

# Get the set of elements from the Universe (1 to maxVal) that may be included into the filter.
# m is the number of positions (counters) of the CBF
# k is the number of hash functions
# cbf is the Counting Bloom Filter
# p is the P array with all the elements from the universe that returned positive from CBF
def peeling(m, k, cbf, p):
    # Retrieve the hash function used
    hashf = cbf.get_hash()
    # T array  with m positions to store the elements that are mapped to them
    elements = [None] * m
    # Count of elements mapped to each position
    count = [0] * m

    # For all the positions in p
    for i in range(len(p)):
        # for the k hash functions
        for j in range(k):
            # Get the position mapped for the element p[i] and the jth hash function
            pos = hashf.getbit_idx(p[i], j)
            # Retrieve the position pos of the T array
            list_pos = elements[pos]
            # If no elements are assigned to that position, create a list and assign it
            if list_pos is None:
                list_pos = list()
                elements[pos] = list_pos
            # Include the element into the list of elements mapped to the position
            list_pos.append(p[i])
            # Increase the count of elements mapped to the position
            count[pos] += 1

    # Set that will store the positives that were extracted from the filter
    positives = set()
    # Values for the CBF counters
    counters = cbf.get_counters()

    while True:
        # If we found an element that could be extracted in this iteration
        found = False
        # Go through the m positions
        # TODO: exclude positions that were 0 in the previous iterations to speed up the process
        for i in range(m):
            # when the ith position does not have values, go to next position
            if count[i] == 0:
                continue
            # the position has values => it is not empty
            # we only want those positions where we have the same number of elements in T and CBF
            if count[i] != counters[i]:
                continue
            # we found at least one position
            found = True
            # add the elements of the ith position to the positives
            positives.update(elements[i])
            # all these elements must be removed as well, but not from elements
            removers = elements[i].copy()
            # call the function that clears the removers and related false positives
            # pass True as last parameter as they are real positives
            clear_positions(elements, removers, counters, count, hashf, k, True)
        # if no new positives were found in the iteration, we should finish the algorithm
        if not found:
            break

    # return elements that were retrieved from the CBF
    return positives

# Function that decides whether a posible true positive (by removing it and doing
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
        positives_set = positives_set - {element}
        return True
    # Otherwise, we see what happens if we add all the elements that disappeared
    # from the filter not taking the element itself into account
    elif len(diff) > 1:
        diff = diff - {element}
        temp_added = []
        for e in list(diff):
            bf.add(e)
            temp_added.append(e)
            # If it's not in the filter, then we have a hash collision
            if not bf.check(e):
                for r in temp_added:
                    bf.remove(r)
                bf.add(element)
                return False
        # If x isn't in the filter, then it is a true positive
        if not bf.check(element):
            for e in list(diff):
                bf.remove(e)
                removals.add(e)
                positives_set = positives_set - {e}
            removals.add(element)
            positives_set = positives_set - {element}
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

# Function that decides whether a pair of posible true positives (by removing it and doing
# some checkings with the rest of the positives) is really a pair of true positives or unknown
# Returns True when it is a pair of tps, False if unknown
# pos1 is the first element of the pair
# pos2 is the second element of the pair
# bf is the Counting Bloom Filter
# positives_set is the set of all positives accepted by the filter
def test_pairs(pos1, pos2, bf, positives_set, removals):
    # We remove the pair of elements to be tested
    bf.remove(pos1)
    bf.remove(pos2)
    # And obtain the difference between the original positives and the new ones
    new_positives = set(find_p_set(bf, list(positives_set)))
    diff = positives_set - new_positives
    # If the only difference is the pair of elements, they are true positives
    if len(diff) == 2:
        removals.add(pos1)
        removals.add(pos2)
        return True
    # Otherwise, we see what happens if we add all the elements that disappeared
    # from the filter not taking the pair into account
    elif len(diff) > 2:
        diff = diff - {pos1}
        diff = diff - {pos2}
        for e in list(diff):
            bf.add(e)
        # If the pair isn't in the filter, then it is a true positive
        if not bf.check(pos1) and not bf.check(pos2):
            for e in list(diff):
                bf.remove(e)
                removals.add(e)
            removals.add(pos1)
            removals.add(pos2)
            return True
        # If one of them is in the filter, we can't say it is a TP for sure
        for e in list(diff):
            bf.remove(e)
        bf.add(pos1)
        bf.add(pos2)
        return False
    # Otherwise, we can't decide whether they are tp or fp
    else:
        bf.add(pos1)
        bf.add(pos2)
        return False

x_axis = []
y1_axis = []
y2_axis = []
y3_axis = []
y4_axis = []
y5_axis = []
y6_axis = []

dots = [(x*n)//10 for x in range(0,51)]

f = open(str(filter_size) + '_' + str(k) + '_' + str(n) + '.results', 'w')
f.write("Start: " + time.ctime(time.time()) + "\n")
for fals in dots:
    avg_blackbox_ind = 0
    worst_blackbox_ind = 100
    avg_blackbox_pairs = 0
    worst_blackbox_pairs = 100
    avg_whitebox = 0
    worst_whitebox = 100
    for _ in range(trials):
        # First we proceed with the blackbox analysis

        # Generate a standard bloom filter with the testing parameters
        bf = CountingBloomFilter(filter_size, k)

        # Fill the filter with random elements
        true_positives = []
        generate_random_elements(n, bf, true_positives, max_val)
        # print(true_positives)

        # Generate a certain number of false positives
        false_positives = generate_random_fp(fals, bf, max_val, true_positives)

        # Run the algorithm
        all_positives = true_positives + false_positives
        all_positives_set = set(all_positives)
        found_tps = []
        new_tp_found = True
        # We try to extract elements with the algorithm one by one
        while new_tp_found:
            new_tp_found = False
            removals = set()
            all_positives_set_temp = all_positives_set.copy()
            for pos in list(all_positives_set):
                # If we have already removed a tp or fp, we don't take it into account
                if pos in removals:
                    continue
                if test_element(pos, bf, all_positives_set_temp, removals):
                    found_tps.append(pos)
                    new_tp_found = True
            all_positives_set = all_positives_set - removals
        # Check that we haven't labeled a FP as a TP
        for z in found_tps:
            if z not in true_positives:
                print("ERROR: Algorithm labeled FP as TP")
                exit(0)
        # When we can't get new TP one by one, we proceed with pairs
        new_tp_found = True
        found_tps_with_pairs = found_tps.copy()
        while new_tp_found:
            new_tp_found = False
            removals = set()
            for pos1 in list(all_positives_set):
                for pos2 in list(all_positives_set):
                    # If we have already removed a tp or fp, we don't take it into account
                    if pos1 in removals or pos2 in removals:
                        continue
                    if test_pairs(pos1, pos2, bf, all_positives_set, removals):
                        found_tps_with_pairs.append(pos1)
                        found_tps_with_pairs.append(pos2)
                        new_tp_found = True
            all_positives_set = all_positives_set - removals
            # If we didn't find a new TP with pairs, we finish the extraction
            if new_tp_found == False:
                break
            iterations = 0
            # If we found a new TP, we try to extract new TPs with the usual method
            while new_tp_found:
                new_tp_found = False
                removals = set()
                for pos in list(all_positives_set):
                    # If we have already removed a tp or fp, we don't take it into account
                    if pos in removals:
                        continue
                    if test_element(pos, bf, all_positives_set, removals):
                        found_tps_with_pairs.append(pos)
                        new_tp_found = True
                all_positives_set = all_positives_set - removals
                iterations += 1
            # If we didn't find a new one, we finish
            if iterations == 1:
                break
            # Otherwise, we run pairs again
            new_tp_found = True

        # Record the results
        prct_obtained_ind = (len(found_tps)/len(true_positives)) * 100
        avg_blackbox_ind += prct_obtained_ind/trials
        if prct_obtained_ind < worst_blackbox_ind:
            worst_blackbox_ind = prct_obtained_ind

        prct_obtained_pairs = (len(found_tps_with_pairs)/len(true_positives)) * 100
        avg_blackbox_pairs += prct_obtained_pairs/trials
        if prct_obtained_pairs < worst_blackbox_pairs:
            worst_blackbox_pairs = prct_obtained_pairs

        # Then, we carry out the whitebox analysis
        # Generate a new CFB with the same data
        bf = CountingBloomFilter(filter_size, k)
        for posit in true_positives:
            bf.add(posit)
        found_tps = peeling(filter_size, k, bf, all_positives)
        prct_obtained = (len(found_tps)/len(true_positives)) * 100
        avg_whitebox += prct_obtained/trials
        if prct_obtained < worst_whitebox:
            worst_whitebox = prct_obtained

    print("FPs:", fals, "Avg Whitebox:", avg_whitebox, "Worst Whitebox:", worst_whitebox, "Avg Blackbox Ind:", avg_blackbox_ind, "Worst Blackbox Ind:", worst_blackbox_ind, "Avg Blackbox Pairs:", avg_blackbox_pairs, "Worst Blackbox Pairs:", worst_blackbox_pairs)

    x_axis.append(fals/n)
    y1_axis.append(avg_whitebox)
    y2_axis.append(avg_blackbox_ind)
    y3_axis.append(avg_blackbox_pairs)
    y4_axis.append(worst_whitebox)
    y5_axis.append(worst_blackbox_ind)
    y6_axis.append(worst_blackbox_pairs)


f.write("End: " + time.ctime(time.time()) + "\n")
f.write("Proportion FP/TP\n")
f.write(str(x_axis) + "\n")
f.write("Avg Whitebox\n")
f.write(str(y1_axis) + "\n")
f.write("Avg Blackbox Ind\n")
f.write(str(y2_axis) + "\n")
f.write("Avg Blackbox Pairs\n")
f.write(str(y3_axis) + "\n")
f.write("Worst Whitebox\n")
f.write(str(y4_axis) + "\n")
f.write("Worst Blackbox Ind\n")
f.write(str(y5_axis) + "\n")
f.write("Worst Blackbox Pairs\n")
f.write(str(y6_axis) + "\n")
f.close()

plt.xlabel('Ratio False positives/True positives')
plt.ylabel('% of successfully obtained elements')
plt.title('m='+str(filter_size)+', '+'k='+str(k)+', '+'n='+str(n))
plt.yticks([0,20,40,60,80,100])
plt.ylim(bottom=0)
plt.xlim(left=0)
plt.xlim(right=5.25)
plt.grid()
lab1 = "Avg Whitebox"
lab2 = "Avg Blackbox Ind"
lab3 = "Avg Blackbox Pairs"
lab4 = "Worst Whitebox"
lab5 = "Worst Blackbox Ind"
lab6 = "Worst Blackbox Pairs"
plt.plot(x_axis, y1_axis, label=lab1)
plt.plot(x_axis, y2_axis, label=lab2)
plt.plot(x_axis, y3_axis, label=lab3)
plt.plot(x_axis, y4_axis, label=lab4)
plt.plot(x_axis, y5_axis, label=lab5)
plt.plot(x_axis, y6_axis, label=lab6)
plt.legend()
plt.savefig(str(filter_size) + '_' + str(k) + '_' + str(n) + '.png')





