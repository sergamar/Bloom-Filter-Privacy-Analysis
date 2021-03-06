import random
import sys
import getopt
from CountingBloomFilterNoCol import CountingBloomFilterNoCol
from GenericHashFunctionsSHA512 import GenericHashFunctionsSHA512
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
parser.add_argument("-t", dest="t", type=int, help="Number of iterations (default 10)", default=10)
parser.add_argument("-k", dest="k", type=int, help="Number of hashes (default 3)", default=3)
args = parser.parse_args()
filter_size = args.m
n = args.n
trials = args.t
k = args.k
pairs = 1

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
# m is the size of the filter
# elements is the T array
# positives is the list of elements to be removed
# count_cbf is the list of counters from the CBF
# count is the list of counters from T
# hashf is the hash function used in the CBF
# k is the number of positions
# is_positive indicates if it is a real positive (true) or a false positive (false)
# nocol is a boolean indicating whether no collision is activated or not
def clear_positions(m, elements, positives, count_cbf, count, hashf, k, is_positive, nocol):
    # Additional elements to be removed
    additional = list()
    # and iterate over them
    num = len(positives)

    # Traverse the positives list
    for i in range(num):
        # get next element to be removed
        next_positive = positives[i]
        # for the k hash functions
        hashes = []
        for j in range(k):
            # Get the position mapped for the element and the jth hash function
            jpos = hashf.getbit_idx(next_positive, j)
            if nocol:
                while jpos in hashes:
                    jpos = (jpos + 1) % m
                hashes.append(jpos)
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
        clear_positions(m, elements, additional, count_cbf, count, hashf, k, False, nocol)

    return

# Get the set of elements from the Universe (1 to maxVal) that may be included into the filter.
# It won't extract elements in counters with 3 or more elements (for better comparison with black box analysis)
# m is the number of positions (counters) of the CBF
# k is the number of hash functions
# cbf is the Counting Bloom Filter
# p is the P array with all the elements from the universe that returned positive from CBF
# nocol is a boolean indicating whether no collision is activated or not
def peeling(m, k, cbf, p, nocol):
    # Retrieve the hash function used
    hashf = cbf.get_hash()
    # T array  with m positions to store the elements that are mapped to them
    elements = [None] * m
    # Count of elements mapped to each position
    count = [0] * m

    # For all the positions in p
    for i in range(len(p)):
        hashes = []
        # for the k hash functions
        for j in range(k):
            # Get the position mapped for the element p[i] and the jth hash function
            pos = hashf.getbit_idx(p[i], j)
            # If no collision is activated, we recalculate the value of the hash
            if nocol:
                while pos in hashes:
                    pos = (pos + 1) % m
                hashes.append(pos)
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
            # and where that number is less than 3
            if count[i] != counters[i] or counters[i] >= 3:
                continue
            # we found at least one position
            found = True
            # add the elements of the ith position to the positives
            positives.update(elements[i])
            # all these elements must be removed as well, but not from elements
            removers = elements[i].copy()
            # call the function that clears the removers and related false positives
            # pass True as last parameter as they are real positives
            clear_positions(m, elements, removers, counters, count, hashf, k, True, nocol)
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
    # If pos2 is negative after we remove pos1 it is not useful for us
    if not bf.check(pos2):
        bf.add(pos1)
        return False
    # We also check the reciprocal
    bf.add(pos1)
    bf.remove(pos2)
    if not bf.check(pos1):
        bf.add(pos2)
        return False
    bf.remove(pos1)
    # Check that both a now negatives
    if  bf.check(pos1) or bf.check(pos2):
        bf.add(pos1)
        bf.add(pos2)
        return False

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
        temp_added = []
        for e in list(diff):
            bf.add(e)
            temp_added.append(e)
            if not bf.check(e):
                for r in temp_added:
                    bf.remove(r)
                bf.add(pos1)
                bf.add(pos2)
                return False
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

for fals in dots:
    avg_blackbox_ind = 0
    worst_blackbox_ind = 100
    avg_blackbox_pairs = 0
    worst_blackbox_pairs = 100
    avg_whitebox = 0
    worst_whitebox = 100
    for _ in range(trials):
        # First we proceed with the blackbox analysis

        # Generate a standard CBF with the testing parameters
        # We create a no colision CBF since we are performing pair extraction
        bf = CountingBloomFilterNoCol(filter_size, k)

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
                    all_positives_set_temp = all_positives_set - removals
            all_positives_set = all_positives_set - removals
        # Check that we haven't labeled a FP as a TP
        for z in found_tps:
            if z not in true_positives:
                print("ERROR: Algorithm labeled FP as TP in simple filtering")
                exit(0)
        # Check if we want to extract pairs
        if pairs:
            # When we can't get new TP one by one, we proceed with pairs
            new_tp_found = True
            found_tps_with_pairs = found_tps.copy()
            while new_tp_found:
                new_tp_found = False
                removals = set()
                all_positives_set_temp = all_positives_set.copy()
                for pos1 in list(all_positives_set):
                    for pos2 in list(all_positives_set):
                        if pos1 == pos2:
                            continue
                        # If we have already removed a tp or fp, we don't take it into account
                        if pos1 in removals or pos2 in removals:
                            continue
                        if test_pairs(pos1, pos2, bf, all_positives_set_temp, removals):
                            found_tps_with_pairs.append(pos1)
                            found_tps_with_pairs.append(pos2)
                            new_tp_found = True
                            all_positives_set_temp = all_positives_set - removals
                all_positives_set = all_positives_set - removals
                # If we didn't find a new TP with pairs, we finish the extraction
                if new_tp_found == True:
                    continue
                iterations = 0
                new_tp_found = True
                # If we found a new TP, we try to extract new TPs with the usual method
                while new_tp_found:
                    new_tp_found = False
                    removals = set()
                    all_positives_set_temp = all_positives_set.copy()
                    for pos in list(all_positives_set):
                        # If we have already removed a tp or fp, we don't take it into account
                        if pos in removals:
                            continue
                        if test_element(pos, bf, all_positives_set_temp, removals):
                            found_tps_with_pairs.append(pos)
                            new_tp_found = True
                            all_positives_set_temp = all_positives_set - removals
                    all_positives_set = all_positives_set - removals
                    iterations += 1
                # Check that we haven't labeled a FP as a TP
                for z in found_tps_with_pairs:
                    if z not in true_positives:
                        print("ERROR: Algorithm labeled FP as TP in pair filtering")
                        exit(0)
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

        if pairs:
            prct_obtained_pairs = (len(found_tps_with_pairs)/len(true_positives)) * 100
            avg_blackbox_pairs += prct_obtained_pairs/trials
            if prct_obtained_pairs < worst_blackbox_pairs:
                worst_blackbox_pairs = prct_obtained_pairs

        # Then, we carry out the whitebox analysis
        # Generate a new CBF with the same data
        bf = CountingBloomFilterNoCol(filter_size, k)
        for posit in true_positives:
            bf.add(posit)
        found_tps = peeling(filter_size, k, bf, all_positives, pairs)
        prct_obtained = (len(found_tps)/len(true_positives)) * 100
        avg_whitebox += prct_obtained/trials
        if prct_obtained < worst_whitebox:
            worst_whitebox = prct_obtained
        if set(found_tps_with_pairs) != found_tps:
            # print(set(found_tps_with_pairs))
            # print(found_tps)
            # print(len(found_tps_with_pairs))
            # print(len(found_tps))
            # print(found_tps - set(found_tps_with_pairs))
            # for ee in list(found_tps - set(found_tps_with_pairs)):
            #     print("-------------", ee, "-------------")
            #     hashes = []
            #     for i in range(bf.nhash):
            #         idx = bf.hash.getbit_idx(ee, i)
            #         while idx in hashes:
            #             idx = (idx + 1) % m
            #         print(idx)
            print("ERROR: Different number of elements extracted")
            exit(0)

    print("Tested with " + str(fals) + " false positives")

print("Blackbox with pairs and whitebox limited to two elements have extracted the same elements.")




