import random
import sys
import getopt
from BloomFilter import BloomFilter
from GenericHashFunctionsSHA512 import GenericHashFunctionsSHA512


# Function to generate the random set of elements.
# Current version uses strings
# num is the number of elements to generate
# if a BloomFilter bf is passed, elements are stored there
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

# Function to find all elements from the universe that returns a positive from CBF
# bf is the Bloom Filter
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

# Function that clears the element from its positions
# elements is the T array
# positives is the list of elements to be removed
# count is the list of counters from T
# hashf is the hash function used in the BF
# k is the number of positions
def clear_positions(elements, positives, count_bf, count, hashf, k):
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
    return


# Function that tries to extract the elements of a given bloom filter
# m is the number of positions (counters) of the BF
# k is the number of hash functions
# bf is the Bloom Filter
# p is the P array with all the elements from the universe that returned positive from BF
def probabilistic_peeling(m, k, bf, p):
    # Retrieve the hash function used
    hashf = bf.get_hash()
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
    safe_positives = set()
    # Set that will store the positives that are not 100% certain
    unsafe_positives = set()

    # Values for the BF counters
    counters = bf.get_counters()

    # We first extract the direct true positives (lists with only one element)
    # while True:
    # If we found an element that could be extracted in this iteration
    found = False
    # Go through the m positions
    for i in range(m):
        # when the ith position does not have values, go to next position
        if count[i] == 0:
            continue
        # the position has values => it is not empty
        # we only want those positions where we have the same number of elements in T and BF
        if count[i] != counters[i]:
            continue
        # we found at least one position
        found = True
        # add the elements of the ith position to the positives
        safe_positives.update(elements[i])

        # TODO: Try to get more elements using heuristics

    # return elements that were retrieved from the BF
    return safe_positives


# Testing parameters
filter_size = 1024
hashes = 3
elements = 10
max_val = 100000

# Generate a standard bloom filter with the testing parameters
bf = BloomFilter(filter_size, hashes)

# Fill the filter with random elements
true_positives = []
generate_random_elements(elements, bf, true_positives, max_val)
print(true_positives)

# Get all the elements from the universe that are identified by the filter
positives = find_p(bf, max_val)
print(positives)

# Run the algorithm
found = probabilistic_peeling(filter_size, hashes, bf, positives)
print(found)

print("Found percentage: ", len(found)/len(true_positives)*100, "%")






