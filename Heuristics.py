from BloomFilter import BloomFilter

# Function that returns the true positives guessed by the first heuristic (starting with
# already found true positives and adding elements until reaching the expected number of true
# positives in the filter).
# This heuristic accesses each of the slots of one element given by each hash function
# and sums all the quantities. The lesser this heuristic is, the more probable an element
# is a true positive.
# bf is the BloomFilter we want to extract the true positives from.
# p is the array of elements returned positive by the bf excluding the true positives already found.
# tp is the set of elements we know they are in the bf
# e is the array of elements associated to each position of the filter by one hash function
# n is the expected number of elements in the bf
def h_number_of_neighbours(bf, p, tp, e, n):
    # If we already have more than the expected number of elements, we don't do anything
    if len(tp) >= n:
        return set()
    # We calculate the value of the heuristic for each element in p
    heuristics = {}
    for element in list(p):
        value = 0
        for i in range(bf.nhash):
            value += len(e[bf.get_hash().getbit_idx(element, i)]) - 1
        heuristics[element] = value
    # We sort them and get the n-len(tp) elements with less heuristic that will be the predicted tp by the heuristic
    sorted_elements = sorted(heuristics, key=heuristics.get)
    return sorted_elements[:n-len(tp)]

# Function that returns the true positives guessed by the second heuristic (starting with
# already found true positives and adding elements until reaching the expected number of true
# positives in the filter).
# This heuristic counts the number of already found true positives in any of the positions of each element.
# The lesser this heuristic is, the more probable an element is a true positive.
# bf is the BloomFilter we want to extract the true positives from.
# p is the array of elements returned positive by the bf excluding the true positives already found
# tp is the set of elements we know they are in the bf
# e is the array of elements associated to each position of the filter by one hash function
# n is the expected number of elements in the bf
def h_number_of_tp(bf, p, tp, e, n):
    # If we already have more than the expected number of elements, we don't do anything
    if len(tp) >= n:
        return tp
    # We calculate the value of the heuristic for each element in tp, but in this heuristic
    # we should iterate through the tp and modifying the values of the heuristic of the rest
    # of the elements
    heuristics = {}
    for element in p:
        heuristics[element] = 0
    for true_positive in list(tp):
        for i in range(bf.nhash):
            for element in e[bf.get_hash().getbit_idx(true_positive, i)]:
                # We don't take the already found tp into account
                if element in tp:
                    continue
                heuristics[element] += 1
    # We sort them and get the n-len(tp) elements with less heuristic that will be the predicted tp by the heuristic
    sorted_elements = sorted(heuristics, key=heuristics.get)
    return sorted_elements[:n-len(tp)]