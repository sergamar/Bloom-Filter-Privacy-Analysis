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
    fps = []
    new_tps = []
    # We initialize two counters for each position, one for keeping track of true positives and
    # another to keep track of the false positives. Also, we update it with the tp already found
    # by Bianchi's method
    tp_counters = [0] * len(e)
    fp_counters = [0] * len(e)
    total_counters = [0] * len(e)
    for i, elements in enumerate(e):
        if not elements:
            continue
        total_counters[i] = len(elements)
    for found_tp in tp:
        for i, elements in enumerate(e):
            if not elements:
                continue
            if found_tp in elements:
                tp_counters[i] += 1
                elements.remove(found_tp)
                total_counters[i] -= 1

    # We take the elements that most true positives neighbours in all their positons have as false positives
    # Then, we try to extract new tps by Bianchi's method. We iterate until we find the required number of tps
    while len(new_tps) < n - len(tp) and not all(c == 0 for c in total_counters):
        heuristics = {}
        for i, elements in enumerate(e):
            if not elements:
                continue
            for element in elements:
                if element not in heuristics:
                    heuristics[element] = tp_counters[i]
                else:
                    heuristics[element] += tp_counters[i]
        maximum = max(heuristics.values())
        to_be_fp = [fp for fp,v in heuristics.items() if v == maximum]
        for fp in to_be_fp:
            for i, elements in enumerate(e):
                if not elements:
                    continue
                if fp in elements:
                    elements.remove(fp)
                    fp_counters[i] += 1
                    total_counters[i] -= 1
        fp_counters += to_be_fp
        to_be_tps = []
        for i, elements in enumerate(e):
            if not elements:
                continue
            if len(elements) == 1:
                to_be_tps.append(elements[0])
                elements.remove(elements[0])
                total_counters[i] = 0
        for to_be_tp in to_be_tps:
            for i, elements in enumerate(e):
                if not elements:
                    continue
                if to_be_tp in elements:
                    elements.remove(to_be_tp)
                    tp_counters[i] = 1
                    total_counters[i] -= 1
        new_tps += to_be_tps

    return new_tps