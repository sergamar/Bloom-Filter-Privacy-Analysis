# Black-box Privacy Analysis of Counting Bloom Filters
Repository that contains all the code related to the author's Master Thesis for the Master in Cybersecurity UC3M. The results are also part of the paper:

Sergio Galán, Pedro Reviriego, Stefan Walzer, Alfonso Sánchez-Macian, Shanshan Liu and Fabrizio Lombardi, "On the Privacy of Counting Bloom Filters for Black-box Attackers", IEEE Transactions on Dependable and Secure Computing (in press).

The code related to Counting Bloom Filters is based on the code from Alfonso Sánchez-Macián (https://github.com/amacian/invertCBF). It has been modified to implement a regular Bloom Filter and a Counting Bloom Filter without collisions.

Exp1 folder contains the code for an experiment that didn't make it to the paper since the results were not as interesting as initially expected. The algorithm tries to extract elements from a regular Bloom Filter in a white-box scenario by using some heuristics.

Exp2 folder contains the implementation of a CBF without collisions and the algorithm that tries to extract as many elements as it can from a CBF. In it, there is also the code for the empirical checking that black-box with pair extraction is equivalent (in the sense they extract the same elements given the same CBF) to a white-box algorithm restricted to counters with values less or equal to 2.

Exp3 folder contains the main experiment, in which we analyze how many elements we can extract from different CBFs (changing the number of False Positives, using filters with low, high, or optimal load) with our black-box algorithm with pair extraction. It uses the white-box algorithm limited to counters with value 2 or less instead of the black-box one for the sake of efficiency (since in the previous experiment we saw they extracted the same elements).
