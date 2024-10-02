# BASS_attack
Implementation of the BASS attack described in (link paper)

The 'pk_functions' folder contains the data structures and functions described in the BASS paper, such as key generation, signing and verification functions.

The 'src' folder contains the data structures and functions described in our paper, along with functions to read the generated public key and parse it. 

The jupyter nootebook file 'create_pk' uses functions from 'pk_functions' to generate the public key, a random message and its corresponding Q polynomial. These polynomials are stored in a file called 'pk.txt' inside the 'src' folder. 'empirical_results' was used to generate the data from which we claim empirical observations in section 5 of our paper, such as the acceptance rate of a trivial signature and the value of delta. 'attack' is the file containing the entire process for generating a forgery described in section 6 of our paper. 
