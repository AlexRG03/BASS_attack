# BASS_attack
Implementation of the BASS attack described in (link paper)

The 'pk_functions' folder contains the data structures and functions described in the BASS paper, such as key generation, signing and verification functions.

The 'src' folder contains the data structures and functions described in our paper, along with functions to read the generated public key and parse it. 

The 'experiments' folder contains all data extracted from the experiments we mention in section 5, together with the code used to generate such data.

The jupyter nootebook file 'create_pk' uses functions from 'pk_functions' to generate the public key, a random message and its corresponding Q polynomial. These polynomials are stored in a file called 'pk.txt' inside the 'src' folder. 'empirical_results' was used to estimate the value of $\delta$ we define in section 5, and also contains code to run fast experiments. 'attack' is the file containing the entire process for generating a forgery described in section 6 of our paper. 
