from __future__ import annotations
from polynomials import *
from typing import Tuple

def precomputation(pk: list[polynomial]) -> Tuple[polynomial, polynomial]:
    # Given the modified public key, corresponding to a list of 6 polynomials: P_i - a_i and phi(P_i) - a_i, 
    # output the polynomials P and phi(P) as defined on the paper (section 4.2).
    start = time.time()
    p = polynomial(pk[0]._idx.copy(), pk[0]._coef.copy())
    p += pk[1]*4
    p += pk[2]*16
    end = time.time()
    print("P calculated in: ", end-start, " seconds, ", (end-start)/60, " minutes.")
    start = time.time()
    phi_p = polynomial(pk[3]._idx.copy(), pk[3]._coef.copy())
    phi_p += pk[4]*4
    phi_p += pk[5]*16
    end = time.time()
    print("phi_P calculated in: ", end-start, " seconds, ", (end-start)/60, " minutes.")
    start = time.time()
    extra_var = np.setdiff1d(np.arange(1,32), p._idx)[:phi_p._idx.size-p._idx.size]
    for var in extra_var:
        p.add_new_idx(var)
    end = time.time()
    print("P expanded in: ", end-start, " seconds, ", (end-start)/60, " minutes.")
    return p, phi_p

def class_computation(p: polynomial, n: int = 64) -> list[list[np.array]]:
    # Compute the sets S_i for all i, given the polynomial p, by inspection.
    result = [[] for _ in range(n)] # List that contains all S_i
    bit_width = (p._coef.size-1).bit_length()
    # Iterate over the array and classify indices, O(2^n)
    for idx, value in enumerate(p._coef): 
        # Convert index to binary array
        c = np.array(list(np.binary_repr(idx, width=bit_width)), dtype=int)
        # Add c to S_i
        result[value].append(c) 
    return result

def fifo(p: polynomial, phi_p: polynomial) -> automorphism:
    # Fifo algorithm, described on the paper as Algorithm 1 
    psi = automorphism(p._idx, phi_p._idx)
    S_i = class_computation(p)
    T_i = class_computation(phi_p)
    for i in range(64):
        size = len(S_i[i])
        for j in range(size):
            psi.set_img(S_i[i][j], T_i[i][j])
    return psi

class prediction:
    # Class of predictions
    x: int
    phi_x: int

    def __init__(self, x: int, phi_x: int) -> None:
        # Initialize the class
        self.x = x
        self.phi_x = phi_x
    
    def __repr__(self) -> str:
        # Procedure to print the prediction
        return "x_"+str(self.x)+" -> x_"+str(self.phi_x)
    
def find_predictions(p: polynomial, phi_p: polynomial, partial_pred: list[prediction] = []) -> list[prediction]:
    # Find all feasible solutions, given a list of already feasible solutions. Algorithm described on the 
    # paper as Algorithm 3 (with some modification)
    n = p._idx.size
    if phi_p._idx.size != n: 
        raise Exception("Wrong call.")
    hist_p = []
    hist_phi_p = []
    predictions = []
    for k in range(n):
        count = np.bincount(p.eval(p._idx[k], 0)._coef, minlength=64)
        hist_p.append(count) #to compute S_i,k,0 for all i
        count = np.bincount(phi_p.eval(phi_p._idx[k], 0)._coef, minlength=64)
        hist_phi_p.append(count) #to compute T_i,k,0 for all i
    for k in range(n):
        for l in range(n):
            pred = prediction(p._idx[k], phi_p._idx[l])
            if pred in partial_pred:
                if np.all(hist_p[k] == hist_phi_p[l]): # Check S_i,k,0 = T_i,l,0 for all i
                    predictions.append(pred)
    return predictions

def find_all_predictions(p: polynomial, phi_p:polynomial, all_pred: list[prediction] = []) -> list[prediction]:
    # Find all feasible and compatible solutions, given a list of already compatible solutions. Algorithm described on the 
    # paper as Algorithm 4 (with some modification)
    if len(all_pred) == 0:
        partial_pred = find_predictions(p, phi_p)
        if len(partial_pred) == 0:
            return []
        all_pred.append(partial_pred[0])
    else:
        partial_pred = all_pred.copy()
    while len(partial_pred) != 0: # Equivalent to a while True
        k = len(all_pred)
        for c in itertools.product([0, 1], repeat=k):
            p_c = p
            phi_p_c = phi_p
            for h in range(k):
                p_c=p_c.eval(all_pred[h].x, c[h])
                phi_p_c=phi_p_c.eval(all_pred[h].phi_x, c[h])
            partial_pred = find_predictions(p_c, phi_p_c, partial_pred)
            if len(partial_pred) == 0:
                return all_pred
        all_pred.append(partial_pred[0])
    return all_pred

def find_psi(p: polynomial, phi_p:polynomial, pred: list[prediction] = []) -> automorphism:
    # Algorithm that finds an automorphism psi as described on the paper as Algorithm 5.
    if p._idx.size < 4:
        return fifo(p, phi_p)
    pred = find_all_predictions(p, phi_p, pred)
    if len(pred) == 0:
        return fifo(p, phi_p)
    p0 = p.eval(pred[0].x, 0)
    phi_p0 = phi_p.eval(pred[0].phi_x, 0)
    psi0 = find_psi(p0, phi_p0, pred[1:])
    p1 = p.eval(pred[0].x, 1)
    phi_p1 = phi_p.eval(pred[0].phi_x, 1)
    if np.all(p0._coef == p1._coef) and np.all(phi_p0._coef == phi_p1._coef):
        psi1 = psi0
    else:
        psi1 = find_psi(p1, phi_p1, pred[1:])
    psi = automorphism(p._idx, phi_p._idx)
    psi.include_partial_psi(psi0, np.array([0]), np.array([0]))
    psi.include_partial_psi(psi1, np.array([1]), np.array([1]))
    return psi
