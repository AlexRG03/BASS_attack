from polynomials import *
import numpy as np
import hashlib
from itertools import islice

def get_p() -> poly:
    # Generate a t-sparse polynomial, with t = 3
    out = poly()
    for _ in range(3): # Define each monomial as described on the BASS paper
        d=np.random.randint(0, 3)+1
        remaining_idx = list(range(0,31))
        i1 = np.random.randint(0, len(remaining_idx))
        idx1 = remaining_idx[i1]
        remaining_idx.remove(idx1)
        i2 = np.random.randint(0, len(remaining_idx))
        idx2 = remaining_idx[i2]
        remaining_idx.remove(idx2)
        i3 = np.random.randint(0, len(remaining_idx))
        idx3 = remaining_idx[i3]
        terms = [idx1, idx2,idx3]
        coef = 2*np.random.randint(0, 2)-1
        m = monomial(terms[0:d], coef)
        out += m
    return out

def get_h(i: int = 0, big: bool = True) -> poly:
    # Gives a polynomial h such that it only depends on idx bigger (or smaller if false) than i.
    if big:
        remaining_idx = list(range(i+1,31))
    else:
        remaining_idx = list(range(0,i))
    if len(remaining_idx) < 3: 
        return poly()
    i1 = np.random.randint(0, len(remaining_idx))
    idx1 = remaining_idx[i1]
    remaining_idx.remove(idx1)
    i2 = np.random.randint(0, len(remaining_idx))
    idx2 = remaining_idx[i2]
    remaining_idx.remove(idx2)
    i3 = np.random.randint(0, len(remaining_idx))
    idx3 = remaining_idx[i3]
    p = poly()
    if np.random.binomial(1,0.5) == 0: # Choose between M or 1-M
        p += 1
        mono = monomial([idx1, idx2], -1)
        p += mono
    else: 
        mono = monomial([idx1, idx2], 1)
        p += mono
    q = poly()
    if np.random.binomial(1,0.5) == 0: # Choose between multiplying by x_i or (1-x_i)
        q += 1
        mono = monomial([idx3], -1)
        q += mono
    else: 
        mono = monomial([idx3], 1)
        q += mono
    p*=q
    return p

def generate_phi() -> automorphism: 
    n = 31
    alpha = automorphism()
    beta = automorphism()
    pi = automorphism()
    for i in range(n): # Define alpha
        p = poly()
        p+=monomial([i])
        if np.random.binomial(1,0.5) == 0:
            h = get_h(i)
            p+=h.copy()
            h*=-2
            h*=monomial([i])
            p+=h
        alpha.set_img(i,p.copy())
    for i in range(n-1, -1, -1): # Define beta
        p = poly()
        p+=monomial([i])
        beta._defined[i] = True
        if np.random.binomial(1,0.5) == 0:
            h = get_h(i, False)
            p+=h.copy()
            h*=-2
            h*=monomial([i])
            p+=h
        beta.set_img(i,p.copy())
    remaining_idx = list(range(31))
    for i in range(n): # Define the permutation pi
        idx = np.random.randint(0, len(remaining_idx))
        p=poly()
        p+=monomial([remaining_idx[idx]])
        pi.set_img(i,p.copy())
        remaining_idx.remove(remaining_idx[idx])
    # Compose the automorphisms: phi = pi(beta(alpha))
    beta*=alpha
    pi*=beta
    return pi

def get_Q(s: str) -> poly:
    # Get the polynomial Q from a message
    sha256_hash = hashlib.sha3_256(s.encode()).digest()
    bitstring = ''.join(format(byte, '08b') for byte in sha256_hash)
    it = iter(bitstring)
    # We divide the string into 32 blocks of 8 bits
    blocks_8_bits = iter(lambda: tuple(islice(it, 8)), ())
    idx = 0
    q = poly()
    xs = list(range(32))
    for block in blocks_8_bits: # For each block, we create a monomial and sum it to q
        coeff_bits, exponents = block[0:5], block[5:8]
        coeff = sum(int(bit) for bit in coeff_bits)%3
        coeff = -1 if coeff == 2 else coeff
        mono_idx = []
        for c in exponents:
            if c == '1':
                mono_idx.append(xs[idx])
            idx += 1
            if idx >= 32:
                idx = 0 
        q += monomial(mono_idx, coeff)
    return q