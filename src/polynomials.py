from __future__ import annotations
import numpy as np
import time as time
import itertools

class polynomial:
    # Class of polynomials in orthogonal basis, as described on the paper.
    # It verifies that r = len(self._idx) and 2^r = len(self._coef)

    _idx: np.array # Sorted list of idx for which the polynomial depends
    _coef: np.array # List of coefficients. 

    def __init__(self, idx: np.array, coef: np.array = None) -> None:
        # Function to initiallize the class
        if coef is None:
            self._coef = np.zeros(2**idx.size, dtype=np.int64)
        elif 2**idx.size == coef.size: 
            self._coef = coef
        else:
            raise Exception("The coef array does not correspond to idx array.")
        self._idx = np.sort(idx)

    def add_bc(self, mono_idx: np.array, coef: int = 1) -> None:
        # Add the monomial in standard basis to the polynomial, expressed in the orthogonal basis. This 
        # algorithm substitutes Algorithm 2 of the paper.
        if coef == 0:
            return
        if not np.all(np.in1d(mono_idx, self._idx)): # Check that the polynomial is defined in the same space
            for idx in np.setdiff1d(mono_idx, self._idx):
                self.add_new_idx(idx)
        if not np.all(np.in1d(mono_idx, self._idx)):
            raise Exception("function error.")
        fixed_bits = int("".join(map(str, np.in1d(self._idx, mono_idx).astype(int))), 2) # Condition c_i = 1
        indices = np.arange(self._coef.size)
        mask = (indices & fixed_bits) == fixed_bits # Indices of the array corresponding to all c_i = 1
        self._coef[mask] += coef
        return self
    
    def add_ec(self, c: np.array, coef: int = 1) -> None:
        # Adds a polynomial expressed on the orthogonal basis. 
        self._coef[int(''.join(c.astype(str)), 2)]+=coef
        return self
    
    def add_new_idx(self, idx: int) -> None:
        # Represents the same polynomial with an extra index.
        p = polynomial(self._idx, self._coef)
        extension = extend(idx, p, p)
        self._idx = extension._idx
        self._coef = extension._coef
        return self
        
    def __iadd__(self, other) -> polynomial:
    # Performs addition
        if isinstance(other, polynomial):
            for idx in np.setdiff1d(other._idx, self._idx):
                self.add_new_idx(idx)
            for idx in np.setdiff1d(self._idx, other._idx):
                other.add_new_idx(idx)
            self._coef += other._coef
            return self
        if isinstance(other, int):
            self._coef += other # Elementwise sum
            return self
        raise Exception("Operation not supported")
    
    def __mul__(self, other) -> polynomial:
        # Performs multiplication
        if isinstance(other, polynomial):
            for idx in np.setdiff1d(other._idx, self._idx):
                self.add_new_idx(idx)
            for idx in np.setdiff1d(self._idx, other._idx):
                other.add_new_idx(idx)
            self._coef *= other._coef # Elementwise product
            return self
        if isinstance(other, int):
            self._coef *= other
            return self
        raise Exception("Operation not supported")
    
    def eval(self, idx: int, value: int) -> polynomial:
        # Evaluates the polynomial in the given idx to 0 or 1 by just selecting the correct coefficients 
        # from the coefficients list.
        if idx in self._idx and value in [0,1]:
            i = np.where(self._idx == idx)[0][0]
            l = 2**(self._idx.size-i-1)
            m = 2**i
            new_idx = np.delete(self._idx, i)
            new_coef = self._coef.reshape(m, 2 * l)
            if value == 0:
                coef = new_coef[:, :l].flatten()
            elif value == 1:
                coef = new_coef[:, l:].flatten()
            return polynomial(new_idx, coef)
        if not idx in self._idx:
            raise Exception("Invalid evaluation: idx = "+str(idx)+", poly defined over "+str(self._idx))
        if not value in [0,1]:
            raise Exception("Invalid evaluation: you tried to eval in "+str(value)+" but only evaluations in 0 and 1 are supported.")

 
    def print(self, terms: str = "") -> str:
        # Returns a string that represents the polynomial in the standard basis. 
        # Function described on the paper as Algorithm 7.
        word = "0"
        last_idx = self._idx[-1]
        if self._idx.size == 1:
            a = self._coef[0]
            b = self._coef[1]-self._coef[0]
            if a == 0 and b != 0:
                word = str(b)+"x_"+str(last_idx)+terms
            elif a != 0 and b == 0:
                word = str(a)+terms
            elif a != 0 and b != 0:
                word = str(a)+terms+" + "+str(b)+"x_"+str(last_idx)+terms
            return word
        p0 = self.eval(last_idx, 0)
        terms0 = p0.print(terms)
        p0 *= -1
        p1 = self.eval(last_idx, 1)
        p1 += p0
        terms1 = p1.print("x_"+str(last_idx)+terms)
        if terms0 == "0" and terms1 != "0":
            word = terms1
        elif terms0 != "0" and terms1 == "0":
            word = terms0
        elif terms0 != "0" and terms1 != "0":
            word = terms0+" + "+terms1
        return word
    
    def __repr__(self) -> str:
        # Procedure to print the polynomial
        return self.print()

class automorphism:
    # Class automorphism, defined as a list of images (d) for each phi(e_c).
    _domain_idx: np.array
    _image_idx: np.array
    _image: np.array # List of images

    def __init__(self, domain_idx: np.array, image_idx: np.array) -> None:
        # Initialize the class
        if domain_idx.size > image_idx.size:
            extend_img_idx = np.sort(np.setdiff1d(np.arange(1, 33), image_idx))
            self._image_idx = np.sort(np.concatenate((image_idx, extend_img_idx[:domain_idx.size - image_idx.size])))
            self._domain_idx = np.sort(domain_idx)
        else:
            extend_domain_idx = np.sort(np.setdiff1d(np.arange(1, 33), domain_idx))
            self._domain_idx = np.sort(np.concatenate((domain_idx, extend_domain_idx[:image_idx.size - domain_idx.size])))
            self._image_idx = np.sort(image_idx)
        n = self._image_idx.size
        m = 2**(n)
        self._image = np.zeros((m,n), dtype=domain_idx.dtype)
    
    def set_img(self, c: np.array, d: np.array) -> None: 
        # Set psi(e_c) = e_d
        c = np.array(c, dtype=int)
        d = np.array(d, dtype=int)
        self._image[int(''.join(c.astype(str)), 2)]=d
        return self

    def apply(self, poly: polynomial) -> polynomial:
        # Computes psi(poly), which is a polynomial in the orthogonal basis. Part of the signing algorithm.
        for idx in np.setdiff1d(poly._idx, self._domain_idx): # psi is extended in domain (phi(P) variables) to all Q's domain
            self.extend_psi(idx)
        for idx in np.setdiff1d(self._domain_idx, poly._idx): # Q is extended to all P variables
            poly.add_new_idx(idx)
        result = polynomial(self._image_idx)
        n = 2**self._domain_idx.size
        for i in range(n):
            d = self._image[i]
            result.add_ec(d, poly._coef[i])
        return result
    
    # c and d can be integers and then simplify the function (which gives wrong outputs)
    def include_partial_psi(self, psi_c: automorphism, c: np.array, d: np.array) -> None:
        # Given that psi(e_c)=e_d and we found a solution psi_c to a smaller problem, join this solution to our actual problem.
        # Use the relation: psi(e_{c||c'})=e_d·psi_c(e_{c'}), as described on the paper, Algorithm 5.
        domain_extension_idx = np.setdiff1d(self._domain_idx, psi_c._domain_idx)
        image_extension_idx = np.setdiff1d(self._image_idx, psi_c._image_idx)
        if domain_extension_idx.size + psi_c._domain_idx.size != self._domain_idx.size:
            raise Exception("Wrong domain extension.")
        if image_extension_idx.size + psi_c._image_idx.size != self._image_idx.size:
            raise Exception("Wrong image extension.")
        domain_idx_order = np.argsort(np.concatenate((psi_c._domain_idx, domain_extension_idx)))
        image_idx_order = np.argsort(np.concatenate((psi_c._image_idx, image_extension_idx)))
        if c.size != domain_extension_idx.size:
            raise Exception("The given evaluation doesn't correspond to the reduction on idx.")
        n = psi_c._domain_idx.size
        i=0
        for binary_tuple in itertools.product([0, 1], repeat=n):
            c_vec = np.concatenate((np.array(binary_tuple), c))[domain_idx_order]
            d_vec = np.concatenate((psi_c._image[i], d))[image_idx_order]
            self.set_img(c_vec, d_vec)
            i += 1
        return self
    
    def extend_psi(self, new_idx: int) -> automorphism:
        # Extend the definition of psi to a new index.
        if new_idx in self._domain_idx: 
            print("Error call")
            return self
        psi = automorphism(np.concatenate((self._domain_idx, np.array([new_idx]))), self._image_idx)
        psi.include_partial_psi(self, np.array([0]), np.array([0]))
        psi.include_partial_psi(self, np.array([1]), np.array([1]))
        self._domain_idx = psi._domain_idx
        self._image_idx = psi._image_idx
        self._image = psi._image
        return self

def extend(new_idx: int, p0: polynomial, p1: polynomial) -> polynomial:
    # Define a new polynomial corresponding to p = x_l·p1+(1-x_l)·p0, where l = new_idx.
    if new_idx in p0._idx:
        raise Exception("The new idx is already in the polynomial.")
    if np.all(p0._idx != p1._idx):
        raise Exception("The polynomials used to extend should be defined in the same space.")
    idx = np.sort(np.append(p0._idx, new_idx))
    i = np.where(idx == new_idx)[0][0]
    l = 2**(idx.size-i-1)
    m = 2**i
    new_p0_coef = p0._coef.reshape(m, l)
    new_p1_coef = p1._coef.reshape(m, l)
    coef = np.empty((2 * m, l), dtype=new_p0_coef.dtype)
    coef[0::2] = new_p0_coef
    coef[1::2] = new_p1_coef
    coef = coef.flatten()
    return polynomial(idx, coef)
