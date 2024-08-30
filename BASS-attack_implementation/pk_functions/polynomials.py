from __future__ import annotations

class monomial:
    # Class to be used only to define a poly
    _degree: int
    _idx: list[int]
    _coef: int

    def _in_order(self) -> None:
        # Sort the idx of the monomial
        self._idx.sort()

    def __init__(self, idx: list[int] = [], coef: int = 1) -> None:
        # Initialize the class
        self._degree = len(idx)
        self._idx = idx
        self._in_order()
        self._coef = coef
    
    def sum_coef(self, c: int) -> None:
        # Add a constant to the coefficient
        self._coef+=c
    
    def coef(self) -> int:
        # Obtain the coefficient
        out = self._coef
        return out
    
    def idx(self) -> list[int]:
        # Return a copy of the indices that appear on the monomial
        out = []
        [out.append(x) for x in self._idx]
        return out
    

    def copy(self) -> monomial:
        # Returns a copy of the monomial
        return monomial(self.idx(), self.coef())

    def print(self) -> str:
        # Returns a string that represents the monomial
        if len(self._idx) == 0: return str(self._coef) #If it is an independent term
        word = str(self._coef)+"·" if self._coef != 1 else ""
        first = True
        for term in self._idx:
            if not first: word+="·"
            first=False
            word+="x_"+str(term)
        return word
    
    def __mul__(self, other: monomial) -> monomial:
        # Monomial multiplication
        self._coef*=other.coef()
        for idx in other.idx():
            if idx not in self.idx():
                self._idx.append(idx)
        self._in_order()
        return self

    def __repr__(self) -> str:
        # Procedure to print the monomial
        return self.print()

    def __eq__(self, mono: monomial) -> bool:
        # Compare two monomials and return true if they are equal. Two monomials are equal if they
        # depend on the same variables, regardless of its coefficient.
        if self._degree != mono._degree: return False
        if self._idx != mono._idx: return False
        return True
    
    def __lt__(self, mono: monomial) -> bool:
        # Compare two monomials and returnt true if the second one is less than the first one. 
        # Procedure used on sorting. 
        if self._degree < mono._degree: return True
        if self._degree > mono._degree: return False
        for i in range(self._degree):
            if self._idx[i] < mono._idx[i]: return True
            if self._idx[i] > mono._idx[i]: return False
        return False

class poly:
    # Poly in the standard basis, implemented as a list of monomials
    _monomials: list[monomial]
    _indexes: list[int]

    def _in_order(self) -> None:
        # Orders the vector of monomials and indexes
        self._monomials.sort()
        self._indexes.sort()

    def find_idx(self) -> None:
        # Computes the vector of indexes
        self._indexes = []
        for mono in self._monomials:
            for i in mono.idx():
                if i not in self._indexes:
                    self._indexes.append(i)
            
    def idx(self) -> list[int]:
        # Return a list of the indices for which this polynomial depends on
        return self._indexes.copy()
    
    def _sum_repeated(self) -> None:
        # Sum coefficient of those repeated monomials
        self._in_order()
        res = [monomial([], 0)]
        for i in range(len(self._monomials)):
            mono=self._monomials[i]
            if mono.coef() != 0: 
                if mono == res[-1]: 
                    res[-1].sum_coef(mono.coef())
                    if res[-1].coef() == 0 and len(res)>1: 
                        res.pop() #maybe erasing a variable
                else: 
                    res.append(mono.copy())
        self._monomials=res
        self.find_idx()
        return

    def __init__(self, p: poly = None) -> None:
        # Init class
        self._indexes = []
        self._monomials = [monomial([], 0)]
        self._ec_base = []
        self._ec_idx = []

    def copy(self) -> poly:
        # Returns a copy of this polynomial
        out = poly()
        out._indexes = self._indexes.copy()
        out._monomials = []
        [out._monomials.append(x.copy()) for x in self._monomials]
        return out

    def __iadd__(self, other) -> poly:
        # Performs a sum
        if isinstance(other, poly): 
            # Add poly to itself
            self._monomials += other._monomials
            self._indexes = list(set(self.idx())&set(other.idx()))
        elif isinstance(other, monomial): 
            # Add monomial to itself
            self._monomials.append(other)
            self._indexes = list(set(self.idx())&set(other.idx()))
        elif isinstance(other, int): 
            # Add constant term monomial to itself
            if len(self._monomials[0].idx()) != 0: 
                raise Exception("Poly no longer has an independent term")
            self._monomials[0].sum_coef(other)
        else:
            raise Exception("Operation + not supported between polynomial and "+str(type(other)))
        self._sum_repeated()
        return self

    def __mul__(self, other) -> poly:
        # Performs a multiplication. First convert the given type to a polynomial and then perform 
        # standard polynomial multiplication (O(m^2))
        if isinstance(other, poly): 
            p=other          
        elif isinstance(other, monomial): 
            p = poly()
            p += other
        elif isinstance(other, int): 
            p = poly()
            p += other
        else:
            raise Exception("Operation not supported.")
        mul = poly()
        mul._ec_base = self._ec_base
        mul._ec_idx = self._ec_idx
        for mono1 in p._monomials:
            for mono2 in self._monomials:
                m1=mono1.copy()
                m2=mono2.copy()
                m3=m1*m2
                mul += m3     
        mul._sum_repeated()         
        return mul

    def print(self, small_print: bool = True) -> str:
        # Returns a string that represents the polynomial
        if len(self._monomials) > 20 and small_print: # If it is too large, don't print it all
            return "(poly with "+str(len(self._monomials))+" terms)"
        first = True
        word=""
        for mono in self._monomials:
            if not first: word+=" + "
            first=False
            word+=mono.print()
        return word
    
    def __repr__(self) -> str:
        # Procedure to print the polynomial
        return self.print()


class automorphism:
    # Class automorphism, defined as a list of images (polynomials) for each phi(x_i).
    # This class is only used to generate the public key, the attack uses the implementation 
    # described on the paper and can be found on /src/polynomials.py
    _defined: list[bool] # True at position i if the automorphism is defined at x_i
    _image: list[poly] # Image of each x_i

    def __init__(self) -> None:
        # Initiaize the class
        self._defined = [False]*32
        self._image = []
        for i in range(32):
            self._image.append(poly())
        self._ec_img = []
    
    def copy(self) -> automorphism:
        # Return a copy of that automorphism
        out = automorphism()
        for i in range(32):
            if self._defined[i]:
                out.set_img(i, self._image[i].copy())
        return out
    
    def set_img(self, idx: int, img: poly = poly()) -> None:
        # Set the image of x_{idx} to the given polynomial
        self._defined[idx] = True
        self._image[idx] = img
    
    def apply(self, p: poly) -> poly:
        # Compute phi(p)
        image = poly()
        for mono in p._monomials:
            if mono.coef()==0:
                if len(mono.idx())==0:
                    pass
                else:
                    raise Exception("Poly with a coefficient 0 in a monomial.")
            mul = poly()
            mul += monomial([], mono.coef())
            for idx in mono.idx():
                if self._defined[idx]:
                    mul *= self._image[idx]
                else:
                    raise Exception("This automorphism can't be applied to this poly since it is not defined on some used variable.")
            image += mul
        return image
    
    def print(self) -> str:
        # Returns a string that represents the automorphism. Also indicates the number of simple images.
        word = ""
        count = 0
        for i in range(32):
            if self._defined[i]:
                if len(self._image[i]._monomials) == 2:
                    count += 1
                word += "x_"+str(i)+" -> "+self._image[i].print()+"\n"
        word += "Number of simple images: "+str(count)+"\n"
        return word

    def __repr__(self) -> str:
        # Procedure to print the automorphism
        return self.print()
    
    def __imul__(self, a: automorphism) -> automorphism:
        # Automorphism composition. First applies a, then self.
        self._defined = a._defined
        image = self.copy()
        for idx in range(32):
            if a._defined[idx]:
                self._image[idx] = image.apply(a._image[idx])
        return self        
