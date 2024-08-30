from polynomials import *
from typing import Tuple
import re

def read_python_pk(file_path: str) -> Tuple[list[polynomial], polynomial]:
    # Returns a list of 6 polynomials, corresponding to p_i - a_i, phi(p_i) - a_i and the polynomial Q.
    with open(file_path, 'r') as file:
        content = file.read()
    lines = content.splitlines()
    start = time.time()
    p1, a1 = find_polynomial_python(lines[0])
    p2, a2 = find_polynomial_python(lines[1], 0)
    p3, a3 = find_polynomial_python(lines[2], 0)
    end = time.time()
    print("All p_i calculated in: ", end-start, " seconds.")
    start = time.time()
    phi_p1, _ = find_polynomial_python(lines[3], a1)
    end = time.time()
    print("phi_p1 calculated in: ", end-start, " seconds. Size: 2^", phi_p1._idx.size)
    start = time.time()
    phi_p2, _ = find_polynomial_python(lines[4], a2)
    end = time.time()
    print("phi_p2 calculated in: ", end-start, " seconds. Size: 2^", phi_p2._idx.size)
    start = time.time()
    phi_p3, _ = find_polynomial_python(lines[5], a3)
    end = time.time()
    print("phi_p3 calculated in: ", end-start, " seconds. Size: 2^", phi_p3._idx.size)
    start = time.time()
    q, _ = find_polynomial_python(lines[6], 0)
    end = time.time()
    print("Q calculated in: ", end-start, " seconds. Size: 2^", q._idx.size)
    p1 += a1
    p2 += a2
    p3 += a3
    return [p1,p2,p3,phi_p1,phi_p2,phi_p3], q

def find_polynomial_python(poly_str: str, add: int = 0) -> Tuple[polynomial, int]:
    # Given the polynomial in string format, convert it to the class 'polynomial' in orthogonal basis.
    # Return the number of negative terms as well
    parts = poly_str.split(" + ")
    # Find all idx on the second monomial as the first one is the constant term
    if len(parts) < 2:
        raise Exception("Constant polynomial cannot be read: "+poly_str)
    mono = parts[1]+"·"
    matches = re.findall(r'x_(\d+)(?=\·)', mono)
    idx = list(set(int(match) for match in matches))

    # Declare the output
    poly = polynomial(np.array(idx))
    poly += add
    poly += int(parts[0])
    a = 0

    # Explore all monomials 
    for m in parts[1:]:
        coef_str = m.split("x")[0]
        if coef_str == "":
            coef = 1
        elif coef_str[0] == '-':
            a += 1
            coef = -int(coef_str[1:-1])
        else:
            coef = int(coef_str[:-1])
        m += "·"
        m = m[len(coef_str):]
        # Convert string to array of idx 
        matches = re.findall(r'x_(\d+)(?=\·)', m)
        idx = [int(match) for match in matches]
        poly.add_bc(np.array(idx), coef)
    return poly, a

