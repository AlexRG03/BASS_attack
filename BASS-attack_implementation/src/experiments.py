from polynomials import *

def validate(pk, q, signature):
    # Verification algorithm
    montecarlo_it = 3000
    n = 32
    # Define the positions on the vector of coefficients u such that the monomial
    # depends on x_i
    idx1 = np.array([1,5,6,7,11,12,13,15])
    idx2 = np.array([2,5,8,9,11,12,14,15])
    idx3 = np.array([3,6,8,10,11,13,14,15])
    idx4 = np.array([4,7,9,10,12,13,14,15])
    r_count = 0
    s_count = 0
    # Sample 3000 tuples
    random_binary_tuples = np.random.randint(2, size=(montecarlo_it, n))
    u = np.random.choice([-2, -1, 0, 1, 2], size=16)
    for c in random_binary_tuples:
        # Compute R, S evaluated on each tuple
        p0_c = pk[0]._coef[int(''.join(c[pk[0]._idx].astype(str)), 2)]
        p1_c = pk[1]._coef[int(''.join(c[pk[1]._idx].astype(str)), 2)]
        p2_c = pk[2]._coef[int(''.join(c[pk[2]._idx].astype(str)), 2)]
        phi_p0_c = pk[3]._coef[int(''.join(c[pk[3]._idx].astype(str)), 2)]
        phi_p1_c = pk[4]._coef[int(''.join(c[pk[4]._idx].astype(str)), 2)]
        phi_p2_c = pk[5]._coef[int(''.join(c[pk[5]._idx].astype(str)), 2)]
        q_c = q._coef[int(''.join(c[q._idx].astype(str)), 2)]
        signature_c = signature._coef[int(''.join(c[signature._idx].astype(str)), 2)]
        r = u.copy()
        r[idx1] *= p0_c
        r[idx2] *= p1_c
        r[idx3] *= p2_c
        r[idx4] *= q_c
        if np.sum(r) > 0: r_count += 1
        s = u.copy()
        s[idx1] *= phi_p0_c
        s[idx2] *= phi_p1_c
        s[idx3] *= phi_p2_c
        s[idx4] *= signature_c
        if np.sum(s) > 0: s_count += 1
    # Compute the difference and check if it is smaller than a 3%
    diff = abs(r_count-s_count)
    return diff*100 < 3*montecarlo_it

def get_coef_list(p: polynomial, idx_pq: np.array) -> list[np.array]:
    # Return a list of polynomials that corresponds to p evaluated on every boolean tuple of coefficients
    # given in idx_pq
    k = idx_pq.size
    new_idx = p._idx
    coef_list_p = [p._coef]
    for h in range(k): # eval in idx_pq[h] and save both results, for each array of coef
        idx = idx_pq[h]
        i = np.where(new_idx == idx)[0][0]
        l = 2**(new_idx.size-i-1)
        m = 2**i
        new_idx = np.delete(new_idx, i)
        new_coef_list = [[] for i in range(2**(h+1))]
        i=0
        for coef in coef_list_p:
            new_coef = coef.reshape(m, 2 * l)
            new_coef_list[i] = new_coef[:, :l].flatten()
            i += 1
            new_coef_list[i] = new_coef[:, l:].flatten()
            i += 1
        coef_list_p = new_coef_list.copy()
    return coef_list_p

def get_hist(p: polynomial, q: polynomial) -> np.array:
    idx_pq = np.intersect1d(p._idx, q._idx)
    k = idx_pq.size
    n = 32
    scale = 2**(n-p._idx.size-q._idx.size+idx_pq.size)
    hist_p = np.zeros((64, 64))
    coef_list_p = get_coef_list(p, idx_pq)
    coef_list_q = get_coef_list(q, idx_pq)
    for h in range(2**k): # For each boolean tuple we evaluated, compute the distribution of values
        count_p = np.bincount(coef_list_p[h], minlength=64)
        count_q = np.bincount(coef_list_q[h], minlength=64)
        # The histogram corresponds to the count of every c such that restricted to idx_pq has given values. 
        hist_p += np.outer(count_p, count_q)
    hist_p*=scale
    return hist_p.astype(int)

