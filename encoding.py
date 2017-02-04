from liblll import *
from math import sqrt
import sys

def swap_columns(B, i, j):
    if i != j:
        for k in range(len(B)):
            tmp = B[k][i]
            B[k][i] = B[k][j]
            B[k][j] = tmp

# Returns a matrix B' having the same columns of B, but in increasing order of l2_norm.
# So, if i <= j, then l2_norm(B'[i]) <= l2_norm(B'[j])
def sort_colunms_by_l2norm(B):
    m = len(B[0]) # number of columns
    for i in range(m-1):
        min_pos = i
        min_l2_norm = norml2(get_vector(B,i))
        for j in range(i+1, m):
            current_l2_norm = norml2(get_vector(B,j))
            if current_l2_norm < min_l2_norm:
                min_pos = j
                min_l2_norm = current_l2_norm
        swap_columns(B, i, min_pos)

def remove_three_last_entries(v):
    return [v[i] for i in range(len(v) - 3)]

#   Return the first value j such that the j-th column of
# B has a nonzero value at the position n-1 (so, the n-th position,
# since the indices start from 0).
#   If no such j exists, the value -1 is returned.
def find_first_column_with_nonzero_at_nth_position(B):
    n_minus_one = len(B) - 4
    for j in range(len(B[0])):
        if B[n_minus_one][j] != 0:
            return j
    return -1

def print_as_polynomial(vec):
    null_poly = True
    for i in range(len(vec)-1, -1, -1):
        if vec[i] > 0:
            sys.stdout.write(' +%dZ^%d' % (vec[i], i))
            null_poly = False
        elif vec[i] < 0:
            sys.stdout.write(' %dZ^%d' % (vec[i], i))
            null_poly = False
    if null_poly:
        print "0"
    else:
        print ""


def vector_norm(v):
    return sqrt(norml2(v))

def transpose(A):
    return [list(i) for i in zip(*A)]

def root_of_unity(n):
    return math.e**(math.pi*1j/n)

#   C: the scaling factor
#   n: number of coefficients
def scalled_real_coefficients(C, n):
    s = root_of_unity(n)
    return [ int(round((C*s**k).real)) for k in range(n)]

#   C: the scaling factor
#   n: number of coefficients
def scalled_complex_coefficients(C, n):
    s = root_of_unity(n)
    return [ int(round((C*s**k).imag)) for k in range(n)]

#   alpha: the complex number to be approximated
def construct_lattice_basis(alpha, C, T, n):
    B = [ [ 0 for i in range(n+1) ] for j in range(n+3) ]
    for i in range(n):
        B[i][i] = 1
    B[n][n] = T
    a_is = scalled_real_coefficients(C, n)
    b_is = scalled_complex_coefficients(C, n)
    a = int(round((C*alpha).real))
    b = int(round((C*alpha).imag))
    for k in range(n):
        B[n+1][k] = a_is[k]
    B[n+1][n] = -a
    for k in range(n):
        B[n+2][k] = b_is[k]
    B[n+2][n] = -b
    return B


def encode(alpha, C, T, n):
    B = construct_lattice_basis(alpha, C, T, n)
    R = lll_reduction(B)
    if not islll(R):
        print "ERROR! the following basis could not be reduced!"
        print_mat(B)
        sys.exit(0)

    sort_colunms_by_l2norm(R)
    l = find_first_column_with_nonzero_at_nth_position(R)
    if l == -1:
        print "ERROR! no vector in the reduced basis has n-th position != 0"
        sys.exit(0)
    
    v = remove_three_last_entries(get_vector(R, l))
    return v

def decode(vec, N):
    root_unity = math.e**(math.pi * 1j / N)
    pow_root_unity = 1
    val = 0j
    for i in range(len(vec)):
        val += vec[i] * pow_root_unity 
        pow_root_unity *= root_unity
    return val


C = 10**10
T = 10
n = 16
alpha = 0.655981733221013+0.923883055400882j  # example from the paper
print "Encoding ", alpha

v = encode(alpha, C, T, n)

print "Encoded to: "
print_as_polynomial(v)

print "Decoded to:"
decoded_value = decode(v, n)
print decoded_value

print "<Original value>  -  <Decoded value>"
print alpha - decoded_value

print "<Original value> /  <Decoded value>"
print alpha / decoded_value
