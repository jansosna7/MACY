from GF_arithmetic import *


def rs_shiozaki_gao_decode(codeword, nsym):
    n = len(codeword)

    R = [0] * n
    for i in range(n):
        R = gf_poly_add(R, gf_poly_scale([1], codeword[i]))
        R = gf_poly_mul(R, [1, gf_pow(2, i)])  # α^i w GF(2^m)

    Q = gf_poly_mul(R, [1] + [0] * nsym)

    x_n = [1] + [0] * n
    _, remainder = gf_poly_div(Q, x_n)

    E = [1]
    M, _ = gf_poly_div(Q, E)

    return M[:-(nsym)], M[-(nsym):]
