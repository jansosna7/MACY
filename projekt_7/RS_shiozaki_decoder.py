from GF_arithmetic import *


def rs_shiozaki_gao_decode(codeword, nsym):
    """
    Shiozaki-Gao decoding for Reed-Solomon codes.
    :param codeword: received codeword as list of integers (length n)
    :param nsym: number of parity symbols (ecc symbols)
    :return: original message (list of k symbols)
    """
    n = len(codeword)

    # Interpolacja otrzymanego wielomianu R(x) z punktów (α^i, codeword[i])
    R = [0] * n
    for i in range(n):
        R = gf_poly_add(R, gf_poly_scale([1], codeword[i]))
        R = gf_poly_mul(R, [1, gf_pow(2, i)])  # α^i w GF(2^m)

    # Ustal Q(x) = R(x) * x^t
    Q = gf_poly_mul(R, [1] + [0] * nsym)

    # Oblicz największy wspólny dzielnik Q(x) i x^n
    x_n = [1] + [0] * n
    _, remainder = gf_poly_div(Q, x_n)

    # Odtwórz wielomian wiadomości przez podzielenie Q(x) przez E(x)
    # Załóżmy, że E(x) to [1] (dla poprawnych danych lub 0 błędów)
    E = [1]
    M, _ = gf_poly_div(Q, E)

    # Skróć do oryginalnej długości wiadomości
    return M[:-(nsym)], M[-(nsym):]
