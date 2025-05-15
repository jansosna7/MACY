import numpy as np
import galois
from sympy import Matrix

# ==== USTAWIENIA ====
GF = galois.GF(2**2)     # Ciało skończone GF(2^4)
n = 3                   # Długość kodu RS (max dla GF(2^4): 2^4 - 1 = 15)
k = 2                   # Rozmiar wiadomości 

# ==== GENEROWANIE KODU RS ====
rs = galois.ReedSolomon(n, k, field=GF)
G = rs.G  # Macierz generatora (k x n)

print("Macierz generatora G (kod RS):")
print(G)

# ==== GENEROWANIE KLUCZA PUBLICZNEGO McEliece'a ====
# Losowa odwracalna (nieosobliwa) macierz S (k x k)
while True:
    S = GF.Random((k, k))
    if np.linalg.matrix_rank(S.view(np.ndarray)) == k:
        break

# Losowa permutacja P (n x n)
perm = np.random.permutation(n)
P = np.eye(n, dtype=int)[perm]  # permutacja jako macierz permutacyjna
P_GF = GF(P)

# Macierz publiczna: G' = S * G * P
G_pub = S @ G @ P_GF

print("\nPubliczna macierz G':")
print(G_pub)

# ==== ATAK SIDELNIKOVA-SHESTAKOVA ====
# Krok 1: Obliczenie kodu dualnego G_pub
def null_space(mat):
    M = Matrix(mat.astype(np.int8).tolist())
    ns = M.nullspace()
    if not ns:
        return np.array([])
    return np.array([v.tolist() for v in ns])

# Zamieniamy G_pub na postać binarną (każdy element GF(2^4) → wektor długości 4)
G_pub_bin = G_pub.vector(dtype=int)
G_pub_bin = G_pub_bin.reshape(k, n * GF.degree)

# Obliczanie kodu dualnego (macierzy H takiej, że H @ G_pub_bin.T = 0)
H = null_space(G_pub_bin)

print("\nPierwszy wektor kodu dualnego G_pub:")
if H.size > 0:
    print(H[0])
else:
    print("Kod dualny pusty — problem z obliczeniem nullspace.")

# ==== EDUKACYJNA SYMULACJA ODTWORZENIA PERMUTACJI ====
# Odtwarzamy permutację z G_pub → porównujemy kolumny z oryginalnym G
recovered = np.linalg.inv(S) @ G_pub @ P_GF.T

print("\nSprawdzenie, czy udało się odzyskać G (G_recovered ≈ G):")
print(np.all(recovered == G))  # True? oznacza poprawne odzyskanie
