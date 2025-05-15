import numpy as np

def gram_schmidt(A):
    n = len(A)
    A_star = np.zeros_like(A, dtype=float) 
    mu = np.zeros((n, n))
    norms = np.zeros(n)
    for i in range(n):
        A_star[i] = A[i]
        for j in range(i):
            mu[i][j] = np.dot(A[i], A_star[j]) / np.dot(A_star[j], A_star[j])
            A_star[i] -= mu[i][j] * A_star[j]
        norms[i] = np.dot(A_star[i],A_star[i])
    return A_star, mu, norms

def lll_reduction(A, delta=0.75):
    A = np.array(A, dtype=float)
    n = len(A)
    A_star, mu, _ = gram_schmidt(A)

    k = 1
    while k < n:
        for j in range(k-1, -1, -1):
            q = round(mu[k][j])
            if q != 0:
                A[k] -= q*A[j]
                A_star, mu, _ = gram_schmidt(A)
        if np.dot(A_star[k], A_star[k]) >= (delta - mu[k][k-1] ** 2) * np.dot(A_star[k-1], A_star[k-1]):
            k += 1
        else:
            A[[k, k-1]] = A[[k-1, k]]
            A_star, mu, _ = gram_schmidt(A)
            k = max(k-1, 1)
    return np.array(A, dtype=int)


def aks_reduction(B, num_samples=1000, sigma=3.0):
    B = np.array(B)
    n = B.shape[1]
    shortest_vector = None
    min_norm = float('inf')

    for _ in range(num_samples):
        v = np.random.normal(loc=0, scale=sigma, size=n)
        coeffs = np.linalg.lstsq(B.T, v, rcond=None)[0]
        x = np.dot(B.T, np.round(coeffs))
        norm = np.linalg.norm(x)

        if 0 < norm < min_norm: 
            min_norm = norm
            shortest_vector = x

    return shortest_vector.astype(int), min_norm


def random_lattice_vector(B, coeff_range=5):
    n = B.shape[0]
    coeffs = np.random.randint(-coeff_range, coeff_range + 1, size=n)
    return np.dot(coeffs, B)

def hill_climb_svp(B, start_vec=None, k=1, max_steps=1000):
    n = B.shape[0]
    if start_vec is None:
        v = B[0].copy()  
    else:
        v = start_vec.copy()
    
    step = 0
    while step < max_steps:
        improved = False
        for i in range(n):
            for c in range(-k, k + 1):
                if c == 0:
                    continue
                candidate = v + c * B[i]
                if np.linalg.norm(candidate) < np.linalg.norm(v):
                    v = candidate
                    improved = True
        if not improved:
            break
        step += 1
    return v


def enum_svp(B, radius=None):
    B = np.array(B, dtype=float)
    n = len(B)
    B_star, mu, norms = gram_schmidt(B)

    initial_radius = np.sum(np.linalg.norm(B, axis=1) ** 2)
    best = {'vec': None, 'norm2': initial_radius}

    coeffs = np.zeros(n)
    y = np.zeros(n)
    c = np.zeros(n)
    dist = np.zeros(n+1)
    dist[n] = 0

    def search(k):
        if k == -1:
            v = sum(int(round(coeffs[i])) * B[i] for i in range(n))
            norm2 = np.dot(v, v)
            if 0 < norm2 < best['norm2']:
                best['vec'] = v
                best['norm2'] = norm2
            return

        r = int(np.sqrt((radius if radius else best['norm2']) - dist[k+1]) / np.sqrt(norms[k])) + 1
        c[k] = 0
        for z in range(-r, r + 1):
            coeffs[k] = z
            y[k] = z - c[k]
            dist[k] = dist[k + 1] + y[k] ** 2 * norms[k]
            if dist[k] < (radius if radius else best['norm2']):
                search(k - 1)

    search(n - 1)
    return best['vec'], np.sqrt(best['norm2'])


def gauss_sieve(B, max_iter=1000, coeff_range=5):
    B = np.array(B, dtype=float)
    n = B.shape[0]
    L = []

    for _ in range(max_iter):
        v = random_lattice_vector(B, coeff_range)
        changed = True
        while changed:
            changed = False
            for i, w in enumerate(L):
                if np.linalg.norm(v - w) < np.linalg.norm(v):
                    v = v - w
                    changed = True
            if np.allclose(v, 0):
                break

        if not np.allclose(v, 0):
            # Also try to reduce existing list
            for i, w in enumerate(L):
                if np.linalg.norm(w - v) < np.linalg.norm(w):
                    L[i] = w - v
            L.append(v)

    if L:
        shortest = min(L, key=lambda x: np.linalg.norm(x))
        return shortest, np.linalg.norm(shortest)
    else:
        return None, None

def printing(A):
    print("A = ")
    print(A)
    print("\n")
    
    A_reduced = lll_reduction(A)
    print("=== LLL ===\n")
    print("Zredukowana baza (LLL):")
    print(A_reduced)
    print("\n")


    shortest_vec, norm = aks_reduction(A)
    print("=== Ajtai–Kumar–Sivakumar ===\n")
    print("Przybliżenie najkrótszego wektora:", shortest_vec)
    print("Norma:", norm)
    print("\n")


    print("=== Hill Climbing SVP ===\n")
    for i in range(10):
        v0 = random_lattice_vector(A, coeff_range=4)
        result = hill_climb_svp(A, start_vec=v0, k=1)
        print(f"Próba {i + 1}: Norma = {np.linalg.norm(result):.2f}, Wektor = {result}")
    print("\n")


    print("=== ENUM (enumeracja) ===\n")
    e_v, e_norm = enum_svp(A)
    print("Najkrótszy wektor:", e_v)
    print("Norma:", e_norm)
    print("\n")


    print("=== Gauss Sieve ===\n")
    g_v, g_norm = gauss_sieve(A)
    print("Najkrótszy wektor:", g_v)
    print("Norma:", g_norm)
    print("\n")


if __name__ == "__main__":
    A = np.array([[4, 1], [2, 3]])
    printing(A)
    B = np.array([
        [4, 1, 3, 2],
        [2, 1, 1, 0],
        [0, 2, 1, 1],
        [1, 0, 2, 3]
    ])
    printing(B)
    
