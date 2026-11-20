"""Verificación independiente en precisión extendida (fractions / mpmath).

* Los datos A y b se leen como decimales EXACTOS (Fraction), no como floats.
* El residuo r = b - A·x se calcula exacto en racionales: no hay redondeo en r.
* El refinamiento iterativo reutiliza L y U (sin refactorizar).
* La solución exacta racional se usa SOLO como control para contar dígitos correctos.
"""
from fractions import Fraction

from mpmath import mp, mpf, log10

import gauss

mp.dps = 60


def a_fracciones(matriz_texto):
    return [[Fraction(v) for v in fila] for fila in matriz_texto]


def vector_fracciones(vector_texto):
    return [Fraction(v) for v in vector_texto]


def a_mp(q):
    """Fraction -> mpf (60 dígitos)."""
    return mpf(q.numerator) / mpf(q.denominator)


def residuo_exacto(A_f, b_f, x):
    """r = b - A·x en aritmética racional exacta. x puede ser float o Fraction."""
    n = len(A_f)
    xf = [Fraction(v) for v in x]  # un float se convierte a Fraction sin pérdida
    return [b_f[i] - sum(A_f[i][j] * xf[j] for j in range(n)) for i in range(n)]


def norma_inf(v):
    return max(abs(t) for t in v)


def solucion_exacta(A_f, b_f):
    """Eliminación con Fraction (pivote no nulo): solución exacta del sistema dado. Control."""
    n = len(A_f)
    M = [A_f[i][:] + [b_f[i]] for i in range(n)]
    for k in range(n):
        p = next(i for i in range(k, n) if M[i][k] != 0)
        M[k], M[p] = M[p], M[k]
        for i in range(k + 1, n):
            m = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= m * M[k][j]
    x = [Fraction(0)] * n
    for i in range(n - 1, -1, -1):
        s = M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))
        x[i] = s / M[i][i]
    return x


def refinar(fact, A_f, b_f, x0, iteraciones=8):
    """Refinamiento iterativo: x_{k+1} = x_k + resolver(LU, r_k), r_k exacto.

    x se mantiene en Fraction; solo la corrección se calcula en doble precisión.
    Devuelve el historial (una entrada por pasada, más la del punto de partida).
    """
    x = [Fraction(v) for v in x0]
    historial = []
    for it in range(iteraciones + 1):
        r = residuo_exacto(A_f, b_f, x)
        entrada = {"iter": it, "x": x[:], "r_inf": norma_inf(r), "delta_inf": None, "delta": None}
        if it > 0:
            entrada["delta_inf"] = norma_inf(delta)
            entrada["delta"] = delta
        historial.append(entrada)
        if it == iteraciones:
            break
        d = gauss.resolver(fact, [float(t) for t in r])  # r redondeado a double solo aquí
        delta = [Fraction(t) for t in d]
        x = [x[i] + delta[i] for i in range(len(x))]
    # Reordenar: la corrección de la pasada it se guarda en la entrada it (ya hecho arriba)
    return historial


def norma1_matriz(A):
    n = len(A)
    return max(sum(abs(A[i][j]) for i in range(n)) for j in range(n))


def estimar_norma1_inversa(resolver_B, resolver_Bt, n, max_iter=5):
    """Estimador de Hager de ||B^-1||_1 usando solo resoluciones (no forma la inversa)."""
    x = [1.0 / n] * n
    est = 0.0
    for _ in range(max_iter):
        y = resolver_B(x)
        est = sum(abs(t) for t in y)
        xi = [1.0 if t >= 0 else -1.0 for t in y]
        z = resolver_Bt(xi)
        zmax = max(abs(t) for t in z)
        if zmax <= sum(z[i] * x[i] for i in range(n)):
            break
        j = max(range(n), key=lambda i: abs(z[i]))
        x = [0.0] * n
        x[j] = 1.0
    return est


def estimar_kappa1(fact, A, equilibrada=False):
    """κ1(A) ≈ ||A||_1 · est(||A^-1||_1), reutilizando L y U.

    Con equilibrada=True estima κ1(D·A), D = diag(1/max_j|a_ij|): la matriz con filas
    normalizadas. Como (D·A)^-1 = A^-1·D^-1 no se necesita otra factorización.
    """
    n = len(A)
    if not equilibrada:
        nA = norma1_matriz(A)
        est = estimar_norma1_inversa(
            lambda v: gauss.resolver(fact, v),
            lambda v: gauss.resolver_transpuesta(fact, v),
            n,
        )
        return nA * est
    s = [max(abs(v) for v in fila) for fila in A]
    DA = [[A[i][j] / s[i] for j in range(n)] for i in range(n)]
    # (DA) z = c  ->  A z = D^-1 c = s·c
    res_B = lambda c: gauss.resolver(fact, [s[i] * c[i] for i in range(n)])
    # (DA)ᵀ z = c ->  Aᵀ (D z) = c -> w = A^-T c ; z = D^-1 w = s·w
    def res_Bt(c):
        w = gauss.resolver_transpuesta(fact, c)
        return [s[i] * w[i] for i in range(n)]
    return norma1_matriz(DA) * estimar_norma1_inversa(res_B, res_Bt, n)


def digitos_correctos(x_aprox, x_exacto):
    """floor(-log10(|x_aprox - x_exacto| / |x_exacto|)); None si el error es exactamente 0."""
    err = abs(Fraction(x_aprox) - x_exacto)
    if err == 0:
        return None
    rel = a_mp(err) / abs(a_mp(x_exacto))
    return int(mp.floor(-log10(rel)))


def casi_dependencia(A, i_obj, i_a, i_b):
    """Ajusta la fila i_obj (normalizada a norma 2 unitaria) como α·fila_a + β·fila_b por
    mínimos cuadrados (ecuaciones normales 2x2) y devuelve (α, β, residuo relativo)."""
    n = len(A)

    def unit(i):
        nrm = sum(t * t for t in A[i]) ** 0.5
        return [t / nrm for t in A[i]]

    f, u, v = unit(i_obj), unit(i_a), unit(i_b)
    dot = lambda p, q: sum(p[j] * q[j] for j in range(n))
    a11, a12, a22 = dot(u, u), dot(u, v), dot(v, v)
    r1, r2 = dot(u, f), dot(v, f)
    det = a11 * a22 - a12 * a12
    alfa = (r1 * a22 - r2 * a12) / det
    beta = (a11 * r2 - a12 * r1) / det
    res = [f[j] - alfa * u[j] - beta * v[j] for j in range(n)]
    return alfa, beta, dot(res, res) ** 0.5
