"""Eliminación gaussiana densa con contador de operaciones instrumentado.

Resuelve P·A = L·U una sola vez y reutiliza L y U para cualquier lado derecho.
Estrategias de pivoteo:
    'ninguno'  : sin intercambio de filas (pivote = a(k,k)).
    'parcial'  : pivote = mayor |a(i,k)| entre las filas i >= k.
    'escalado' : pivote = mayor |a(i,k)| / s_i, con s_i = max_j |a(i,j)| de la matriz ORIGINAL.

Convención de conteo (se declara también en resultados/conteo_flops.md):
    * 1 flop = 1 operación aritmética elemental (+, -, *, /) en doble precisión.
    * La anulación de a(i,k) bajo la diagonal NO se calcula: se asigna 0.0 por construcción.
    * Comparaciones, divisiones de escalado e intercambios de filas se cuentan APARTE
      (sobrecosto de estabilidad); no se mezclan con los flops del núcleo.
    * El cálculo de |.| y de máximos para diagnóstico (factor de crecimiento) no se cuenta.
"""
import csv

N = 12


class Contador:
    """Acumula operaciones por etapa. Etapas: ('A',k), ('b',k), ('sust',0), ('piv',k)."""

    CAMPOS = ("muldiv", "addsub", "comparaciones", "divisiones_escalado", "intercambios")

    def __init__(self):
        self.datos = {}

    def sumar(self, etapa, campo, cantidad=1):
        fila = self.datos.setdefault(etapa, dict.fromkeys(self.CAMPOS, 0))
        fila[campo] += cantidad

    def total(self, campo, tipos=None):
        """Suma de un campo sobre las etapas cuyo tipo esté en `tipos` (None = todas)."""
        return sum(
            fila[campo]
            for (tipo, _), fila in self.datos.items()
            if tipos is None or tipo in tipos
        )

    def flops_nucleo(self):
        """Flops de la eliminación densa: A, actualización de b y sustitución regresiva."""
        return self.total("muldiv", ("A", "b", "sust")) + self.total("addsub", ("A", "b", "sust"))

    def sobrecosto_flops(self):
        """Divisiones |a_ik|/s_i del pivoteo escalado (operaciones aritméticas extra)."""
        return self.total("divisiones_escalado")


def leer_matriz_texto(ruta):
    """Lee el CSV como TEXTO (sin truncar) para poder convertir a float o a Fraction."""
    with open(ruta, newline="", encoding="utf-8") as f:
        return [fila for fila in csv.reader(f)]


def leer_vector_texto(ruta):
    with open(ruta, encoding="utf-8") as f:
        return [linea.strip() for linea in f if linea.strip()]


def a_float(matriz_texto):
    return [[float(v) for v in fila] for fila in matriz_texto]


def factorizar(A, estrategia, contador, trazar=False):
    """Factoriza P·A = L·U. Devuelve un diccionario con U, L, perm y estadísticas."""
    n = len(A)
    U = [fila[:] for fila in A]
    L = [[0.0] * n for _ in range(n)]
    perm = list(range(n))
    traza = []
    max_inicial = max(abs(v) for fila in A for v in fila)  # diagnóstico
    max_visto = max_inicial
    mult_max = 0.0
    mult_esc_max = 0.0  # diagnóstico: máx |m(i,k)|·s_k/s_i (multiplicador del sistema con filas equilibradas)
    esc_diag = [max(abs(v) for v in fila) for fila in A]  # diagnóstico: no se cuenta como flop

    escala = None
    if estrategia == "escalado":
        # s_i = max_j |a_ij| de la matriz original: (n-1) comparaciones por fila.
        escala = []
        for i in range(n):
            s = abs(U[i][0])
            for j in range(1, n):
                contador.sumar(("piv", 0), "comparaciones")
                if abs(U[i][j]) > s:
                    s = abs(U[i][j])
            escala.append(s)

    for k in range(n - 1):
        paso = k + 1

        # ---- elección del pivote (sobrecosto de estabilidad) ----
        p = k
        if estrategia in ("parcial", "escalado"):
            mejor = abs(U[k][k])
            if estrategia == "escalado":
                mejor = mejor / escala[k]
                contador.sumar(("piv", paso), "divisiones_escalado")
            for i in range(k + 1, n):
                valor = abs(U[i][k])
                if estrategia == "escalado":
                    valor = valor / escala[i]
                    contador.sumar(("piv", paso), "divisiones_escalado")
                contador.sumar(("piv", paso), "comparaciones")
                if valor > mejor:
                    mejor = valor
                    p = i
        if p != k:
            U[k], U[p] = U[p], U[k]
            L[k], L[p] = L[p], L[k]  # multiplicadores ya calculados (columnas < k)
            perm[k], perm[p] = perm[p], perm[k]
            if escala is not None:
                escala[k], escala[p] = escala[p], escala[k]
            esc_diag[k], esc_diag[p] = esc_diag[p], esc_diag[k]
            contador.sumar(("piv", paso), "intercambios")

        pivote = U[k][k]
        if pivote == 0.0:
            raise ZeroDivisionError(f"pivote nulo en el paso {paso}")

        # ---- eliminación: un multiplicador por fila, barrido desde la columna k+1 ----
        mult_paso = 0.0
        for i in range(k + 1, n):
            m = U[i][k] / pivote
            contador.sumar(("A", paso), "muldiv")
            L[i][k] = m
            if abs(m) > mult_paso:
                mult_paso = abs(m)
            m_esc = abs(m) * esc_diag[k] / esc_diag[i]  # diagnóstico
            if m_esc > mult_esc_max:
                mult_esc_max = m_esc
            U[i][k] = 0.0  # anulado por construcción: no se calcula
            for j in range(k + 1, n):
                U[i][j] = U[i][j] - m * U[k][j]
                contador.sumar(("A", paso), "muldiv")
                contador.sumar(("A", paso), "addsub")
        mult_max = max(mult_max, mult_paso)
        max_visto = max(max_visto, max(abs(v) for fila in U for v in fila))  # diagnóstico

        if trazar:
            traza.append(
                {
                    "paso": paso,
                    "fila_pivote_original": perm[k],
                    "pivote": pivote,
                    "mult_max": mult_paso,
                    "perm": perm[:],
                    "U": [fila[:] for fila in U],
                }
            )

    if U[n - 1][n - 1] == 0.0:
        raise ZeroDivisionError("pivote nulo en el último paso")

    for i in range(n):
        L[i][i] = 1.0
    return {
        "U": U, "L": L, "perm": perm, "estrategia": estrategia,
        "mult_max": mult_max, "mult_esc_max": mult_esc_max, "crecimiento": max_visto / max_inicial,
        "max_inicial": max_inicial, "max_visto": max_visto, "traza": traza,
    }


def resolver(fact, b, contador=None):
    """Resuelve A·x = b con L y U ya calculadas (sin refactorizar).

    Sustitución progresiva por columnas: en el paso k se actualiza b con el
    multiplicador m(i,k), exactamente lo que haría la eliminación sobre [A|b].
    """
    U, L, perm = fact["U"], fact["L"], fact["perm"]
    n = len(U)
    y = [b[perm[i]] for i in range(n)]
    for k in range(n - 1):
        for i in range(k + 1, n):
            y[i] = y[i] - L[i][k] * y[k]
            if contador is not None:
                contador.sumar(("b", k + 1), "muldiv")
                contador.sumar(("b", k + 1), "addsub")
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = y[i]
        for j in range(i + 1, n):
            s = s - U[i][j] * x[j]
            if contador is not None:
                contador.sumar(("sust", 0), "muldiv")
                contador.sumar(("sust", 0), "addsub")
        x[i] = s / U[i][i]
        if contador is not None:
            contador.sumar(("sust", 0), "muldiv")
    return x


def resolver_transpuesta(fact, c):
    """Resuelve Aᵀ·z = c con la misma factorización (Aᵀ = Uᵀ Lᵀ P). Solo para estimar κ."""
    U, L, perm = fact["U"], fact["L"], fact["perm"]
    n = len(U)
    w = [0.0] * n
    for i in range(n):  # Uᵀ w = c  (triangular inferior)
        s = c[i]
        for j in range(i):
            s -= U[j][i] * w[j]
        w[i] = s / U[i][i]
    v = [0.0] * n
    for i in range(n - 1, -1, -1):  # Lᵀ v = w  (triangular superior, diagonal unitaria)
        s = w[i]
        for j in range(i + 1, n):
            s -= L[j][i] * v[j]
        v[i] = s
    z = [0.0] * n
    for i in range(n):  # z = Pᵀ v
        z[perm[i]] = v[i]
    return z
