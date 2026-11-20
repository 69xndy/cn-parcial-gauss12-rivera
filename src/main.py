"""Punto de entrada único: reproduce todo el contenido de resultados/.

Uso:  python src/main.py
"""
import csv
import os
import platform
import sys

import numpy as np
from mpmath import mp, mpf, log, log10

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import conteo  # noqa: E402
import gauss  # noqa: E402
import verificacion as ver  # noqa: E402

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATOS = os.path.join(RAIZ, "datos")
RES = os.path.join(RAIZ, "resultados")
N = gauss.N
ESTRATEGIAS = ("ninguno", "parcial", "escalado")
NOMBRES = {"ninguno": "sin pivoteo", "parcial": "pivoteo parcial", "escalado": "pivoteo parcial escalado"}
PASADAS_REFINAMIENTO = 6
mp.dps = 60


def fmt(x, cifras=6):
    return mp.nstr(mpf(x), cifras)


def error_relativo(x, x_ex):
    """(norma-inf normwise, máximo componente a componente) del error relativo vs la solución exacta."""
    ex = [ver.a_mp(q) for q in x_ex]
    ap = [mpf(v) if not hasattr(v, "numerator") else ver.a_mp(v) for v in x]
    dif = [abs(ap[i] - ex[i]) for i in range(N)]
    return max(dif) / max(abs(t) for t in ex), max(dif[i] / abs(ex[i]) for i in range(N))


def bits_perdidos(err_rel):
    """Bits de la mantisa (53) que no sobreviven: 53 - bits correctos = 53 + log2(err)."""
    return 53 + log(err_rel, 2)


def escribir_traza(ruta, A, fact):
    def bloques(M, perm):
        out = []
        for c0 in (0, 4, 8):
            out.append("    fila  | " + " ".join(f"{'c' + str(j + 1):>24}" for j in range(c0, c0 + 4)))
            for i in range(N):
                vals = " ".join(f"{M[i][j]:+24.17e}" for j in range(c0, c0 + 4))
                out.append(f"    F{perm[i] + 1:<4} | {vals}")
            out.append("")
        return out

    L = ["TRAZA DE LA ELIMINACION (pivoteo parcial escalado, doble precision)",
         "Cada valor con 17 cifras significativas. 'F#' = fila ORIGINAL que ocupa esa posicion.",
         "Los ceros bajo la diagonal son asignados por construccion (no calculados).", ""]
    L.append("=== Paso 0: matriz original A ===")
    L += bloques(A, list(range(N)))
    for t in fact["traza"]:
        L.append(f"=== Paso {t['paso']} de {N - 1}: pivote = {t['pivote']:+.17e} "
                 f"(fila original F{t['fila_pivote_original'] + 1}); multiplicador maximo = {t['mult_max']:.6e} ===")
        L.append("    permutacion (filas originales en orden): " + " ".join(f"F{p + 1}" for p in t["perm"]))
        L += bloques(t["U"], t["perm"])
    bajo = sum(1 for i in range(N) for j in range(i) if fact["U"][i][j] == 0.0)
    L.append(f"Verificacion final: {bajo} de {N * (N - 1) // 2} elementos bajo la diagonal son exactamente 0.0")
    with open(ruta, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    return bajo


def main():
    os.makedirs(RES, exist_ok=True)
    A_txt = gauss.leer_matriz_texto(os.path.join(DATOS, "matriz_a.csv"))
    b_txt = gauss.leer_vector_texto(os.path.join(DATOS, "vector_b.csv"))
    assert len(A_txt) == N and all(len(f) == N for f in A_txt) and len(b_txt) == N
    A = gauss.a_float(A_txt)
    b = [float(v) for v in b_txt]
    A_f = ver.a_fracciones(A_txt)
    b_f = ver.vector_fracciones(b_txt)
    x_ex = ver.solucion_exacta(A_f, b_f)  # control exacto (racional)
    nA_inf = max(sum(abs(t) for t in fila) for fila in A_f)

    corridas, facts, sols, stats = {}, {}, {}, {}
    for est in ESTRATEGIAS:
        c = gauss.Contador()
        fact = gauss.factorizar(A, est, c, trazar=(est == "escalado"))
        x = gauss.resolver(fact, b, c)
        corridas[est], facts[est], sols[est] = c, fact, x
        e_norm, e_comp = error_relativo(x, x_ex)
        r = ver.residuo_exacto(A_f, b_f, x)
        stats[est] = {
            "mult_max": fact["mult_max"], "crecimiento": fact["crecimiento"],
            "err_norm": e_norm, "err_comp": e_comp, "bits": bits_perdidos(e_norm),
            "r_inf": ver.norma_inf(r),
            "eta": ver.a_mp(ver.norma_inf(r)) / (ver.a_mp(nA_inf) * max(abs(mpf(v)) for v in x)
                                                  + max(abs(ver.a_mp(t)) for t in b_f)),
            "digitos": [ver.digitos_correctos(x[i], x_ex[i]) for i in range(N)],
        }

    # ---------------- resultados/traza_eliminacion.txt ----------------
    bajo = escribir_traza(os.path.join(RES, "traza_eliminacion.txt"), A, facts["escalado"])
    assert bajo == N * (N - 1) // 2

    # ---------------- resultados/conteo_flops.md ----------------
    resumen_flops = conteo.escribir_reporte(os.path.join(RES, "conteo_flops.md"), corridas, N)

    # ---------------- refinamiento iterativo (sobre la corrida escalada) ----------------
    fact = facts["escalado"]
    hist = ver.refinar(fact, A_f, b_f, sols["escalado"], PASADAS_REFINAMIENTO)
    x0 = sols["escalado"]
    x1 = hist[1]["x"]
    delta1 = hist[1]["delta"]
    afirmados = []
    for i in range(N):
        rel = abs(ver.a_mp(delta1[i])) / abs(ver.a_mp(x1[i]))
        afirmados.append(int(mp.floor(-log10(rel))))
    reales = stats["escalado"]["digitos"]

    # ---------------- resultados/solucion.csv ----------------
    xf = hist[-1]["x"]
    dig_ref = []
    for i in range(N):
        d = abs(ver.a_mp(hist[-1]["delta"][i])) / abs(ver.a_mp(xf[i]))
        dig_ref.append(int(mp.floor(-log10(d))) if d != 0 else 40)
    DIG_REF_REPORTADOS = min(30, min(dig_ref))
    with open(os.path.join(RES, "solucion.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["i", "x_reportada", "digitos_afirmados", "x_doble_17_cifras",
                    "digitos_correctos_vs_exacta", "x_refinada", "digitos_refinada"])
        for i in range(N):
            w.writerow([i + 1, mp.nstr(mpf(x0[i]), afirmados[i]), afirmados[i],
                        f"{x0[i]:.17e}", reales[i],
                        mp.nstr(ver.a_mp(xf[i]), DIG_REF_REPORTADOS, strip_zeros=False), DIG_REF_REPORTADOS])

    # ---------------- κ(A) ----------------
    k1 = ver.estimar_kappa1(fact, A)
    k1_eq = ver.estimar_kappa1(fact, A, equilibrada=True)
    Af_np = np.array(A)
    s = np.abs(Af_np).max(axis=1)
    DA = Af_np / s[:, None]
    k1_np, k2_np = np.linalg.cond(Af_np, 1), np.linalg.cond(Af_np, 2)
    k1_eq_np, k2_eq_np = np.linalg.cond(DA, 1), np.linalg.cond(DA, 2)
    eps = 2.0 ** -52

    # ---------------- control con biblioteca (solo control) ----------------
    x_np = np.linalg.solve(Af_np, np.array(b))
    e_np, _ = error_relativo([float(v) for v in x_np], x_ex)

    # ---------------- resultados/residuo.txt ----------------
    L = ["RESIDUO r = b - A x EN PRECISION EXTENDIDA (fractions: aritmetica racional exacta)", ""]
    L.append("A y b se leen como decimales exactos; x (double) se convierte a fraccion sin perdida.")
    L.append("")
    L.append("Norma infinito del residuo por estrategia:")
    for est in ESTRATEGIAS:
        L.append(f"  {NOMBRES[est]:<28} ||r||_inf = {fmt(ver.a_mp(stats[est]['r_inf']), 6)}"
                 f"   error relativo vs exacta = {fmt(stats[est]['err_norm'], 6)}"
                 f"   error inverso normwise = {fmt(stats[est]['eta'], 3)}")
    L.append("")
    L.append("Componentes del residuo (pivoteo parcial escalado, antes de refinar), 12 cifras:")
    r0 = ver.residuo_exacto(A_f, b_f, x0)
    for i in range(N):
        L.append(f"  r[{i + 1:>2}] = {fmt(ver.a_mp(r0[i]), 12):>22}")
    L.append("")
    L.append("Refinamiento iterativo (r exacto; correccion con L y U reutilizadas, sin refactorizar):")
    L.append(f"  {'pasada':>6} {'||r||_inf':>14} {'||delta||_inf':>14} {'error rel. vs exacta':>22}")
    for h in hist:
        e_n, _ = error_relativo(h["x"], x_ex)
        dl = "-" if h["delta_inf"] is None else fmt(ver.a_mp(h["delta_inf"]), 4)
        L.append(f"  {h['iter']:>6} {fmt(ver.a_mp(h['r_inf']), 4):>14} {dl:>14} {fmt(e_n, 4):>22}")
    L.append("")
    L.append(f"Control con numpy.linalg.solve (solo control): error relativo vs exacta = {fmt(e_np, 4)}")
    L.append("")
    L.append("ANEXO A. Estabilidad de las tres corridas (mismo codigo, misma factorizacion reutilizada)")
    L.append(f"  {'estrategia':<28} {'mult max':>11} {'mult max esc.':>13} {'crecimiento':>12} {'err norma':>10} {'err comp.':>10} {'bits perd.':>10} {'digitos (min-max)':>18}")
    for est in ESTRATEGIAS:
        s_ = stats[est]
        dg = s_["digitos"]
        L.append(f"  {NOMBRES[est]:<28} {s_['mult_max']:>11.3e} {facts[est]['mult_esc_max']:>13.3e} "
                 f"{s_['crecimiento']:>12.3e} {fmt(s_['err_norm'], 3):>10} {fmt(s_['err_comp'], 3):>10} "
                 f"{fmt(s_['bits'], 4):>10} {min(dg):>9}-{max(dg):<8}")
    L.append("  (mult max esc. = max |m(i,k)|*s_k/s_i: el multiplicador del sistema con filas equilibradas)")
    L.append("")
    L.append("ANEXO B. Condicionamiento")
    L.append(f"  kappa_1(A)   estimado (Hager, reutiliza L y U) = {k1:.4e}   control numpy = {k1_np:.4e}")
    L.append(f"  kappa_2(A)   control numpy                     = {k2_np:.4e}")
    L.append(f"  kappa_1(D A) estimado (Hager)                  = {k1_eq:.4e}   control numpy = {k1_eq_np:.4e}")
    L.append(f"  kappa_2(D A) control numpy                     = {k2_eq_np:.4e}   (D = diag(1/max_j |a_ij|))")
    L.append(f"  kappa_2(A) * eps = {k2_np * eps:.3e}   kappa_1(A) * eps = {k1 * eps:.3e}   (eps = 2^-52)")
    L.append(f"  kappa_2(DA) * eps = {k2_eq_np * eps:.3e}   kappa_1(DA) * eps = {k1_eq * eps:.3e}")
    L.append("")
    L.append("ANEXO C. Casi dependencia lineal (filas normalizadas a norma 2 = 1, ajuste por minimos cuadrados)")
    for (i, a_, b_) in ((6, 2, 4), (4, 2, 6), (2, 4, 6)):
        al, be, rr = ver.casi_dependencia(A, i, a_, b_)
        L.append(f"  F{i + 1} ~ {al:+.4f}*F{a_ + 1} {be:+.4f}*F{b_ + 1}   residuo relativo = {rr:.3e}")
    sv = np.linalg.svd(DA, compute_uv=False)
    L.append(f"  Valores singulares de D A (control numpy): mayor = {sv[0]:.4e}, menor = {sv[-1]:.4e}")
    with open(os.path.join(RES, "residuo.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")

    # ---------------- resumen por pantalla ----------------
    print("== Costo (nucleo):", resumen_flops["nucleo"], "flops;",
          "mult/div", resumen_flops["muldiv"], "sumas/restas", resumen_flops["addsub"])
    for est in ESTRATEGIAS:
        s_ = stats[est]
        c = corridas[est]
        print(f"== {NOMBRES[est]:<26} mult_max={s_['mult_max']:.3e} crecimiento={s_['crecimiento']:.3e} "
              f"err_rel={fmt(s_['err_norm'], 4)} (comp {fmt(s_['err_comp'], 4)}) bits_perdidos={fmt(s_['bits'], 4)} "
              f"||r||={fmt(ver.a_mp(s_['r_inf']), 4)} flops={c.flops_nucleo()} comp={c.total('comparaciones')} "
              f"divesc={c.total('divisiones_escalado')} swaps={c.total('intercambios')}")
        print("   digitos correctos:", s_["digitos"])
    print("== digitos afirmados (escalado):", afirmados)
    print("== digitos refinada:", dig_ref)
    print(f"== kappa1 est={k1:.4e} (numpy {k1_np:.4e});  kappa2 numpy={k2_np:.4e}")
    print(f"== equilibrada: kappa1 est={k1_eq:.4e} (numpy {k1_eq_np:.4e}); kappa2 numpy={k2_eq_np:.4e}")
    print(f"== kappa2*eps = {k2_np * eps:.3e}; kappa1*eps = {k1 * eps:.3e}")
    print("== numpy.solve err_rel:", fmt(e_np, 4))
    print("== versiones: python", platform.python_version(), "numpy", np.__version__, "mpmath", mp.__version__ if hasattr(mp, '__version__') else '')
    for h in hist:
        e_n, _ = error_relativo(h["x"], x_ex)
        print("   iter", h["iter"], "r_inf", fmt(ver.a_mp(h["r_inf"]), 4), "delta", None if h["delta_inf"] is None else fmt(ver.a_mp(h["delta_inf"]), 4), "err", fmt(e_n, 4))


if __name__ == "__main__":
    main()
