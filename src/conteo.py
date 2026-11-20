"""Instrumentación y reporte del costo operacional.

Todos los números salen de los contadores incrementados dentro de gauss.py,
operación por operación. Las fórmulas cerradas de este módulo solo sirven para
CONTRASTAR el contador, nunca para producirlo.
"""

REFERENCIA_MINIMA = 1486          # cota del enunciado (incluye operar la columna k)
REFERENCIA_PIVOTEO_COMP = 66      # comparaciones extra de la cota con pivoteo parcial


def formula_cerrada(n):
    """Flops del núcleo si NO se opera la columna k: eliminación + b + sustitución regresiva."""
    elim_A = sum(m * (2 * m + 1) for m in range(1, n))     # m multiplicadores + 2m por fila
    act_b = sum(2 * m for m in range(1, n))
    sust = n * n
    return elim_A + act_b + sust, elim_A, act_b, sust


def tabla_pasos(contador, n):
    """Filas de la tabla paso a paso (11 pasos) + sustitución regresiva."""
    filas = []
    for k in range(1, n):
        a = contador.datos.get(("A", k), {})
        b = contador.datos.get(("b", k), {})
        md_a, as_a = a.get("muldiv", 0), a.get("addsub", 0)
        md_b, as_b = b.get("muldiv", 0), b.get("addsub", 0)
        filas.append((k, md_a, as_a, md_b, as_b, md_a + as_a + md_b + as_b))
    s = contador.datos.get(("sust", 0), {})
    filas.append(("sust. regresiva", s.get("muldiv", 0), s.get("addsub", 0), 0, 0,
                  s.get("muldiv", 0) + s.get("addsub", 0)))
    return filas


def escribir_reporte(ruta, corridas, n):
    """corridas: dict nombre -> Contador. Escribe resultados/conteo_flops.md."""
    esc = corridas["escalado"]
    total_nucleo = esc.flops_nucleo()
    cerrada, elim_A, act_b, sust = formula_cerrada(n)
    L = []
    L.append("# Conteo de operaciones (n = 12)\n")
    L.append("Generado por `src/main.py` a partir de los contadores de `src/gauss.py`.\n")
    L.append("## Convención de conteo\n")
    L.append("- 1 flop = 1 operación aritmética elemental (+, −, ×, ÷) en doble precisión.")
    L.append("- Los 66 elementos bajo la diagonal se anulan por construcción (`U[i][k] = 0.0`): "
             "no se calculan y no cuestan flops.")
    L.append("- Un multiplicador `m(i,k) = a(i,k)/a(k,k)` por fila y por paso, calculado una sola vez.")
    L.append("- Se factoriza una sola vez (P·A = L·U); b se actualiza con los mismos multiplicadores.")
    L.append("- Comparaciones, divisiones de escalado e intercambios de filas se reportan aparte "
             "(sobrecosto de estabilidad).\n")
    L.append("## Núcleo: eliminación densa paso a paso (pivoteo parcial escalado)\n")
    L.append("| Paso | mult/div (A) | sumas/restas (A) | mult/div (b) | sumas/restas (b) | Total |")
    L.append("|---|---:|---:|---:|---:|---:|")
    filas = tabla_pasos(esc, n)
    tot = [0, 0, 0, 0, 0]
    for f in filas:
        L.append(f"| {f[0]} | {f[1]} | {f[2]} | {f[3]} | {f[4]} | {f[5]} |")
        for i in range(5):
            tot[i] += f[i + 1]
    L.append(f"| **Total** | **{tot[0]}** | **{tot[1]}** | **{tot[2]}** | **{tot[3]}** | **{tot[4]}** |")
    L.append("")
    L.append("En la fila «sust. regresiva» las mult/div incluyen las 12 divisiones por la diagonal de U.\n")
    md = tot[0] + tot[2]
    ad = tot[1] + tot[3]
    L.append(f"Resumen: multiplicaciones/divisiones = **{md}**, sumas/restas = **{ad}**, "
             f"total del núcleo = **{total_nucleo}** flops.\n")
    L.append("## Contraste con la fórmula cerrada y con la cota de referencia\n")
    L.append("| Concepto | Valor |")
    L.append("|---|---:|")
    L.append(f"| Contador (núcleo) | {total_nucleo} |")
    L.append(f"| Fórmula cerrada Σ m(2m+1) + Σ 2m + n² = {elim_A} + {act_b} + {sust} | {cerrada} |")
    L.append(f"| 2n³/3 (término dominante) | {2 * n**3 / 3:.0f} |")
    L.append(f"| Cota de referencia del enunciado | {REFERENCIA_MINIMA} |")
    desv = 100.0 * (total_nucleo - REFERENCIA_MINIMA) / REFERENCIA_MINIMA
    L.append(f"| Desviación respecto de la cota | {desv:+.1f} % (aceptado: ±10 %) |")
    L.append("")
    extra = 2 * (n * (n - 1) // 2)
    L.append(f"La cota 1486 se reproduce exactamente sumando {extra} flops = 2 × 66: son las operaciones "
             f"que resultan de *operar* la columna k en cada paso, es decir, calcular unos elementos "
             f"que ya se sabe que valen cero. Aquí no se calculan; por eso el contador da "
             f"{total_nucleo} y la cota se cumple con margen ({total_nucleo} ≤ {REFERENCIA_MINIMA}).\n")

    L.append("## Sobrecosto de estabilidad (aparte del núcleo)\n")
    L.append("| Estrategia | Flops núcleo | Comparaciones | Divisiones de escalado | Intercambios | Flops totales (núcleo + divisiones) |")
    L.append("|---|---:|---:|---:|---:|---:|")
    nombres = {"ninguno": "Sin pivoteo", "parcial": "Pivoteo parcial", "escalado": "Pivoteo parcial escalado"}
    for clave in ("ninguno", "parcial", "escalado"):
        c = corridas[clave]
        L.append(f"| {nombres[clave]} | {c.flops_nucleo()} | {c.total('comparaciones')} | "
                 f"{c.total('divisiones_escalado')} | {c.total('intercambios')} | "
                 f"{c.flops_nucleo() + c.sobrecosto_flops()} |")
    L.append("")
    pre = esc.datos.get(("piv", 0), {}).get("comparaciones", 0)
    L.append(f"Desglose del pivoteo escalado: {pre} comparaciones para los factores de escala "
             f"s_i = max_j |a_ij| (12 filas × 11), "
             f"{esc.total('comparaciones') - pre} comparaciones en la búsqueda de pivote "
             f"(11 + 10 + … + 1) y {esc.total('divisiones_escalado')} divisiones |a_ik|/s_i. "
             f"El pivoteo parcial simple solo necesita {corridas['parcial'].total('comparaciones')} comparaciones "
             f"(cota del enunciado: 1486 + {REFERENCIA_PIVOTEO_COMP}).\n")
    tot_con_todo = (total_nucleo + esc.sobrecosto_flops() + esc.total("comparaciones"))
    ref_p = REFERENCIA_MINIMA + REFERENCIA_PIVOTEO_COMP
    L.append(f"Si se suman TODAS las operaciones del pivoteo escalado (núcleo + divisiones + comparaciones) "
             f"el total es {tot_con_todo}, frente a {ref_p} de la cota con pivoteo parcial: "
             f"{100.0 * (tot_con_todo - ref_p) / ref_p:+.1f} %.\n")
    L.append("## Costo de reutilizar L y U\n")
    L.append("Cada resolución adicional (por ejemplo, una pasada de refinamiento iterativo) cuesta solo "
             f"{act_b + sust} flops (sustitución progresiva + regresiva), sin refactorizar.\n")
    with open(ruta, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    return {"nucleo": total_nucleo, "muldiv": md, "addsub": ad, "cerrada": cerrada}
