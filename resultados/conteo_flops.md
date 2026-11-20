# Conteo de operaciones (n = 12)

Generado por `src/main.py` a partir de los contadores de `src/gauss.py`.

## Convención de conteo

- 1 flop = 1 operación aritmética elemental (+, −, ×, ÷) en doble precisión.
- Los 66 elementos bajo la diagonal se anulan por construcción (`U[i][k] = 0.0`): no se calculan y no cuestan flops.
- Un multiplicador `m(i,k) = a(i,k)/a(k,k)` por fila y por paso, calculado una sola vez.
- Se factoriza una sola vez (P·A = L·U); b se actualiza con los mismos multiplicadores.
- Comparaciones, divisiones de escalado e intercambios de filas se reportan aparte (sobrecosto de estabilidad).

## Núcleo: eliminación densa paso a paso (pivoteo parcial escalado)

| Paso | mult/div (A) | sumas/restas (A) | mult/div (b) | sumas/restas (b) | Total |
|---|---:|---:|---:|---:|---:|
| 1 | 132 | 121 | 11 | 11 | 275 |
| 2 | 110 | 100 | 10 | 10 | 230 |
| 3 | 90 | 81 | 9 | 9 | 189 |
| 4 | 72 | 64 | 8 | 8 | 152 |
| 5 | 56 | 49 | 7 | 7 | 119 |
| 6 | 42 | 36 | 6 | 6 | 90 |
| 7 | 30 | 25 | 5 | 5 | 65 |
| 8 | 20 | 16 | 4 | 4 | 44 |
| 9 | 12 | 9 | 3 | 3 | 27 |
| 10 | 6 | 4 | 2 | 2 | 14 |
| 11 | 2 | 1 | 1 | 1 | 5 |
| sust. regresiva | 78 | 66 | 0 | 0 | 144 |
| **Total** | **650** | **572** | **66** | **66** | **1354** |

En la fila «sust. regresiva» las mult/div incluyen las 12 divisiones por la diagonal de U.

Resumen: multiplicaciones/divisiones = **716**, sumas/restas = **638**, total del núcleo = **1354** flops.

## Contraste con la fórmula cerrada y con la cota de referencia

| Concepto | Valor |
|---|---:|
| Contador (núcleo) | 1354 |
| Fórmula cerrada Σ m(2m+1) + Σ 2m + n² = 1078 + 132 + 144 | 1354 |
| 2n³/3 (término dominante) | 1152 |
| Cota de referencia del enunciado | 1486 |
| Desviación respecto de la cota | -8.9 % (aceptado: ±10 %) |

La cota 1486 se reproduce exactamente sumando 132 flops = 2 × 66: son las operaciones que resultan de *operar* la columna k en cada paso, es decir, calcular unos elementos que ya se sabe que valen cero. Aquí no se calculan; por eso el contador da 1354 y la cota se cumple con margen (1354 ≤ 1486).

## Sobrecosto de estabilidad (aparte del núcleo)

| Estrategia | Flops núcleo | Comparaciones | Divisiones de escalado | Intercambios | Flops totales (núcleo + divisiones) |
|---|---:|---:|---:|---:|---:|
| Sin pivoteo | 1354 | 0 | 0 | 0 | 1354 |
| Pivoteo parcial | 1354 | 66 | 0 | 9 | 1354 |
| Pivoteo parcial escalado | 1354 | 198 | 77 | 8 | 1431 |

Desglose del pivoteo escalado: 132 comparaciones para los factores de escala s_i = max_j |a_ij| (12 filas × 11), 66 comparaciones en la búsqueda de pivote (11 + 10 + … + 1) y 77 divisiones |a_ik|/s_i. El pivoteo parcial simple solo necesita 66 comparaciones (cota del enunciado: 1486 + 66).

Si se suman TODAS las operaciones del pivoteo escalado (núcleo + divisiones + comparaciones) el total es 1629, frente a 1552 de la cota con pivoteo parcial: +5.0 %.

## Costo de reutilizar L y U

Cada resolución adicional (por ejemplo, una pasada de refinamiento iterativo) cuesta solo 276 flops (sustitución progresiva + regresiva), sin refactorizar.

