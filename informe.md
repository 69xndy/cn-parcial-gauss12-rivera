# El pivote — Eliminación gaussiana con mínimo costo operacional y máxima exactitud de mantisa

## 1. Identificación

| Campo | Detalle |
|---|---|
| Estudiante | Andy Rivera Olguin |
| Sección | 2 |
| Asignatura | Computación Numérica (INF-321) · Ingeniería Civil Informática · UCM |
| Docentes | Juan Manuel Cabezas Contardo · Pedro Pesante |
| Fecha | 20/11/2026 |
| Repositorio | https://github.com/69xndy/cn-parcial-gauss12-rivera (`cn-parcial-gauss12-rivera`) |

**Declaración de uso de asistentes de IA.** Se utilizaron herramientas de inteligencia artificial como apoyo durante el desarrollo del proyecto, principalmente para orientar la redacción y organización del informe y del archivo README, así como para revisar, corregir y mejorar algunos fragmentos del código fuente.

## 2. Planteamiento

Se resuelve el sistema denso **A x = b** con n = 12 (144 elementos no nulos, 9 cifras significativas por elemento; b con 20 cifras). Los datos se ingresan sin truncar: A y b se leen como texto y se convierten a `float` para la factorización en doble precisión y a `Fraction` (decimal exacto) para el residuo y los controles.

Tres rasgos hacen costoso el sistema (medidos por `main.py`, ver `resultados/residuo.txt`):

- **Escalas de fila incompatibles.** Las filas conviven entre 10⁻⁵ y 10⁵. Esto hace que κ₂(A) valga **3.08·10¹³**; si cada fila se normaliza por su máximo (matriz D·A), κ₂(D·A) cae a **8.13·10³**: nueve órdenes de magnitud del condicionamiento vienen de las escalas, no de la geometría del sistema.
- **Pivote inicial minúsculo.** a(1,1) = 1.3407·10⁻⁸. Sin pivoteo el multiplicador máximo es **2.75·10¹³** (≈ 2^44.6) y el elemento más grande de la matriz crece por un factor **2.64·10³** (≈ 2^11.4).
- **Casi dependencia lineal oculta entre las filas F3, F5 y F7.** Con las filas normalizadas a norma 2 = 1, F7 se reproduce con residuo relativo de solo **2.3·10⁻³** como combinación de F3 y F5 (F7 ≈ 0.989·F3 + 1.140·F5; el ajuste simétrico de F5 y F3 con las otras dos da residuos similares). El menor valor singular de D·A es **4.83·10⁻⁴** (el mayor, 3.93), que explica el κ₂(D·A) ≈ 8.1·10³ restante.

## 3. Estrategia de eliminación

Se factoriza **una sola vez** P·A = L·U (`gauss.factorizar`) y todo lo demás reutiliza L y U (`gauss.resolver`):

1. **Paso k (k = 1…11):** se elige el pivote; para cada fila i > k se calcula **un solo multiplicador** m(i,k) = a(i,k)/a(k,k) y se reutiliza en toda la fila; el barrido de la fila va de la columna k+1 en adelante y el elemento a(i,k) se **asigna 0.0 sin calcularlo**. Nunca se vuelve a tocar la parte ya anulada.
2. **Lado derecho:** b se actualiza con los mismos multiplicadores (sustitución progresiva por columnas, equivalente a eliminar sobre [A|b]).
3. **Sustitución regresiva** sobre U.
4. **No se calcula la inversa.** El estimador de condicionamiento (Hager, 1984; Higham, 1988) y el refinamiento usan solo resoluciones con L y U; cada resolución adicional cuesta **276 flops**, sin refactorizar.

Esto alcanza el mínimo de la eliminación densa (≈ 2n³/3 flops más términos de orden n²), como se cuenta en la sección siguiente. Para la estabilidad se compararon tres estrategias con el mismo código: sin pivoteo, pivoteo parcial y **pivoteo parcial escalado** (Forsythe y Moler), que elige la fila con mayor |a(i,k)|/s_i, con s_i = max_j |a(i,j)| de la matriz original. La elegida es la escalada (secuencia de pivotes en `resultados/traza_eliminacion.txt`: F8, F2, F1, F6, F9, F11, F12, F5, F4, F10, F3, F7).

## 4. Conteo de operaciones

**Convención.** 1 flop = 1 operación aritmética elemental (+, −, ×, ÷) en doble precisión. Los 66 elementos bajo la diagonal se anulan por construcción y no cuestan flops. Comparaciones, divisiones de escalado e intercambios de filas se reportan aparte. Los números salen de contadores incrementados **dentro del código, operación por operación** (`Contador` en `gauss.py`); la fórmula cerrada de `conteo.py` solo contrasta el resultado.

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
| Sust. regresiva | 78 | 66 | 0 | 0 | 144 |
| **Total** | **650** | **572** | **66** | **66** | **1354** |

Total: 716 multiplicaciones/divisiones + 638 sumas/restas = **1354 flops**. La fórmula Σ m(2m+1) + Σ 2m + n² = 1078 + 132 + 144 da lo mismo. Frente a la cota de referencia (1486) la desviación es **−8.9 %** (dentro de ±10 %). La cota 1486 se reproduce exactamente sumando 2 × 66 = 132 flops, que son los que costaría *operar* la columna k en cada paso (calcular elementos que ya valen cero); aquí no se calculan.

**Sobrecosto de estabilidad (aparte):**

| Estrategia | Flops núcleo | Comparaciones | Divisiones de escalado | Intercambios |
|---|---:|---:|---:|---:|
| Sin pivoteo | 1354 | 0 | 0 | 0 |
| Pivoteo parcial | 1354 | 66 | 0 | 9 |
| Pivoteo parcial escalado | 1354 | 198 (132 de escalas + 66 de búsqueda) | 77 | 8 |

Sumando todo el pivoteo escalado resultan 1629 operaciones, +5.0 % sobre la cota con pivoteo (1486 + 66 = 1552). Detalle completo en `resultados/conteo_flops.md`.

## 5. Análisis de mantisa y pivoteo

Error relativo medido contra la solución exacta racional (control, `fractions`); «bits perdidos» = 53 + log₂(error relativo), es decir, bits de la mantisa que no sobreviven.

| Corrida | Mult. máx. | Mult. máx. escalado* | Factor de crecimiento | Error rel. (norma ∞) | Error rel. (peor componente) | Bits perdidos (norma / componente) | Dígitos correctos (mín–máx) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Sin pivoteo | 2.753·10¹³ | 6.92·10³ | 2.641·10³ | 2.08·10⁻¹⁰ | 7.51·10⁻⁸ | 20.8 / 29.3 | 7–10 |
| Pivoteo parcial | 0.620 | 8.23 | 1.000 | 5.52·10⁻¹³ | 2.00·10⁻¹⁰ | 12.3 / 20.8 | 9–12 |
| Parcial escalado | 2.466·10⁹ | 0.962 | 1.783 | 2.36·10⁻¹³ | 8.54·10⁻¹¹ | 11.1 / 19.6 | 10–13 |

\* max |m(i,k)|·s_k/s_i: el multiplicador que tendría el sistema con filas equilibradas. El factor de crecimiento es máx_k máx_ij |a_ij⁽ᵏ⁾| / máx_ij |a_ij⁽⁰⁾|.

**Por qué el pivote a(1,1) es inaceptable.** Con a(1,1) = 1.34·10⁻⁸ el multiplicador de la fila 6 es 3.69·10⁵ / 1.34·10⁻⁸ ≈ 2.75·10¹³ ≈ 2^44.6. Al hacer a(i,j) − m·a(1,j), el término m·a(1,j) (del orden de 10⁷–10⁹) arrasa con el valor original de la fila i (10⁴–10⁶): al alinear exponentes, los bits bajos de a(i,j) quedan fuera de los 53 bits de mantisa y se pierden de forma irrecuperable. Ese crecimiento (×2641 ≈ 2^11.4) significa que los elementos más grandes llegan a ser unas 2 600 veces mayores que los datos y su error absoluto de redondeo crece en la misma proporción; esos errores entran a todas las etapas siguientes. Medido: se pierden **20.8 bits** (6.3 dígitos decimales) en norma y **29.3 bits** en el peor componente. Con pivoteo se evita elegir a(1,1) (en el paso 1 se elige F8, cuyo multiplicador máximo es 4.4), el crecimiento baja a 1.0–1.8 y la pérdida cae a **11–12 bits** en norma: el error disminuye en un factor ≈ 880 (2.08·10⁻¹⁰ → 2.36·10⁻¹³).

**Observaciones honestas.** (i) El enunciado anticipaba 6–7 dígitos sin pivoteo; aquí se obtuvieron 7–10 (el error en norma es 2.1·10⁻¹⁰, pero los componentes pequeños x₃ y x₆, de magnitud ≈ 10⁻³, pierden más cifras relativas). (ii) La ventaja del escalado sobre el pivoteo parcial simple es real pero modesta en este caso: ≈ 1.2 bits (error 5.5·10⁻¹³ → 2.4·10⁻¹³, un dígito más en varias incógnitas), a cambio de 77 divisiones y 132 comparaciones adicionales. (iii) El pivoteo escalado puede elegir pivotes pequeños en valor absoluto (p. ej. 8.8·10⁻⁵ de F1 en el paso 3) y tener multiplicadores grandes (2.5·10⁹); eso es inocuo porque en el sistema equilibrado el multiplicador máximo es 0.96 ≤ 1.

## 6. Condicionamiento

Estimaciones (`resultados/residuo.txt`, anexo B), con ε = 2⁻⁵² = 2.22·10⁻¹⁶:

| Magnitud | Valor |
|---|---:|
| κ₁(A) estimado con Hager (reutilizando L y U; sin inversa) | 4.4675·10¹³ (control numpy: igual) |
| κ₂(A) (control numpy) | 3.0847·10¹³ |
| κ₂(A)·ε | 6.85·10⁻³ |
| κ₁(D·A) estimado con Hager | 3.3967·10⁴ (control numpy: igual) |
| κ₂(D·A) (control numpy) | 8.132·10³ |
| κ₂(D·A)·ε | 1.81·10⁻¹² |

**Conclusión sobre dígitos.** La cota clásica ‖δx‖/‖x‖ ≲ κ₂(A)·ε ≈ 6.8·10⁻³ no garantiza ni dos dígitos. Sin embargo, la corrida con pivoteo escalado dio un error en norma de 2.4·10⁻¹³ (**12.6 dígitos**).

**Por qué la cota es pesimista aquí.** (1) κ(A) mide la sensibilidad ante perturbaciones arbitrarias de A y b en norma, pero buena parte de su tamaño es artificio de escala: multiplicar filas por constantes cambia κ(A) pero no la solución. Con filas equilibradas κ₂(D·A) = 8.1·10³, y la cota κ₂(D·A)·ε = 1.8·10⁻¹² (11.7 dígitos) **sí** es coherente con lo observado (2.4·10⁻¹³ < 1.8·10⁻¹²). El pivoteo escalado hace que el algoritmo se comporte como si resolviera el sistema equilibrado. (2) La cota supone el peor caso, con la perturbación alineada con el menor valor singular; es razonable que el error real de redondeo no se alinee con la dirección de la casi dependencia F3–F5–F7. (3) El residuo es quien decide: el error inverso normwise de la solución escalada es 1.1·10⁻¹⁷ (< u = 1.1·10⁻¹⁶), es decir, x̂ resuelve exactamente un sistema A + δA, con δA al nivel del redondeo.

**Matiz por componente.** La cota en norma habla de ‖δx‖/‖x‖ y las incógnitas pequeñas (x₃ = 0.0044, x₆ = −0.0063, frente a ‖x‖ = 2.85) heredan un error relativo propio unas 450–650 veces mayor; por eso x₃ y x₆ tienen 10 dígitos y x₁ tiene 13.

## 7. Resultados

**Criterio de «dígitos correctos»:** d = ⌊−log₁₀(|x̂ᵢ − xᵢ| / |xᵢ|)⌋, es decir, el error es menor que una unidad en la d-ésima cifra significativa. **Los dígitos afirmados no se toman de la solución exacta** (que solo se usa como control): se estiman con el tamaño de la primera corrección del refinamiento, δ₁ ≈ x̂ − x, y se verificó después que coinciden con los dígitos efectivamente correctos en las 12 incógnitas.

| i | x reportada (doble precisión, d cifras) | d afirmados | d reales (vs exacta) | x refinada (30 cifras) |
|---:|---:|---:|---:|---:|
| 1 | 0.7384915626999 | 13 | 13 | 0.738491562700000000000000000000 |
| 2 | -1.2059374168 | 12 | 12 | -1.20593741680000000000000000000 |
| 3 | 0.0043871952 | 10 | 10 | 0.00438719520000000000000000000000 |
| 4 | -0.0917264538 | 11 | 11 | -0.0917264538000000000000000000000 |
| 5 | 2.8461093725 | 12 | 12 | 2.84610937250000000000000000000 |
| 6 | -0.006281540899 | 10 | 10 | -0.00628154090000000000000000000000 |
| 7 | 1.4738205916 | 12 | 12 | 1.47382059160000000000000000000 |
| 8 | -0.0384927164 | 11 | 11 | -0.0384927164000000000000000000000 |
| 9 | 0.9152648037 | 12 | 12 | 0.915264803700000000000000000000 |
| 10 | -2.6071394852 | 12 | 12 | -2.60713948520000000000000000000 |
| 11 | 0.0518736294 | 11 | 11 | 0.0518736294000000000000000000000 |
| 12 | -1.0849371625 | 12 | 12 | -1.08493716250000000000000000000 |

Las 12 incógnitas de la corrida con pivoteo escalado tienen entre **10 y 13 dígitos correctos** (C4: ≥ 8). La solución refinada (columna «x refinada») coincide con la solución exacta del sistema dado con un error relativo menor que 10⁻⁶⁰ (resolución de los 60 dígitos de trabajo); resulta ser un decimal de 10 cifras significativas (a 30 cifras el resto son ceros), lo que es coherente con que b se construyó con aritmética racional exacta a partir de esa solución.

## 8. Verificación

- **Residuo en precisión extendida** (`fractions`, aritmética racional exacta, sin redondeo en r): ‖r‖∞ = **2.291·10⁻¹⁰** con pivoteo escalado (la componente dominante es r₆ = 2.29·10⁻¹⁰, en la fila de escala 10⁵–10⁶, es decir, 2.4·10⁻¹⁶ en términos relativos a su fila); 5.373·10⁻¹⁰ con pivoteo parcial y 1.772·10⁻⁷ sin pivoteo. Un residuo pequeño no implica error pequeño (κ es grande); por eso el residuo se complementa con el error contra la solución exacta.
- **Refinamiento iterativo** con residuo exacto y corrección calculada con L y U reutilizadas (276 flops por pasada):

| Pasada | ‖r‖∞ | ‖δ‖∞ | Error relativo vs exacta |
|---:|---:|---:|---:|
| 0 | 2.291·10⁻¹⁰ | — | 2.357·10⁻¹³ |
| 1 | 1.366·10⁻²² | 6.709·10⁻¹³ | 6.174·10⁻²⁶ |
| 2 | 3.928·10⁻³⁵ | 1.757·10⁻²⁵ | 1.614·10⁻³⁸ |
| 3 | 1.383·10⁻⁴⁸ | 4.593·10⁻³⁸ | 4.224·10⁻⁵¹ |
| 4 | 5.374·10⁻⁶¹ | 1.202·10⁻⁵⁰ | 1.068·10⁻⁶³ |
| 5 | 6.298·10⁻⁷³ | 3.14·10⁻⁶³ | < 10⁻⁶⁰ (resolución de trabajo) |
| 6 | 7.963·10⁻⁸⁶ | 8.22·10⁻⁷⁶ | < 10⁻⁶⁰ (resolución de trabajo) |

  Cada pasada gana ≈ 12.6 dígitos (la contracción del error es ≈ 2.6·10⁻¹³ por pasada), coherente con κ₂(D·A)·ε ≪ 1. El refinamiento **sí mejora**: de 12–13 a más de 60 dígitos.
- **Control con biblioteca** (solo control; no sustituye la implementación propia): `numpy.linalg.solve` da un error relativo de 5.23·10⁻¹³, del mismo orden que el pivoteo parcial propio (5.52·10⁻¹³) y 2.2 veces peor que el pivoteo escalado propio (2.36·10⁻¹³).

## 9. Conclusiones

- **Compromiso elegido.** El núcleo de la eliminación cuesta 1354 flops con cualquier estrategia; no hay trabajo redundante que ahorrar (un multiplicador por fila, una sola factorización, sin inversa). El pivoteo escalado añade 198 comparaciones, 77 divisiones y 8 intercambios (≈ 20 % de operaciones extra frente al núcleo: 1629 vs. 1354 si se cuenta todo, +5.0 % sobre la cota de referencia con pivoteo) y reduce el error ≈ 880 veces respecto de no pivotear. Es el compromiso razonable: pagar un sobrecosto pequeño y conocido por 9.8 bits de mantisa. La estrategia más barata (sin pivoteo) no es la más estable, y aquí se verificó con cifras.
- **Lo que no compra el pivoteo.** La exactitud final está acotada por κ(D·A)·ε; para ir más allá de 13 dígitos hace falta el refinamiento con residuo en precisión extendida, que cuesta 276 flops por pasada más el residuo exacto.
- **Con n = 1000.** El costo 2n³/3 ≈ 6.7·10⁸ flops obliga a abandonar los bucles de Python y usar una LU por bloques de biblioteca (nivel 3 de BLAS), con pivoteo parcial (el escalado explícito pasa a ser equilibrar filas una vez, O(n²)) y estimación de κ con Hager en O(n²) por iteración. El residuo exacto con `Fraction` sería inviable (O(n²) operaciones con enteros enormes); se reemplazaría por residuo con productos compensados o *double-double*, y la factorización podría hacerse en menor precisión con refinamiento en mayor precisión. También habría que revisar la estructura (si la matriz fuera dispersa, un solver disperso).

**Referencias.** Golub, G. H. y Van Loan, C. F., *Matrix Computations*; Higham, N. J., *Accuracy and Stability of Numerical Algorithms*; Forsythe, G. E. y Moler, C. B., *Computer Solution of Linear Algebraic Systems* (pivoteo escalado); Hager, W. W. (1984), «Condition estimates», *SIAM J. Sci. Stat. Comput.*; Higham, N. J. (1988), «FORTRAN codes for estimating the one-norm of a real or complex matrix», *ACM TOMS*; enunciado de la prueba (INF-321, UCM).

## 10. Reproducibilidad
- Python 3.13.2; numpy 2.3.1; mpmath 1.3.0 (`requirements.txt`; `fractions` y `csv` son de la biblioteca estándar).
- Comando: `pip install -r requirements.txt && python src/main.py` (reproduce todo `resultados/` sin intervención manual; los archivos generados son idénticos entre corridas).
- No se usa ninguna rutina de biblioteca para resolver: `numpy.linalg.solve`, `numpy.linalg.cond` y `numpy.linalg.svd` se usan solo como control y diagnóstico en `main.py`.

