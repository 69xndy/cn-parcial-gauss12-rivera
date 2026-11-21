# cn-parcial-gauss12-rivera

Resuelve un sistema lineal denso de 12 incógnitas por eliminación gaussiana propia (una sola factorización P·A = L·U), con conteo instrumentado de operaciones, pivoteo parcial escalado y verificación en precisión extendida.

**Autoría:** <Andy Rivera Olguin> · Computación Numérica (INF-321), Ingeniería Civil Informática, UCM · 

## Requisitos

- Python 3.10 o superior (probado con 3.13.2).
- Instalación: `pip install -r requirements.txt` (numpy 2.3.1 y mpmath 1.3.0).

## Ejecución

```bash
python src/main.py
```

Un solo comando, sin intervención manual. Genera en `resultados/`:

- `solucion.csv`: las 12 incógnitas con los dígitos afirmados y la solución refinada.
- `traza_eliminacion.txt`: la matriz tras cada uno de los 11 pasos.
- `conteo_flops.md`: tabla de operaciones por paso.
- `residuo.txt`: residuo r = b − A·x (precisión extendida), norma infinito, refinamiento y anexos de estabilidad y condicionamiento.

## Estructura de carpetas

```
datos/       A (matriz_a.csv) y b (vector_b.csv) exactos, sin truncar
src/         gauss.py (eliminación + contador), conteo.py, verificacion.py, main.py
resultados/  salidas generadas por main.py
informe.md   informe técnico y análisis de error
```

## Resultado resumido

| i | x (10 cifras significativas) |
|---:|---:|
| 1 | 0.7384915627 |
| 2 | -1.205937417 |
| 3 | 0.004387195200 |
| 4 | -0.09172645380 |
| 5 | 2.846109373 |
| 6 | -0.006281540900 |
| 7 | 1.473820592 |
| 8 | -0.03849271640 |
| 9 | 0.9152648037 |
| 10 | -2.607139485 |
| 11 | 0.05187362940 |
| 12 | -1.084937163 |

(Valores de la solución refinada, redondeados a 10 cifras significativas.)

- **Flops del núcleo:** 1354 (+ 198 comparaciones, 77 divisiones y 8 intercambios del pivoteo escalado).
- **‖r‖∞ en precisión extendida:** 2.29·10⁻¹⁰ antes de refinar; 1.37·10⁻²² tras una pasada de refinamiento.

## Informe

Ver [informe.md](informe.md).

## Licencia

Trabajo académico de la asignatura INF-321 (UCM). Sin licencia de reutilización.
