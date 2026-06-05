# Optimización de Cross Docking — LogiFast CR

Proyecto desarrollado para el curso de Ingeniería Industrial, Universidad de Costa Rica — I Semestre 2026.

## Descripción

Este proyecto resuelve el problema de programación de camiones en un centro de distribución tipo *cross docking* con almacenamiento temporal, usando un modelo de Programación Entera Mixta (MIP) basado en:

> Yu, W., & Egbelu, P. J. (2008). Scheduling of inbound and outbound trucks in cross docking systems with temporary storage. *European Journal of Operational Research*, 184(1), 377–396.

El objetivo es encontrar el orden óptimo de atención de los camiones de entrada y salida para minimizar el tiempo total de operación (makespan).

## Archivos

| Archivo | Descripción |
|---|---|
| `app.py` | Aplicación web en Streamlit |
| `cross_docking_pulp.py` | Modelo MIP resuelto con Python y PuLP |
| `cross_docking.mod` | Modelo formal en lenguaje AMPL |
| `TS5.txt` | Instancia de datos de prueba |

## Cómo correr la aplicación

```bash
pip install -r requirements.txt
python -m streamlit run app.py
```

Luego subís el archivo `TS5.txt` (o cualquier archivo en el mismo formato) y la app resuelve el modelo y muestra los resultados.

## Formato del archivo de datos

```
i 5   o 3   n 8
r  2  1  6
r  2  2  6
...
s  1  1  75
s  1  2  12
...
```

Donde:
- `i` = número de camiones de entrada
- `o` = número de camiones de salida  
- `n` = número de tipos de producto
- `r camión producto cantidad` = camión de entrada con su carga
- `s camión producto cantidad` = camión de salida con su pedido
