# ============================================================
# Modelo MIP: Programación de Camiones en Cross Docking
# Basado en: Yu & Egbelu (2008)
# Universidad de Costa Rica - Ingeniería Industrial
# ============================================================

# ── Conjuntos ────────────────────────────────────────────────
param I >= 1;           # número de camiones de entrada
param J >= 1;           # número de camiones de salida
param K >= 1;           # número de tipos de producto

set IN  := 1..I;
set OUT := 1..J;
set PROD := 1..K;

# ── Parámetros ───────────────────────────────────────────────
param a{IN, PROD}  >= 0;   # unidades producto k en camión entrada i
param b{OUT, PROD} >= 0;   # unidades producto k requeridas por camión salida j
param p  := 1;             # minutos por unidad (carga/descarga)
param h  := 5;             # minutos por lote (traslado interno)
param sc := 10;            # minutos de cambio entre camiones
param M  := 10000;         # constante grande (Big-M)

# Tiempos de proceso por camión
param pi_i{i in IN}  := p * sum{k in PROD} a[i,k];
param pi_j{j in OUT} := p * sum{k in PROD} b[j,k];

# ── Variables de decisión ────────────────────────────────────
# Continuas no negativas
var x{IN, OUT, PROD} >= 0;  # unidades de prod k transferidas DIRECTAMENTE de i a j
var y{IN, OUT, PROD} >= 0;  # unidades de prod k transferidas vía ALMACENAMIENTO TEMPORAL de i a j
var A{IN}  >= 0;             # tiempo inicio descarga camión entrada i
var D{IN}  >= 0;             # tiempo fin descarga camión entrada i
var E{OUT} >= 0;             # tiempo inicio carga camión salida j
var F{OUT} >= 0;             # tiempo fin carga/salida camión salida j
var C      >= 0;             # makespan (tiempo total de operación)

# Enteras binarias
var z{IN, OUT}   binary;    # 1 si hay transferencia (directa o temporal) de i a j
var delta{i in IN,  ii in IN  : i <> ii} binary;  # 1 si camión entrada i precede a i'
var sigma{j in OUT, jj in OUT : j <> jj} binary;  # 1 si camión salida j precede a j'

# ── Función objetivo ─────────────────────────────────────────
minimize Makespan: C;

# ── Restricciones ────────────────────────────────────────────

# (R1) El makespan es mayor o igual al tiempo de salida del último camión de salida
subject to R1 {j in OUT}:
    C >= F[j];

# (R2) Conservación de oferta: todo lo que llega en camión i se distribuye
subject to R2 {i in IN, k in PROD: a[i,k] > 0}:
    sum{j in OUT} (x[i,j,k] + y[i,j,k]) = a[i,k];

# (R3) Satisfacción de demanda: el camión j recibe exactamente lo que necesita
subject to R3 {j in OUT, k in PROD: b[j,k] > 0}:
    sum{i in IN} (x[i,j,k] + y[i,j,k]) = b[j,k];

# (R4) Relación entre variables de flujo y variable binaria z
subject to R4 {i in IN, j in OUT}:
    sum{k in PROD} (x[i,j,k] + y[i,j,k]) <= M * z[i,j];

# (R5) Tiempo de proceso del camión de entrada i
subject to R5 {i in IN}:
    D[i] = A[i] + pi_i[i];

# (R6) Tiempo de proceso del camión de salida j
subject to R6 {j in OUT}:
    F[j] = E[j] + pi_j[j];

# (R7) Secuenciación de camiones de entrada: si i precede a i', i' empieza después que i termina
subject to R7 {i in IN, ii in IN: i <> ii}:
    A[ii] >= D[i] + sc - M * (1 - delta[i,ii]);

# (R8) Orden total para camiones de entrada: exactamente uno precede al otro
subject to R8 {i in IN, ii in IN: i < ii}:
    delta[i,ii] + delta[ii,i] = 1;

# (R9) Un camión de entrada no puede precederse a sí mismo
subject to R9 {i in IN}:
    sum{ii in IN: i <> ii} delta[i,i] = 0;
    # Nota: delta[i,i] no está definido; esta restricción la maneja el conjunto {i <> ii}

# (R10) Secuenciación de camiones de salida
subject to R10 {j in OUT, jj in OUT: j <> jj}:
    E[jj] >= F[j] + sc - M * (1 - sigma[j,jj]);

# (R11) Orden total para camiones de salida
subject to R11 {j in OUT, jj in OUT: j < jj}:
    sigma[j,jj] + sigma[jj,j] = 1;

# (R12) Un camión de salida no puede precederse a sí mismo
# (garantizado por la definición del conjunto {j <> jj})

# (R13) Si hay transferencia de i a j, el camión j no puede salir antes de que i termine + manejo
subject to R13 {i in IN, j in OUT}:
    F[j] >= D[i] + h - M * (1 - z[i,j]);
