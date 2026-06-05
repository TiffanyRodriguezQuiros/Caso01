import streamlit as st
import pulp
from itertools import combinations

st.set_page_config(page_title="Cross Docking - LogiFast CR", layout="wide")

st.title("Optimización de Cross Docking")
st.markdown("**Universidad de Costa Rica — Ingeniería Industrial**")
st.markdown("Modelo MIP basado en Yu & Egbelu (2008)")
st.divider()

# ── Parsear archivo ────────────────────────────────────────────────────────────
def parse_ts(contenido):
    tokens = contenido.split()
    idx = 0
    assert tokens[idx] == 'i'; idx += 1; I = int(tokens[idx]); idx += 1
    assert tokens[idx] == 'o'; idx += 1; J = int(tokens[idx]); idx += 1
    assert tokens[idx] == 'n'; idx += 1; K = int(tokens[idx]); idx += 1

    a = {(i, k): 0 for i in range(1, I+1) for k in range(1, K+1)}
    b = {(j, k): 0 for j in range(1, J+1) for k in range(1, K+1)}

    while idx < len(tokens):
        t = tokens[idx]; idx += 1
        num = int(tokens[idx]); idx += 1
        prod = int(tokens[idx]); idx += 1
        qty  = int(tokens[idx]); idx += 1
        if t == 'r':
            a[(num, prod)] = qty
        elif t == 's':
            b[(num, prod)] = qty

    return I, J, K, a, b


# ── Modelo MIP ────────────────────────────────────────────────────────────────
def resolver(I, J, K, a, b):
    IN   = range(1, I+1)
    OUT  = range(1, J+1)
    PROD = range(1, K+1)

    p, h, sc, M = 1, 5, 10, 10000

    pi_i = {i: p * sum(a[i, k] for k in PROD) for i in IN}
    pi_j = {j: p * sum(b[j, k] for k in PROD) for j in OUT}

    mdl = pulp.LpProblem("CrossDocking", pulp.LpMinimize)

    x = {(i,j,k): pulp.LpVariable(f"x_{i}_{j}_{k}", lowBound=0) for i in IN for j in OUT for k in PROD}
    y = {(i,j,k): pulp.LpVariable(f"y_{i}_{j}_{k}", lowBound=0) for i in IN for j in OUT for k in PROD}
    z     = {(i,j):  pulp.LpVariable(f"z_{i}_{j}",  cat='Binary') for i in IN for j in OUT}
    delta = {(i,ii): pulp.LpVariable(f"d_{i}_{ii}", cat='Binary') for i in IN for ii in IN if i != ii}
    sigma = {(j,jj): pulp.LpVariable(f"s_{j}_{jj}", cat='Binary') for j in OUT for jj in OUT if j != jj}
    A = {i: pulp.LpVariable(f"A_{i}", lowBound=0) for i in IN}
    D = {i: pulp.LpVariable(f"D_{i}", lowBound=0) for i in IN}
    E = {j: pulp.LpVariable(f"E_{j}", lowBound=0) for j in OUT}
    F = {j: pulp.LpVariable(f"F_{j}", lowBound=0) for j in OUT}
    C = pulp.LpVariable("C", lowBound=0)

    mdl += C

    for j in OUT: mdl += C >= F[j]

    for i in IN:
        for k in PROD:
            if a[i, k] > 0:
                mdl += pulp.lpSum(x[i,j,k] + y[i,j,k] for j in OUT) == a[i,k]
            else:
                for j in OUT:
                    mdl += x[i,j,k] == 0
                    mdl += y[i,j,k] == 0

    for j in OUT:
        for k in PROD:
            if b[j, k] > 0:
                mdl += pulp.lpSum(x[i,j,k] + y[i,j,k] for i in IN) == b[j,k]

    for i in IN:
        for j in OUT:
            mdl += pulp.lpSum(x[i,j,k] + y[i,j,k] for k in PROD) <= M * z[i,j]

    for i in IN: mdl += D[i] == A[i] + pi_i[i]
    for j in OUT: mdl += F[j] == E[j] + pi_j[j]

    for i in IN:
        for ii in IN:
            if i != ii:
                mdl += A[ii] >= D[i] + sc - M*(1 - delta[i,ii])

    for i, ii in combinations(IN, 2):
        mdl += delta[i,ii] + delta[ii,i] == 1

    for j in OUT:
        for jj in OUT:
            if j != jj:
                mdl += E[jj] >= F[j] + sc - M*(1 - sigma[j,jj])

    for j, jj in combinations(OUT, 2):
        mdl += sigma[j,jj] + sigma[jj,j] == 1

    for i in IN:
        for j in OUT:
            mdl += F[j] >= D[i] + h - M*(1 - z[i,j])

    mdl.solve(pulp.PULP_CBC_CMD(msg=0, timeLimit=300))

    if pulp.LpStatus[mdl.status] != 'Optimal':
        return None

    return {
        "makespan": pulp.value(C),
        "A": {i: pulp.value(A[i]) for i in IN},
        "D": {i: pulp.value(D[i]) for i in IN},
        "E": {j: pulp.value(E[j]) for j in OUT},
        "F": {j: pulp.value(F[j]) for j in OUT},
        "z": {(i,j): pulp.value(z[i,j]) for i in IN for j in OUT},
        "x": {(i,j,k): pulp.value(x[i,j,k]) for i in IN for j in OUT for k in PROD},
        "y": {(i,j,k): pulp.value(y[i,j,k]) for i in IN for j in OUT for k in PROD},
        "pi_i": pi_i,
        "pi_j": pi_j,
        "I": I, "J": J, "K": K,
    }


# ── Interfaz ──────────────────────────────────────────────────────────────────
archivo = st.file_uploader("Subí el archivo de datos (.txt en formato TS)", type="txt")

if archivo:
    contenido = archivo.read().decode("utf-8")

    try:
        I, J, K, a, b = parse_ts(contenido)
    except Exception:
        st.error("El archivo no tiene el formato correcto.")
        st.stop()

    st.success(f"Archivo cargado: {I} camiones de entrada, {J} camiones de salida, {K} tipos de producto")

    with st.spinner("Resolviendo el modelo MIP... esto puede tardar unos segundos"):
        res = resolver(I, J, K, a, b)

    if res is None:
        st.error("No se encontró solución óptima.")
        st.stop()

    IN   = range(1, I+1)
    OUT  = range(1, J+1)
    PROD = range(1, K+1)

    # ── Makespan ──────────────────────────────────────────────────────────────
    st.divider()
    col1, col2, col3 = st.columns(3)
    col1.metric("Makespan óptimo", f"{res['makespan']:.0f} min")
    horas = int(res['makespan'] // 60)
    minutos = int(res['makespan'] % 60)
    col2.metric("En horas", f"{horas} h {minutos} min")
    col3.metric("Estado", "Óptimo ✅")

    # ── Secuencia entrada ─────────────────────────────────────────────────────
    st.divider()
    st.subheader("Secuencia de camiones de entrada (muelle de recepción)")

    entrada_ordenada = sorted(IN, key=lambda i: res["A"][i])
    filas_entrada = []
    for pos, i in enumerate(entrada_ordenada, 1):
        filas_entrada.append({
            "Posición": f"{pos}°",
            "Camión": f"Entrada {i}",
            "Inicio (min)": int(res["A"][i]),
            "Fin (min)": int(res["D"][i]),
            "Duración (min)": int(res["pi_i"][i])
        })

    import pandas as pd
    st.dataframe(pd.DataFrame(filas_entrada), use_container_width=True, hide_index=True)

    secuencia_entrada = " → ".join([f"Camión {i}" for i in entrada_ordenada])
    st.markdown(f"**Orden:** {secuencia_entrada}")

    # ── Secuencia salida ──────────────────────────────────────────────────────
    st.divider()
    st.subheader("Secuencia de camiones de salida (muelle de despacho)")

    salida_ordenada = sorted(OUT, key=lambda j: res["E"][j])
    filas_salida = []
    for pos, j in enumerate(salida_ordenada, 1):
        filas_salida.append({
            "Posición": f"{pos}°",
            "Camión": f"Salida {j}",
            "Inicio (min)": int(res["E"][j]),
            "Fin (min)": int(res["F"][j]),
            "Duración (min)": int(res["pi_j"][j])
        })

    st.dataframe(pd.DataFrame(filas_salida), use_container_width=True, hide_index=True)

    secuencia_salida = " → ".join([f"Camión {j}" for j in salida_ordenada])
    st.markdown(f"**Orden:** {secuencia_salida}")

    # ── Flujos ────────────────────────────────────────────────────────────────
    st.divider()
    st.subheader("Flujos de transferencia")

    filas_flujo = []
    for i in IN:
        for j in OUT:
            if res["z"][(i,j)] and res["z"][(i,j)] > 0.5:
                directo = sum(res["x"][(i,j,k)] or 0 for k in PROD)
                temporal = sum(res["y"][(i,j,k)] or 0 for k in PROD)
                if directo > 0.01 or temporal > 0.01:
                    filas_flujo.append({
                        "Camión entrada": f"Entrada {i}",
                        "Camión salida": f"Salida {j}",
                        "Directo (u)": int(round(directo)),
                        "Vía almac. temporal (u)": int(round(temporal)),
                        "Total (u)": int(round(directo + temporal))
                    })

    if filas_flujo:
        st.dataframe(pd.DataFrame(filas_flujo), use_container_width=True, hide_index=True)
    else:
        st.info("No hay flujos de transferencia registrados.")
