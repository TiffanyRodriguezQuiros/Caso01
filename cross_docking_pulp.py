"""
LogiFast CR - Cross Docking MIP Scheduler
Basado en Yu & Egbelu (2008)
Universidad de Costa Rica - Ingeniería Industrial

Versión PuLP (solver CBC, sin licencia comercial)
Uso:  python cross_docking_pulp.py TS5.txt
"""
import sys
from itertools import combinations
import pulp

def parse_ts(path):
    """Lee archivo de datos en formato TS."""
    with open(path) as f:
        tokens = f.read().split()
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


def solve(data_path):
    I, J, K, a, b = parse_ts(data_path)
    IN   = range(1, I+1)
    OUT  = range(1, J+1)
    PROD = range(1, K+1)

    # Parámetros operativos
    p, h, sc, M = 1, 5, 10, 10000

    pi_i = {i: p * sum(a[i, k] for k in PROD) for i in IN}
    pi_j = {j: p * sum(b[j, k] for k in PROD) for j in OUT}

    mdl = pulp.LpProblem("CrossDocking_LogiFast", pulp.LpMinimize)

    # ── Variables ─────────────────────────────────────────────
    x = {(i,j,k): pulp.LpVariable(f"x_{i}_{j}_{k}", lowBound=0) for i in IN for j in OUT for k in PROD}
    y = {(i,j,k): pulp.LpVariable(f"y_{i}_{j}_{k}", lowBound=0) for i in IN for j in OUT for k in PROD}
    z     = {(i,j):   pulp.LpVariable(f"z_{i}_{j}",   cat='Binary') for i in IN for j in OUT}
    delta = {(i,ii):  pulp.LpVariable(f"d_{i}_{ii}",  cat='Binary') for i in IN for ii in IN if i != ii}
    sigma = {(j,jj):  pulp.LpVariable(f"s_{j}_{jj}",  cat='Binary') for j in OUT for jj in OUT if j != jj}
    A = {i: pulp.LpVariable(f"A_{i}", lowBound=0) for i in IN}
    D = {i: pulp.LpVariable(f"D_{i}", lowBound=0) for i in IN}
    E = {j: pulp.LpVariable(f"E_{j}", lowBound=0) for j in OUT}
    F = {j: pulp.LpVariable(f"F_{j}", lowBound=0) for j in OUT}
    C = pulp.LpVariable("C", lowBound=0)

    # ── Objetivo ──────────────────────────────────────────────
    mdl += C, "Minimizar_Makespan"

    # ── Restricciones ─────────────────────────────────────────
    # R1: makespan
    for j in OUT:
        mdl += C >= F[j], f"R1_j{j}"

    # R2: conservación oferta inbound
    for i in IN:
        for k in PROD:
            if a[i, k] > 0:
                mdl += pulp.lpSum(x[i,j,k] + y[i,j,k] for j in OUT) == a[i,k], f"R2_i{i}_k{k}"
            else:
                for j in OUT:
                    mdl += x[i,j,k] == 0, f"R2x0_i{i}_j{j}_k{k}"
                    mdl += y[i,j,k] == 0, f"R2y0_i{i}_j{j}_k{k}"

    # R3: satisfacción demanda outbound
    for j in OUT:
        for k in PROD:
            if b[j, k] > 0:
                mdl += pulp.lpSum(x[i,j,k] + y[i,j,k] for i in IN) == b[j,k], f"R3_j{j}_k{k}"

    # R4: relación z con flujos
    for i in IN:
        for j in OUT:
            mdl += pulp.lpSum(x[i,j,k] + y[i,j,k] for k in PROD) <= M * z[i,j], f"R4_i{i}_j{j}"

    # R5: tiempo proceso entrada
    for i in IN:
        mdl += D[i] == A[i] + pi_i[i], f"R5_i{i}"

    # R6: tiempo proceso salida
    for j in OUT:
        mdl += F[j] == E[j] + pi_j[j], f"R6_j{j}"

    # R7: secuenciación entrada
    for i in IN:
        for ii in IN:
            if i != ii:
                mdl += A[ii] >= D[i] + sc - M*(1 - delta[i,ii]), f"R7_i{i}_ii{ii}"

    # R8: orden total entrada
    for i, ii in combinations(IN, 2):
        mdl += delta[i,ii] + delta[ii,i] == 1, f"R8_i{i}_ii{ii}"

    # R9: no auto-precedencia entrada — garantizada por definición de delta

    # R10: secuenciación salida
    for j in OUT:
        for jj in OUT:
            if j != jj:
                mdl += E[jj] >= F[j] + sc - M*(1 - sigma[j,jj]), f"R10_j{j}_jj{jj}"

    # R11: orden total salida
    for j, jj in combinations(OUT, 2):
        mdl += sigma[j,jj] + sigma[jj,j] == 1, f"R11_j{j}_jj{jj}"

    # R12: no auto-precedencia salida — garantizada por definición de sigma

    # R13: tiempo de transferencia
    for i in IN:
        for j in OUT:
            mdl += F[j] >= D[i] + h - M*(1 - z[i,j]), f"R13_i{i}_j{j}"

    # ── Resolver ──────────────────────────────────────────────
    solver = pulp.PULP_CBC_CMD(msg=0, timeLimit=300)
    mdl.solve(solver)

    print(f"\n{'='*55}")
    print(f"  Estado: {pulp.LpStatus[mdl.status]}")
    print(f"  MAKESPAN OPTIMO: {pulp.value(C):.0f} minutos")
    print(f"{'='*55}")

    A_v  = {i: pulp.value(A[i]) for i in IN}
    D_v  = {i: pulp.value(D[i]) for i in IN}
    E_v  = {j: pulp.value(E[j]) for j in OUT}
    F_v  = {j: pulp.value(F[j]) for j in OUT}
    z_v  = {(i,j): pulp.value(z[i,j]) for i in IN for j in OUT}

    print("\nSecuencia camiones de ENTRADA (orden de atencion en muelle):")
    pos = 1
    for i in sorted(IN, key=lambda i: A_v[i]):
        print(f"  Posicion {pos}: Camion entrada {i}  |  t={A_v[i]:.0f} min -> {D_v[i]:.0f} min  (duracion {pi_i[i]} min)")
        pos += 1

    print("\nSecuencia camiones de SALIDA (orden de atencion en muelle):")
    pos = 1
    for j in sorted(OUT, key=lambda j: E_v[j]):
        print(f"  Posicion {pos}: Camion salida  {j}  |  t={E_v[j]:.0f} min -> {F_v[j]:.0f} min  (duracion {pi_j[j]} min)")
        pos += 1

    print("\nFlujos de transferencia:")
    for i in IN:
        for j in OUT:
            if z_v[(i,j)] and z_v[(i,j)] > 0.5:
                direct = sum(pulp.value(x[i,j,k]) or 0 for k in PROD)
                via_ts = sum(pulp.value(y[i,j,k]) or 0 for k in PROD)
                tag = []
                if direct > 0.01: tag.append(f"directo={direct:.0f}")
                if via_ts > 0.01: tag.append(f"almac.temporal={via_ts:.0f}")
                print(f"  Camion entrada {i} -> Camion salida {j}: {', '.join(tag)} unidades")


if __name__ == '__main__':
    path = sys.argv[1] if len(sys.argv) > 1 else 'TS5.txt'
    solve(path)
