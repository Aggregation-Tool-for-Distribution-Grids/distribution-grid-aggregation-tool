import numpy as np
import pandapower as pp
import pandapower.networks as pn
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from pandapower.pd2ppc import _pd2ppc_recycle
from pandapower.pypower.makeYbus import makeYbus
from pandapower.plotting import simple_plot

import time
from plot_results import plot_results

def opti_kron_three_phase(network, threshold, show_network=True):

    ## run power flow for 3 phase:
    # - take network parameters
    # - makes calculation
    # - update network
    pp.runpp_3ph(network)

    # show network
    if show_network:
        simple_plot(network)

    bus_lookup = network._pd2ppc_lookups['bus']
    num_nodes = len(network.bus)

    ### calculate Ybus for every phase
    _, ppci0 = _pd2ppc_recycle(network, 0, recycle=None)
    _, ppci1 = _pd2ppc_recycle(network, 1, recycle=None)
    _, ppci2 = _pd2ppc_recycle(network, 2, recycle=None)

    y_0_bus, _, _ = makeYbus(ppci0["baseMVA"], ppci0["bus"], ppci0["branch"])
    y_1_bus, _, _ = makeYbus(ppci1["baseMVA"], ppci1["bus"], ppci1["branch"])
    y_2_bus, _, _ = makeYbus(ppci2["baseMVA"], ppci2["bus"], ppci2["branch"])

    ### Transform: seq (0,1,2) -> phases (A,B,C)

    # complex rotation operator
    a = np.exp(1j * 2 * np.pi / 3)

    # matrix of Fortescue symmetrical components
    # of size 3 * 3
    A_fort = np.array([[1, 1, 1],
                       [1, a ** 2, a],
                       [1, a, a ** 2]])

    # Fortescue inverse matrix for the reverse transition
    A_inv = np.linalg.inv(A_fort)

    # creating projectors [E]
    E0 = sp.csc_matrix([[1, 0, 0], [0, 0, 0], [0, 0, 0]])
    E1 = sp.csc_matrix([[0, 0, 0], [0, 1, 0], [0, 0, 0]])
    E2 = sp.csc_matrix([[0, 0, 0], [0, 0, 0], [0, 0, 1]])

    Y_seq = (sp.kron(y_0_bus, E0)
             + sp.kron(y_1_bus, E1)
             + sp.kron(y_2_bus, E2)).tocsc()

    T = sp.kron(
        sp.identity(num_nodes),
        sp.csc_matrix(A_fort),
        format='csc'
    )

    Tinv = sp.kron(
        sp.identity(num_nodes),
        sp.csc_matrix(A_inv),
        format='csc')

    # Fortescue transformation for node admittance matrices
    Y_3ph = (T @ Y_seq @ Tinv).tocsc()

    ### voltages from pandapower 3ph results
    V_abc = np.zeros((num_nodes, 3), dtype=complex)
    node_phases = {}

    for pd_idx in network.bus.index:
        ppc_idx = bus_lookup[pd_idx]
        r = network.res_bus_3ph.loc[pd_idx]

        vm_a = r.vm_a_pu if np.isfinite(r.vm_a_pu) else 0.0
        vm_b = r.vm_b_pu if np.isfinite(r.vm_b_pu) else 0.0
        vm_c = r.vm_c_pu if np.isfinite(r.vm_c_pu) else 0.0

        va_a = r.va_a_degree if np.isfinite(r.va_a_degree) else 0.0
        va_b = r.va_b_degree if np.isfinite(r.va_b_degree) else 0.0
        va_c = r.va_c_degree if np.isfinite(r.va_c_degree) else 0.0

        V_abc[ppc_idx] = [
            vm_a * np.exp(1j * np.deg2rad(va_a)),
            vm_b * np.exp(1j * np.deg2rad(va_b)),
            vm_c * np.exp(1j * np.deg2rad(va_c)),
        ]

        phases = set()
        if vm_a > 0.05: phases.add(0)
        if vm_b > 0.05: phases.add(1)
        if vm_c > 0.05: phases.add(2)
        node_phases[ppc_idx] = phases

    V_hat_3ph = V_abc.reshape(-1)
    I_current_3ph = Y_3ph @ V_hat_3ph
    I_current_copy = I_current_3ph.copy()

    # get slack buses
    unique_buses = {int(bus_lookup[b]) for b in network.ext_grid['bus'].values}
    slack_buses = np.array(sorted(unique_buses))

    # calculates the row and column indices
    # of a three-phase conductance matrix
    slack_ph = np.array(
        [3 * s + k for s in slack_buses for k in range(3)],
        dtype=int
    )

    # determining free nodes
    free = np.setdiff1d(np.arange(3 * num_nodes), slack_ph)

    # block partitioning of Ybus
    Y_ff = Y_3ph[free][:, free].tocsc()
    Y_fs = Y_3ph[free][:, slack_ph].tocsc()
    V_s = V_hat_3ph[slack_ph]

    # LU decomposition
    lu = spla.splu(Y_ff)

    ## finding all possible candidates

    # determine which nodes are directly connected (Adjacency Matrix)
    # - if 1, node i and j connected
    # - else 0

    A_adj = np.zeros((num_nodes, num_nodes), dtype=int)

    for _, line in network.line.iterrows():
        if line.get('in_service', True):
            u = bus_lookup[line['from_bus']]
            v = bus_lookup[line['to_bus']]

            A_adj[u, v] = 1
            A_adj[v, u] = 1

    for _, trafo in network.trafo.iterrows():
        if trafo.get('in_service', True):
            u = bus_lookup[trafo['hv_bus']]
            v = bus_lookup[trafo['lv_bus']]

            A_adj[u, v] = 1
            A_adj[v, u] = 1

    orig_adj = A_adj.copy()
    orig_degrees = np.sum(orig_adj, axis=1)

    step = 0
    adj = A_adj.copy()
    cur_i = I_current_3ph.copy()
    rep = np.arange(num_nodes)

    while True:
        # we are looking for all pairs of candidates
        rows, cols = np.where(adj == 1)

        # we can't use slack bus
        candidates = []
        for s, r in zip(rows, cols):
            if r not in slack_buses and s not in slack_buses:
                if node_phases[r].issubset(node_phases[s]):
                    candidates.append((s, r))

        # the graph is completely reduced
        # or there are no connections left.
        if not candidates:
            break

        ### find the best pair
        best_pair = None
        best_smice = np.inf

        for s, r in candidates:
            i_c = cur_i.copy()
            i_c[3 * s:3 * s + 3] += i_c[3 * r:3 * r + 3]
            i_c[3 * r:3 * r + 3] = 0.0

            try:
                v_c = V_hat_3ph.copy()
                v_c[free] = lu.solve(i_c[free] - Y_fs @ V_s)
                v_c = v_c.reshape(-1, 3)
            except Exception:
                continue

            rep_c = rep.copy()
            rep_c[rep_c == r] = s
            v_c = v_c[rep_c]

            mice = np.abs(np.abs(V_abc) - np.abs(v_c))

            if mice.max() > threshold:
                continue

            node_max = mice.max(axis=1)
            cl_max = np.zeros(num_nodes)
            np.maximum.at(cl_max, rep_c, node_max)
            smice = cl_max[np.unique(rep_c)].sum()

            if smice < best_smice:
                best_smice, best_pair = smice, (s, r)

        if best_pair is None:
            break

        s_best, r_best = best_pair
        step += 1

        print(f"Step {step}: Collapse node {r_best} "
              f"-> node {s_best} | SMICE = {best_smice}")

        cur_i[3 * s_best:3 * s_best + 3] += cur_i[3 * r_best:3 * r_best + 3]
        cur_i[3 * r_best:3 * r_best + 3] = 0.0

        for nb in np.where(adj[r_best] == 1)[0]:
            if nb != s_best:
                adj[s_best, nb] = adj[nb, s_best] = 1

        adj[r_best, :] = 0
        adj[:, r_best] = 0

        rep[rep == r_best] = s_best

    # update cluster
    clusters = {i: [] for i in range(num_nodes)}
    for node in range(num_nodes):
        clusters[rep[node]].append(node)

    ## network radialization

    while True:
        triangles = []
        for i in range(adj.shape[0]):
            for j in range(i + 1, adj.shape[0]):
                if adj[i, j]:
                    for k in range(j + 1, adj.shape[0]):
                        if adj[i, k] and adj[j, k]:
                            triangles.append((i, j, k))

        if not triangles:
            break

        u, v, w = triangles[0]
        restored = False

        for sn in [u, v, w]:
            for member in list(clusters[sn]):
                if orig_degrees[member] >= 3:
                    if member in clusters[sn]:
                        clusters[sn].remove(member)
                    clusters[member] = [member]
                    rep[member] = member

                    for neighbor in np.where(orig_adj[member] == 1)[0]:
                        target_sn = None

                        for s_cand, members in clusters.items():
                            if neighbor in members:
                                target_sn = s_cand
                                break

                        if target_sn is not None and target_sn != member:
                            adj[member, target_sn] = 1
                            adj[target_sn, member] = 1

                    restored = True
                    break
            if restored:
                break

        if not restored:
            adj[u, v] = 0
            adj[v, u] = 0

    ### final error calculation

    I_final = np.zeros(3 * num_nodes, dtype=complex)
    for sn, members in clusters.items():
        if len(members) > 0:
            I_final[3 * sn: 3 * sn + 3] = sum(
                I_current_copy[3 * m: 3 * m + 3] for m in members
            )

    V_final = V_hat_3ph.copy()
    V_final[free] = lu.solve(I_final[free] - Y_fs @ V_s)
    V_final_rs = V_final.reshape(-1, 3)

    for sn, members in clusters.items():
        if len(members) > 0:
            for m in members:
                V_final_rs[m] = V_final_rs[sn]

    final_mice_error = np.abs(np.abs(V_abc) - np.abs(V_final_rs))
    maximum_mice_error = np.max(final_mice_error)
    mean_mice_error = np.mean(final_mice_error)

    final_smice_error = sum(
        np.max(final_mice_error[members])
        for sn, members in clusters.items()
        if len(members) > 0
    )

    active_supernodes = sorted(sn for sn, m in clusters.items() if len(m) > 0)
    reduced_nodes = [n for n in range(num_nodes) if n not in active_supernodes]

    ### kron reduction
    active_3ph = [3 * s + k for s in active_supernodes for k in range(3)]
    reduced_3ph = [3 * r + k for r in reduced_nodes for k in range(3)]

    Y_KK = Y_3ph[np.ix_(active_3ph, active_3ph)].tocsc()

    if len(reduced_3ph) == 0:
        Y_kron = Y_KK

    else:
        Y_KR = Y_3ph[np.ix_(active_3ph, reduced_3ph)].tocsc()
        Y_RK = Y_3ph[np.ix_(reduced_3ph, active_3ph)].tocsc()
        Y_RR = Y_3ph[np.ix_(reduced_3ph, reduced_3ph)].tocsc()

        inv_YRR_YRK = spla.spsolve(Y_RR, Y_RK)
        Y_kron = Y_KK - Y_KR.dot(inv_YRR_YRK)
        Y_kron = sp.csc_matrix(Y_kron)

    ### show result
    print(
        f"{'\nMax ΔV Error:'} {maximum_mice_error:.5f}\n"
        f"{'Mean ΔV Error:'} {mean_mice_error:.5f}\n"
        f"{'SMICE Error:'} {final_smice_error:.5f}\n"
        f"{'Active Supernodes:'} {len(active_supernodes)} / {num_nodes} "
        f"({100 * (num_nodes - len(active_supernodes)) / num_nodes:.2f}% reduced)\n"
    )

    return orig_adj, adj, num_nodes, slack_buses, active_supernodes, clusters, Y_kron


# start time
start = time.time()

## create 3-phase network
network = pn.ieee_european_lv_asymmetric()

# run function
orig_adj, A_adj, num_nodes, slack_buses, active_supernodes, clusters, Y_kron = opti_kron_three_phase(
    network,
    threshold=0.005,
    show_network=False
)

# end time
end = time.time()
print(f"\nTotal runtime of the program is {end - start:.3f} seconds")

# plot result
plot_results(orig_adj, A_adj, num_nodes, slack_buses, active_supernodes, clusters)
