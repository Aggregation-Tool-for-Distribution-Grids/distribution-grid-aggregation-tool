import numpy as np
import pandapower as pp
import pandapower.networks as pn

import scipy.sparse as sp
import scipy.sparse.linalg as spla

from pandapower.plotting import simple_plot

import time
from plot_results import plot_results

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Pt

def opti_kron_single_phase(network, threshold, show_network=True):

    ## run power flow:
    # - take network parameters
    # - makes calculation
    # - update network
    pp.runpp(network)

    # show network
    if show_network:
        simple_plot(network)

    ## take matrices and vectors, for opti-kron calculation

    # conductivity matrix
    Y_matrix = network._ppc['internal']['Ybus'].tocsc()

    num_nodes = Y_matrix.shape[0]
    bus_lookup = network._pd2ppc_lookups['bus']

    # voltage module [p.u] and angle [radians]
    v_mag = np.zeros(num_nodes)
    v_ang = np.zeros(num_nodes)

    for pd_idx, ppc_idx in enumerate(bus_lookup):
        v_mag[ppc_idx] = network.res_bus.vm_pu.loc[pd_idx]
        v_ang[ppc_idx] = np.deg2rad(network.res_bus.va_degree.loc[pd_idx])

    # complex vector voltage
    V_hat = v_mag * np.exp(1j * v_ang)

    # initial current vector [I = Y * V]
    I_current = Y_matrix.dot(V_hat)
    I_current_copy = I_current.copy()

    # get slack buses
    unique_buses = set()
    for bus in network.ext_grid['bus'].values:
        mapped_id = int(bus_lookup[bus])
        unique_buses.add(mapped_id)

    slack_buses = np.array(sorted(unique_buses))

    # determining free nodes
    free = np.setdiff1d(np.arange(num_nodes), slack_buses)

    # block partitioning of Ybus
    Y_ff = Y_matrix[free][:, free].tocsc()
    Y_fs = Y_matrix[free][:, slack_buses].tocsc()
    V_s = V_hat[slack_buses]

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
    clusters = {i: [i] for i in range(num_nodes)}
    cur_adj = A_adj.copy()

    while True:
        # we are looking for all pairs of candidates
        rows, cols = np.where(cur_adj == 1)

        # we can't use slack bus
        candidates = []
        for s, r in zip(rows, cols):
            if r not in slack_buses and len(clusters[s]) > 0 and len(clusters[r]) > 0:
                candidates.append((s, r))

        # the graph is completely reduced
        # or there are no connections left.
        if not candidates:
            break

        ### find the best pair
        best_pair = None
        best_smice = np.inf

        for s, r in candidates:
            # copy the current so as not to change the original results
            I_cand = I_current.copy()

            I_cand[s] += I_cand[r]
            I_cand[r] = 0.0

            # solves a system of linear equations
            try:
                # V_cand = spla.spsolve(Y_matrix, I_cand)
                V_cand = np.zeros(num_nodes, dtype=complex)
                V_cand[slack_buses] = V_s
                V_cand[free] = lu.solve(I_cand[free] - Y_fs @ V_s)
            except Exception:
                continue

            temp_clusters = {}
            for node, members in clusters.items():
                temp_clusters[node] = list(members)

            temp_clusters[s].extend(temp_clusters[r])
            temp_clusters[r] = []

            for sn, members in temp_clusters.items():
                if len(members) > 0:
                    for member in members:
                        V_cand[member] = V_cand[sn]

            # MICE
            mice = np.abs(np.abs(V_hat) - np.abs(V_cand))

            if np.max(mice) > threshold:
                continue

            # SMICE
            smice = 0.0
            for sn, members in temp_clusters.items():
                if len(members) > 0:
                    cluster_mice = np.max(mice[members])
                    smice += cluster_mice

            if smice < best_smice:
                best_smice = smice
                best_pair = (s, r)

        # no candidates with SMICE error <= E
        if best_pair is None or best_smice == np.inf:
            break

        ## merging pairs and updating the graph
        s_best, r_best = best_pair
        step += 1

        print(f"Step {step}: Collapse node {r_best} "
              f"-> node {s_best} | SMICE = {best_smice:.6f}")

        # we fix the transfer of load current
        I_current[s_best] += I_current[r_best]
        I_current[r_best] = 0.0

        # updating the graph topology
        # all neighbors of node r become neighbors of node s
        neighbors_r = np.where(A_adj[r_best] == 1)[0]

        for neighbor in neighbors_r:
            # eliminate the loop s -> s
            if neighbor != s_best:
                A_adj[s_best, neighbor] = 1
                A_adj[neighbor, s_best] = 1

        # completely "disconnect" node r from the graph
        A_adj[r_best, :] = 0
        A_adj[:, r_best] = 0

        # update cluster tracking (to know which nodes entered s)
        clusters[s_best].extend(clusters[r_best])
        clusters[r_best] = []

    ## network radialization

    # finding closed triangles
    while True:
        triangles = []
        for i in range(cur_adj.shape[0]):
            for j in range(i + 1, cur_adj.shape[0]):
                if cur_adj[i, j]:
                    for k in range(j + 1, cur_adj.shape[0]):
                        if cur_adj[i, k] and cur_adj[j, k]:
                            triangles.append((i, j, k))

        # radial network
        if not triangles:
            break

        # we search for a critical node
        # among the reduced members of these supernodes,
        # which had a high degree of connectivity (degree >= 3)
        # in the original network
        u, v, w = triangles[0]
        restored = False

        for sn in [u, v, w]:
            for member in list(clusters[sn]):
                if orig_degrees[member] >= 3:
                    if member in clusters[sn]:
                        clusters[sn].remove(member)
                    clusters[member] = [member]

                    for neighbor in np.where(orig_adj[member] == 1)[0]:
                        target_sn = None

                        for s_candidate, members in clusters.items():
                            if neighbor in members:
                                target_sn = s_candidate
                                break

                        if target_sn is not None and target_sn != member:
                            A_adj[member, target_sn] = 1
                            A_adj[target_sn, member] = 1

                    restored = True
                    break

            if restored:
                break

        if not restored:
            A_adj[u, v] = 0
            A_adj[v, u] = 0

    ### Final error check after radialization
    I_final = np.zeros(num_nodes, dtype=complex)

    for sn, members in clusters.items():
        if len(members) > 0:
            I_final[sn] = sum(I_current_copy[m] for m in members)

    V_final = np.zeros(num_nodes, dtype=complex)
    V_final[slack_buses] = V_s
    V_final[free] = lu.solve(I_final[free] - Y_fs @ V_s)

    for sn, members in clusters.items():
        if len(members) > 0:
            for m in members:
                V_final[m] = V_final[sn]

    final_mice_error = np.abs(np.abs(V_hat) - np.abs(V_final))
    maximum_mice_error = np.max(final_mice_error)
    mean_mice_error = np.mean(final_mice_error)

    final_smice_error = sum(
        np.max(final_mice_error[members])
        for sn, members in clusters.items()
        if len(members) > 0
    )

    active_supernodes = []
    for node, members in clusters.items():
        if len(members) > 0:
            active_supernodes.append(node)

    reduced_nodes = []
    for node in range(num_nodes):
        if node not in active_supernodes:
            reduced_nodes.append(node)

    # calculate the Y_kron matrix
    Y_KK = Y_matrix[np.ix_(active_supernodes, active_supernodes)].tocsc()

    if len(reduced_nodes) == 0:
        Y_kron = Y_KK

    else:
        # Ybus matrix blocks for Kron reduction
        Y_KR = Y_matrix[np.ix_(active_supernodes, reduced_nodes)].tocsc()
        Y_RK = Y_matrix[np.ix_(reduced_nodes, active_supernodes)].tocsc()
        Y_RR = Y_matrix[np.ix_(reduced_nodes, reduced_nodes)].tocsc()

        # Y_RR * X = Y_RK
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

    metrics = {
        "num_nodes": num_nodes,
        "active_supernodes": len(active_supernodes),
        "max_error": maximum_mice_error,
        "mean_error": mean_mice_error,
        "smice_error": final_smice_error,
    }

    return orig_adj, A_adj, num_nodes, slack_buses, active_supernodes, clusters, Y_kron, metrics


# start time
start = time.time()

# create simple radial networks:
# - pn.create_cigre_network_mv(with_der="all")
# - pn.create_kerber_vorstadtnetz_kabel_2()

# network = pn.create_kerber_vorstadtnetz_kabel_2()

# run function
# orig_adj, A_adj, num_nodes, slack_buses, active_supernodes, clusters, Y_kron, metrics = opti_kron_single_phase(
#     network,
#     threshold=0.005,
#     show_network=False
# )

# end time
# end = time.time()
# print(f"\nTotal runtime of the program is {end - start:.3f} seconds")

# plot result
# plot_results(orig_adj, A_adj, num_nodes, slack_buses, active_supernodes, clusters)



networks = [
    pn.create_cigre_network_mv(with_der="all"),
    pn.create_kerber_vorstadtnetz_kabel_2()
]

thresholds = [
    0.0001,
    0.001,
    0.01,
]

rows = []
for id, network in enumerate(networks):
    for threshold in thresholds:
        result = opti_kron_single_phase(
            network,
            threshold=threshold,
            show_network=False
        )

        *_, metrics = result

        original = metrics["num_nodes"]
        reduced = metrics["active_supernodes"]
        reduction_pct = 100.0 * (original - reduced) / original

        rows.append({
            "network": str(id),
            "original": original,
            "reduced": reduced,
            "reduction": reduction_pct,
            "threshold": threshold,
            "max_dv": metrics["max_error"],
            "mean_dv": metrics["mean_error"],
        })


doc = Document()
doc.add_heading("Opti-Kron Reduction Results", level=1)

headers = [
    "Network", "Original Nodes", "Reduced Nodes",
    "Reduction [%]", "Threshold", "Max ΔV", "Mean ΔV",
]

table = doc.add_table(rows=1, cols=len(headers))

hdr_cells = table.rows[0].cells
for i, h in enumerate(headers):
    hdr_cells[i].text = h
    for p in hdr_cells[i].paragraphs:
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        for r in p.runs:
            r.bold = True
            r.font.size = Pt(10)

for row in rows:
    cells = table.add_row().cells
    cells[0].text = row["network"]
    cells[1].text = str(row["original"])
    cells[2].text = str(row["reduced"])
    cells[3].text = f"{row['reduction']:.2f}"
    cells[4].text = f"{row['threshold']:.4f}"
    cells[5].text = f"{row['max_dv']:.5f}"
    cells[6].text = f"{row['mean_dv']:.5f}"

doc.save("reduction_results.docx")







