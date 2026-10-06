# Aggregation Tool for Distribution Grids

A Python implementation of the **Opti-KRON** algorithm

## Method & Reference

The underlying methodology is based on:

> Mokhtari, O., Chevalier, S., and Almassalkhi, M. 2026. *Optimal Kron-based Reduction of Networks (Opti-KRON) for Three-phase Distribution Feeders.* arXiv:2510.19608v2. [https://arxiv.org/abs/2510.19608](https://arxiv.org/abs/2510.19608)

## Repository Structure

* `opti_kron_single_phase.py` - Single-phase distribution network reduction module.
* `opti_kron_three_phase.py` - Extension for unbalanced three-phase power flow models.
* `plot_results.py` - Visualization utilities for plotting network before and after aggregation.
* `create_network_example.py` - A code example of how to create a network from scratch in pandapower.
* `export_and_import_network.py` - Code example for exporting and importing files in various formats.

## How to run a script

1. **Environment Setup**: Clone the repository, create a Python environment, and install the required dependencies: `pip install -r requirements.txt`
2. **Select the Module**: Choose the appropriate reduction script for your analysis:
    * `opti_kron_single_phase.py`: For balanced single-phase networks.
    * `opti_kron_three_phase.py`: For unbalanced three-phase networks.
3. **Prepare the Network Model**: Create or load a pandapower network object (pandapowerNet):
    * **Custom networks**: Build your own model using pandapower elements (refer to `create_network_example.py` for a starter template or check the [pandapower element documentation](https://pandapower.readthedocs.io/en/latest/about.html)).
    * **Benchmark networks**: Load standard test feeders directly from `pandapower.networks` (e.g., `pn.create_cigre_network_mv(with_der="all")` or `pn.create_kerber_vorstadtnetz_kabel_2()`).
4. **Configure Function Parameters**: Call the core function (`opti_kron_single_phase()` or `opti_kron_three_phase()`) passing the required arguments:
    * `network`: The input pandapowerNet instance.
    * `threshold`: Maximum allowable voltage magnitude deviation in p.u.
    * `show_network`: Boolean flag (`True`/`False`) to visually render the grid before reduction.
5. **Execute and Review Results**: Run the script to perform the network reduction:

```
orig_adj, A_adj, num_nodes, slack_buses, active_supernodes, clusters, Y_kron, metrics = opti_kron_single_phase(
    network,
    threshold=0.005,
    show_network=False
)
```

* **Terminal Outputs**: Real-time logging displays step-by-step node collapses and final performance indicators (Max delta V, Mean delta V, SMICE error, and reduction percentage).
* **Visualization & Export**: Optional helper functions let you visualize the reduced topology using `plot_results(...)` or export summary benchmark tables directly to a Word document.

## Benchmarks & Performance Results

### Below are the experimental results for single-phase networks.

| Network ID | Original Nodes | Reduced Nodes | Reduction [%] | Threshold | Max delta V [p.u.] | Mean delta V [p.u.] |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **CIGRE MV Grid** | 18 | 18 | 0.00% | 0.0001 | 0.00000 | 0.00000 |
| **CIGRE MV Grid** | 18 | 16 | 11.11% | 0.0010 | 0.00078 | 0.00010 |
| **CIGRE MV Grid** | 18 | 12 | 33.33% | 0.0100 | 0.00771 | 0.00186 |
| **Kerber Networks** | 290 | 146 | 49.66% | 0.0001 | 0.00009 | 0.00004 |
| **Kerber Networks** | 290 | 144 | 50.34% | 0.0010 | 0.00059 | 0.00005 |
| **Kerber Networks** | 290 | 143 | 50.69% | 0.0100 | 0.00520 | 0.00006 |
