# Aggregation Tool for Distribution Grids

A Python implementation of the **Opti-KRON** algorithm

## Method & Reference

The underlying methodology is based on:

> Mokhtari, O., Chevalier, S., and Almassalkhi, M. 2026. *Optimal Kron-based Reduction of Networks (Opti-KRON) for Three-phase Distribution Feeders.* arXiv:2510.19608v2. [https://arxiv.org/abs/2510.19608](https://arxiv.org/abs/2510.19608)

## Repository Structure

* `opti_kron_single_phase.py` - Single-phase distribution network reduction module.
* `opti_kron_three_phase.py` - Extension for unbalanced three-phase power flow models.
* `plot_results.py` - Visualization utilities for plotting network before and after aggregation 

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
