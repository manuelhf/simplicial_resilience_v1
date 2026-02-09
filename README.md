# WDSA-TDA: Topological Data Analysis for Water Distribution System Resilience

A Python package for analyzing water distribution network resilience using topological data analysis (TDA) methods, specifically simplicial complexes.

## Overview

This package implements methods to:

1. **Build Simplicial Complexes** from water network topology and hydraulic data
   - Structural complex (from graph topology)
   - Functional complex (from hydraulic head correlations)

2. **Compute Topological Invariants**
   - Betti numbers (β₀, β₁, β₂)
   - Functional Cohesion Index (FCI)

3. **Track Operational Changes** during disruption scenarios
   - Time-varying operational complex
   - Correlation with hydraulic performance

4. **Generate Publication-Quality Outputs**
   - Figures (PDF, 300 DPI)
   - LaTeX tables

## Installation

### Requirements

```
wntr >= 0.5.0
networkx >= 2.6
numpy >= 1.20
pandas >= 1.3
matplotlib >= 3.4
scipy >= 1.7
```

### Setup

```bash
# Clone or download the package
cd wdsa_tda_package

# Install dependencies
pip install wntr networkx numpy pandas matplotlib scipy
```

## Quick Start

### Command Line

```bash
# Basic usage with default settings
python main.py

# Specify input file and output directory
python main.py --input MyNetwork.inp --output results/

# Adjust pressure threshold
python main.py --pressure 25.0
```

### Python API

```python
from wdsa_tda_package import run_analysis

# Run complete analysis
results = run_analysis('network.inp', output_dir='outputs')

# Access results
print(f"Redundancy gap: {results['redundancy_gap']:.1f}%")
print(f"Structural triangles: {results['structural_metrics']['n_triangles']}")
```

### Using Individual Modules

```python
import wntr
import networkx as nx

from wdsa_tda_package.simplicial import build_structural_complex, build_functional_complex
from wdsa_tda_package.topology import compute_complex_metrics, compute_betti_numbers

# Load network
wn = wntr.network.WaterNetworkModel('network.inp')
G = nx.Graph(wntr.network.to_graph(wn))
pos = {name: obj.coordinates for name, obj in wn.nodes()}

# Build structural complex
triangles, cycles, n = build_structural_complex(G, pos)
print(f"Found {len(triangles)} structural triangles")

# Compute metrics
metrics = compute_complex_metrics(list(G.nodes()), edges, triangles)
print(f"β₀={metrics['beta0']}, β₁={metrics['beta1']}, FCI={metrics['fci']:.2f}")
```

## Package Structure

```
wdsa_tda_package/
├── __init__.py          # Package initialization and exports
├── simplicial.py        # Simplicial complex construction
├── topology.py          # Betti numbers and FCI computation
├── operational.py       # Time-varying operational complex tracking
├── scenarios.py         # Disruption scenario definition and execution
├── visualization.py     # Publication-quality figure generation
├── tables.py            # LaTeX table generation
├── main.py              # Main entry point script
└── README.md            # This file
```

## Modules

### simplicial.py

Functions for building simplicial complexes:

- `build_structural_complex(G, pos)` - Mesh-loop triangulation from graph
- `build_functional_complex(results, junctions)` - Vietoris-Rips from hydraulic data
- `extract_edges_from_triangles(triangles)` - Extract edge set

### topology.py

Topological invariant computation:

- `compute_betti_numbers(V, E, T)` - Compute β₀, β₁, β₂
- `compute_fci(triangles, overlap, beta0)` - Functional Cohesion Index
- `compute_complex_metrics(V, E, T)` - All metrics in one call

### operational.py

Time-varying analysis:

- `track_operational_complex(wn, results, G, pos)` - Track complex over time
- `compute_hydraulic_metrics(wn, results)` - Service coverage, PDI

### scenarios.py

Scenario simulation:

- `identify_critical_components(wn, results)` - Find trunk/secondary mains
- `define_standard_scenarios(components)` - Create S1, S2 scenario sets
- `run_scenario(inp, name, config, G, pos)` - Execute single scenario
- `run_all_scenarios(inp, scenarios, G, pos)` - Execute all scenarios

### visualization.py

Figure generation:

- `plot_network_comparison(G, pos, struct, func)` - 3-panel baseline figure
- `plot_ensemble_separate_rows(results, prefix, baseline)` - 3×3 ensemble
- `plot_scenario_detail(data, name, baseline)` - 2×2 detail figure

### tables.py

LaTeX table generation:

- `generate_baseline_table(struct, func, eps)` - Baseline comparison
- `generate_scenario_table(results, prefix, baseline)` - Scenario results

## Output Files

The analysis generates:

### Figures (PDF, 300 DPI)

| File | Description |
|------|-------------|
| `Figure1_baseline.pdf` | Network graph + structural + functional complex |
| `Figure2_S1_ensemble.pdf` | S1a/b/c scenarios in separate rows |
| `Figure3_S2_ensemble.pdf` | S2a/b/c scenarios in separate rows |
| `Figure_S1a_detail.pdf` | Individual scenario 2×2 detail |

### Tables (LaTeX)

| File | Description |
|------|-------------|
| `Table1_baseline.tex` | Structural vs functional metrics |
| `Table2_S1_results.tex` | S1 ensemble results |
| `Table3_S2_results.tex` | S2 ensemble results |

### Data

| File | Description |
|------|-------------|
| `results_summary.csv` | All scenario metrics |

## Theory

### Simplicial Complexes

A **simplicial complex** K is a collection of simplices (vertices, edges, triangles, etc.) where:
- Every face of a simplex in K is also in K
- The intersection of any two simplices is a face of both

We construct two types:

1. **Structural Complex (K_str)**: Built from network graph using mesh-loop triangulation
2. **Functional Complex (K_func)**: Built from hydraulic head correlations using Vietoris-Rips

### Betti Numbers

- **β₀**: Number of connected components
- **β₁**: Number of independent loops (1-dimensional holes)
- **β₂**: Number of voids (2-dimensional holes)

### Functional Cohesion Index (FCI)

```
FCI = |Δ₂| × mean_overlap / (β₀ + 1)
```

Where:
- |Δ₂| = number of triangles
- mean_overlap = average triangles per node
- β₀ = connected components

### Redundancy Gap

```
Gap = (|Δ₂|_structural - |Δ₂|_functional) / |Δ₂|_structural × 100%
```

Represents physical loops that don't contribute to hydraulic connectivity.

## Citation

If you use this package in your research, please cite:

```bibtex
@article{wdsa_tda_2024,
  title={Topological Resilience Analysis of Water Distribution Systems},
  author={WDSA Research Team},
  journal={Journal of Water Resources Planning and Management},
  year={2024}
}
```

## License

MIT License

## Contact

For questions or issues, please open a GitHub issue or contact the authors.
