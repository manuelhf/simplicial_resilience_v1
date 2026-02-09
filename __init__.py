"""
================================================================================
WDSA-TDA: Topological Data Analysis for Water Distribution System Resilience
================================================================================

A Python package for analyzing water distribution network resilience using
topological data analysis (TDA) methods, specifically simplicial complexes.

This package provides tools for:

1. SIMPLICIAL COMPLEX CONSTRUCTION
   - Structural complex from network topology (mesh-loop triangulation)
   - Functional complex from hydraulic correlations (Vietoris-Rips)

2. TOPOLOGICAL INVARIANTS
   - Betti numbers (β₀, β₁, β₂)
   - Functional Cohesion Index (FCI)

3. OPERATIONAL TRACKING
   - Time-varying operational complex during simulations
   - Correlation with hydraulic performance metrics

4. SCENARIO ANALYSIS
   - Disruption scenario definition and execution
   - Comparison across demand levels

5. VISUALIZATION
   - Publication-quality figures
   - LaTeX table generation

Quick Start
-----------
>>> from wdsa_tda import run_analysis
>>> results = run_analysis('network.inp', output_dir='outputs')

Or use individual modules:
>>> from wdsa_tda.simplicial import build_structural_complex
>>> from wdsa_tda.topology import compute_betti_numbers
>>> from wdsa_tda.scenarios import run_scenario

Requirements
------------
- wntr >= 0.5.0
- networkx >= 2.6
- numpy >= 1.20
- pandas >= 1.3
- matplotlib >= 3.4
- scipy >= 1.7

References
----------
[1] Edelsbrunner, H., & Harer, J. (2010). Computational Topology.
[2] Ghrist, R. (2008). Barcodes: The persistent topology of data.
[3] WNTR: Water Network Tool for Resilience (US EPA)

Author: WDSA Research Team
License: MIT
Date: 2024
"""

__version__ = "1.0.0"
__author__ = "WDSA Research Team"

# Import main modules
from . import simplicial
from . import topology
from . import operational
from . import scenarios
from . import visualization
from . import tables

# Import commonly used functions for convenience
from .simplicial import (
    build_structural_complex,
    build_functional_complex,
    extract_edges_from_triangles
)

from .topology import (
    compute_betti_numbers,
    compute_fci,
    compute_complex_metrics
)

from .operational import (
    track_operational_complex,
    compute_hydraulic_metrics
)

from .scenarios import (
    identify_critical_components,
    define_standard_scenarios,
    run_scenario,
    run_all_scenarios
)

from .visualization import (
    setup_publication_style,
    plot_network_comparison,
    plot_ensemble_separate_rows,
    plot_scenario_detail
)

from .tables import (
    generate_baseline_table,
    generate_scenario_table
)

# Define what gets imported with "from wdsa_tda import *"
__all__ = [
    # Modules
    'simplicial',
    'topology', 
    'operational',
    'scenarios',
    'visualization',
    'tables',
    
    # Simplicial functions
    'build_structural_complex',
    'build_functional_complex',
    'extract_edges_from_triangles',
    
    # Topology functions
    'compute_betti_numbers',
    'compute_fci',
    'compute_complex_metrics',
    
    # Operational functions
    'track_operational_complex',
    'compute_hydraulic_metrics',
    
    # Scenario functions
    'identify_critical_components',
    'define_standard_scenarios',
    'run_scenario',
    'run_all_scenarios',
    
    # Visualization functions
    'setup_publication_style',
    'plot_network_comparison',
    'plot_ensemble_separate_rows',
    'plot_scenario_detail',
    
    # Table functions
    'generate_baseline_table',
    'generate_scenario_table',
]
