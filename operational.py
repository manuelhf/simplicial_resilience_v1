"""
================================================================================
operational.py - Time-Varying Operational Complex Tracking
================================================================================

This module tracks the evolution of simplicial complexes during hydraulic
simulations. The "operational complex" represents the active portion of the
network where pressure is sufficient for service delivery.

KEY CONCEPT:
    At each timestep, nodes with pressure below the minimum threshold are
    considered "inactive" and removed from the complex. The operational
    complex is then rebuilt from the remaining active nodes and their
    connecting pipes.

This allows tracking how topological redundancy degrades during:
    - Pipe breaks and leaks
    - Pump failures
    - Demand surges
    - Cascading failures

The correlation between topological metrics (triangles, FCI) and hydraulic
performance (service coverage, pressure deficit) demonstrates that topology
can serve as an early warning indicator.

Author: WDSA Research Team
Date: 2024
"""

import numpy as np
import pandas as pd
import networkx as nx

from .simplicial import build_structural_complex
from .topology import compute_node_simplex_overlap, compute_fci


# =============================================================================
# OPERATIONAL COMPLEX TRACKING
# =============================================================================

def track_operational_complex(water_network, simulation_results, 
                               network_graph, node_positions,
                               pressure_threshold=20.0):
    """
    Track the operational simplicial complex over simulation time.
    
    At each timestep, this function:
    1. Identifies "active" nodes with pressure >= threshold
    2. Extracts the subgraph of active nodes
    3. Builds the operational complex from this subgraph
    4. Computes topological metrics (triangles, β₀, FCI)
    
    Parameters
    ----------
    water_network : wntr.network.WaterNetworkModel
        The WNTR water network model.
    simulation_results : wntr.sim.SimulationResults
        Results from WNTR hydraulic simulation.
    network_graph : networkx.Graph
        Full network graph (from wntr.network.to_graph).
    node_positions : dict
        Node coordinates as {node_id: (x, y)}.
    pressure_threshold : float, optional
        Minimum pressure (meters) for a node to be considered active.
        Default is 20.0 m, a common minimum service pressure.
    
    Returns
    -------
    timeline : pandas.DataFrame
        DataFrame with one row per timestep containing:
        - 'timestep': Simulation time in seconds
        - 'time_hours': Simulation time in hours
        - 'n_nodes': Number of active nodes
        - 'n_edges': Number of active edges
        - 'n_triangles_mesh': Number of triangles in operational complex
        - 'beta0': Connected components (≥1 when nodes exist)
        - 'fci': Functional Cohesion Index
    
    Notes
    -----
    - The pressure threshold of 20 m (≈28 psi) is a common engineering
      standard for minimum service pressure in distribution systems.
    - β₀ is set to minimum of 1 when active nodes exist, to avoid
      division issues and reflect that a single component exists.
    - Computational cost scales with simulation duration and network size.
    
    Example
    -------
    >>> results = sim.run_sim()
    >>> G = nx.Graph(wntr.network.to_graph(wn))
    >>> pos = {name: obj.coordinates for name, obj in wn.nodes()}
    >>> timeline = track_operational_complex(wn, results, G, pos)
    >>> print(timeline[['time_hours', 'n_triangles_mesh', 'fci']].head())
    """
    timeline_records = []
    
    # Get junction list for filtering
    junction_names = water_network.junction_name_list
    
    # Iterate over each timestep in the simulation
    for timestep in simulation_results.node['pressure'].index:
        
        # Get pressure at this timestep
        pressure = simulation_results.node['pressure'].loc[timestep]
        
        # Identify active nodes: junctions with pressure >= threshold
        active_nodes = [
            node for node in junction_names
            if node in pressure.index and pressure[node] >= pressure_threshold
        ]
        
        # Build subgraph of active nodes
        active_subgraph = network_graph.subgraph(active_nodes).copy()
        n_active_nodes = len(active_subgraph.nodes())
        n_active_edges = len(active_subgraph.edges())
        
        # Compute connected components
        if n_active_nodes > 0:
            n_components = nx.number_connected_components(active_subgraph)
            # Ensure β₀ is at least 1 when nodes exist
            n_components = max(1, n_components)
        else:
            n_components = 0
        
        # Build operational complex from active subgraph
        if n_active_nodes > 2:
            # Get positions for active nodes
            active_positions = {
                node: node_positions[node] 
                for node in active_subgraph.nodes() 
                if node in node_positions
            }
            
            # Build structural complex on active subgraph
            triangles, _, _ = build_structural_complex(
                active_subgraph, active_positions, max_loop_length=14
            )
            n_triangles = len(triangles)
            
            # Compute FCI
            if n_triangles > 0:
                overlap = compute_node_simplex_overlap(triangles)
                fci = compute_fci(triangles, overlap, n_components)
            else:
                fci = 0.0
        else:
            n_triangles = 0
            fci = 0.0
        
        # Record metrics for this timestep
        timeline_records.append({
            'timestep': timestep,
            'time_hours': timestep / 3600,
            'n_nodes': n_active_nodes,
            'n_edges': n_active_edges,
            'n_triangles_mesh': n_triangles,
            'beta0': n_components,
            'fci': fci
        })
    
    return pd.DataFrame(timeline_records)


def compute_hydraulic_metrics(water_network, simulation_results, 
                               pressure_threshold=20.0):
    """
    Compute hydraulic performance metrics over simulation time.
    
    Parameters
    ----------
    water_network : wntr.network.WaterNetworkModel
        The WNTR water network model.
    simulation_results : wntr.sim.SimulationResults
        Results from WNTR hydraulic simulation.
    pressure_threshold : float, optional
        Minimum pressure (meters) for adequate service.
        Default is 20.0 m.
    
    Returns
    -------
    metrics : pandas.DataFrame
        DataFrame with one row per timestep containing:
        - 'timestep': Simulation time in seconds
        - 'time_hours': Simulation time in hours
        - 'service_pct': Percentage of junctions with adequate pressure
        - 'pdi': Pressure Deficit Index (%)
    
    Metrics Explained
    -----------------
    SERVICE COVERAGE:
        Percentage of junction nodes with pressure >= threshold.
        100% means all customers have adequate pressure.
        
        service_pct = (nodes with P >= P_min) / (total nodes) × 100
    
    PRESSURE DEFICIT INDEX (PDI):
        Normalized measure of total pressure shortfall.
        0% means no deficit, higher values indicate worse conditions.
        
        PDI = Σ max(0, P_min - P_i) / (P_min × N) × 100
        
        where P_i is pressure at node i, and N is total nodes.
    
    Example
    -------
    >>> results = sim.run_sim()
    >>> hydraulic = compute_hydraulic_metrics(wn, results)
    >>> print(f"Min service: {hydraulic['service_pct'].min():.1f}%")
    """
    junction_names = water_network.junction_name_list
    n_junctions = len(junction_names)
    
    metric_records = []
    
    for timestep in simulation_results.node['pressure'].index:
        
        # Get pressures at all junctions
        pressure = simulation_results.node['pressure'].loc[timestep, junction_names]
        
        # Service coverage: % of nodes with adequate pressure
        n_served = (pressure >= pressure_threshold).sum()
        service_pct = n_served / n_junctions * 100
        
        # Pressure Deficit Index
        # Sum of deficits normalized by (threshold × node count)
        deficit = np.maximum(0, pressure_threshold - pressure)
        pdi = deficit.sum() / (pressure_threshold * n_junctions) * 100
        
        metric_records.append({
            'timestep': timestep,
            'time_hours': timestep / 3600,
            'service_pct': service_pct,
            'pdi': pdi
        })
    
    return pd.DataFrame(metric_records)


def merge_topology_hydraulics(topology_timeline, hydraulic_metrics):
    """
    Merge topological and hydraulic timelines for correlation analysis.
    
    Parameters
    ----------
    topology_timeline : pandas.DataFrame
        Output from track_operational_complex().
    hydraulic_metrics : pandas.DataFrame
        Output from compute_hydraulic_metrics().
    
    Returns
    -------
    merged : pandas.DataFrame
        Combined DataFrame with both topological and hydraulic columns,
        merged on 'time_hours'.
    
    Example
    -------
    >>> merged = merge_topology_hydraulics(topology, hydraulic)
    >>> correlation = merged['n_triangles_mesh'].corr(merged['service_pct'])
    """
    return pd.merge(topology_timeline, hydraulic_metrics, on='time_hours')
