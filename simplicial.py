"""
================================================================================
simplicial.py - Simplicial Complex Construction for Water Distribution Networks
================================================================================

This module provides functions to build and analyze simplicial complexes from
water distribution network topology and hydraulic simulation data.

Two types of complexes are constructed:

1. STRUCTURAL COMPLEX (K_str):
   - Built from network graph topology using mesh-loop triangulation
   - Represents physically available redundancy paths
   - Independent of hydraulic conditions

2. FUNCTIONAL COMPLEX (K_func):
   - Built from hydraulic head correlations using Vietoris-Rips construction
   - Represents actual hydraulic connectivity under operating conditions
   - Depends on pressure/flow patterns

The difference between structural and functional complexes reveals the
"redundancy gap" - physical loops that don't contribute to hydraulic resilience.

References:
    - Edelsbrunner, H., & Harer, J. (2010). Computational Topology.
    - Ghrist, R. (2008). Barcodes: The persistent topology of data.

Author: WDSA Research Team
Date: 2024
"""

import numpy as np
import networkx as nx
from itertools import combinations


# =============================================================================
# STRUCTURAL COMPLEX CONSTRUCTION
# =============================================================================

def build_structural_complex(G, pos, max_loop_length=14):
    """
    Build structural 2-simplices using mesh-loop triangulation.
    
    This method identifies all independent cycles (loops) in the network graph,
    then triangulates each bounded cycle using fan triangulation from a fixed
    anchor vertex.
    
    Parameters
    ----------
    G : networkx.Graph
        Network graph where nodes are junctions/tanks/reservoirs and 
        edges are pipes/pumps/valves.
    pos : dict
        Node positions as {node_id: (x, y)} coordinates.
        Required for identifying the outer boundary and computing areas.
    max_loop_length : int, optional
        Maximum number of nodes in a cycle to consider for triangulation.
        Longer cycles are typically the outer network boundary.
        Default is 14.
    
    Returns
    -------
    triangles : list of tuples
        List of 2-simplices as (node_a, node_b, node_c) tuples.
    bounded_cycles : list of lists
        The bounded cycles that were triangulated.
    n_total_cycles : int
        Total number of independent cycles found in the graph.
    
    Algorithm
    ---------
    1. Find all independent cycles using NetworkX cycle_basis()
    2. Normalize cycles to canonical form for deduplication
    3. Compute polygon area for each cycle
    4. Identify outer boundary as largest-area cycle
    5. Triangulate remaining (bounded) cycles using fan triangulation
    
    Example
    -------
    >>> G = nx.Graph()
    >>> G.add_edges_from([(1,2), (2,3), (3,4), (4,1), (2,4)])
    >>> pos = {1: (0,0), 2: (1,0), 3: (1,1), 4: (0,1)}
    >>> triangles, cycles, n = build_structural_complex(G, pos)
    >>> print(f"Found {len(triangles)} triangles from {n} cycles")
    
    Notes
    -----
    - The outer boundary detection assumes the network is planar
    - Fan triangulation may create degenerate triangles for non-convex cycles
    - Cycles with nodes missing from `pos` are skipped
    """
    
    # Step 1: Find all independent cycles in the graph
    # cycle_basis returns a list of cycles forming a basis for the cycle space
    raw_cycles = nx.cycle_basis(G)
    
    def normalize_cycle(cycle):
        """
        Convert cycle to canonical form for comparison.
        
        Cycles can be represented starting from any node and going in
        either direction. This function creates a unique representation
        by starting from the minimum node ID and choosing the smaller
        of the two possible directions.
        
        Parameters
        ----------
        cycle : list
            List of node IDs forming a cycle.
        
        Returns
        -------
        tuple
            Canonical tuple representation of the cycle.
        """
        # Find position of minimum element
        min_pos = min(range(len(cycle)), key=lambda i: cycle[i])
        
        # Rotation starting from minimum
        rotation1 = tuple(cycle[min_pos:] + cycle[:min_pos])
        
        # Reverse direction, then rotate to start from minimum
        reversed_cycle = list(reversed(cycle[:min_pos] + cycle[min_pos:]))
        min_pos_rev = min(range(len(reversed_cycle)), key=lambda i: reversed_cycle[i])
        rotation2 = tuple(reversed_cycle[min_pos_rev:] + reversed_cycle[:min_pos_rev])
        
        # Return lexicographically smaller representation
        return min(rotation1, rotation2)
    
    # Step 2: Deduplicate cycles using canonical form
    seen_cycles = set()
    unique_cycles = []
    
    for cycle in raw_cycles:
        # Skip cycles that are too short or have missing position data
        if len(cycle) < 3:
            continue
        if not all(node in pos for node in cycle):
            continue
        
        canonical = normalize_cycle(cycle)
        if canonical not in seen_cycles:
            seen_cycles.add(canonical)
            unique_cycles.append(list(canonical))
    
    if not unique_cycles:
        return [], [], 0
    
    # Step 3: Compute polygon area for each cycle using shoelace formula
    def compute_polygon_area(node_list):
        """
        Compute area of polygon using the shoelace formula.
        
        Parameters
        ----------
        node_list : list
            List of node IDs forming a closed polygon.
        
        Returns
        -------
        float
            Absolute area of the polygon.
        """
        coords = [pos[n] for n in node_list]
        x_coords, y_coords = zip(*coords)
        n = len(coords)
        
        # Shoelace formula: A = 0.5 * |sum(x_i * y_{i+1} - x_{i+1} * y_i)|
        area = 0.0
        for i in range(n):
            j = (i + 1) % n
            area += x_coords[i] * y_coords[j]
            area -= x_coords[j] * y_coords[i]
        
        return 0.5 * abs(area)
    
    # Compute areas for all cycles
    cycle_areas = [(cycle, compute_polygon_area(cycle)) for cycle in unique_cycles]
    
    # Step 4: Identify outer boundary (largest area cycle)
    outer_boundary = max(cycle_areas, key=lambda x: x[1])[0]
    
    # Step 5: Filter to bounded cycles within length limit
    bounded_cycles = [
        cycle for cycle, area in cycle_areas 
        if cycle != outer_boundary and len(cycle) <= max_loop_length
    ]
    
    # Step 6: Triangulate each bounded cycle using fan triangulation
    # Fan triangulation: pick anchor vertex, connect to all other vertices
    triangles = []
    for cycle in bounded_cycles:
        anchor = cycle[0]
        for i in range(1, len(cycle) - 1):
            triangle = (anchor, cycle[i], cycle[i + 1])
            triangles.append(triangle)
    
    return triangles, bounded_cycles, len(unique_cycles)


# =============================================================================
# FUNCTIONAL COMPLEX CONSTRUCTION
# =============================================================================

def build_functional_complex(simulation_results, junction_list, 
                              target_triangle_range=(60, 300),
                              k_nearest_neighbors=3,
                              epsilon_percentiles=(1, 2, 3, 4, 5, 7, 10, 12, 15, 20)):
    """
    Build functional 2-simplices using Vietoris-Rips construction on hydraulic data.
    
    This method constructs a simplicial complex based on the similarity of
    hydraulic head time series between junction pairs. Junctions with similar
    head profiles (small median absolute difference) are connected.
    
    Parameters
    ----------
    simulation_results : wntr.sim.SimulationResults
        Results from WNTR hydraulic simulation containing node head data.
    junction_list : list
        List of junction node IDs to include in the complex.
    target_triangle_range : tuple of int, optional
        Desired range (min, max) for number of triangles.
        The algorithm searches for epsilon that achieves this.
        Default is (60, 300).
    k_nearest_neighbors : int, optional
        Maximum number of neighbors each node can connect to.
        Prevents overly dense complexes. Default is 3.
    epsilon_percentiles : tuple of float, optional
        Percentiles of pairwise distances to try as epsilon threshold.
        Default searches from 1st to 20th percentile.
    
    Returns
    -------
    epsilon : float
        The selected distance threshold (in meters of head difference).
    edges : set
        Set of edges as frozensets of node pairs.
    triangles : list of tuples
        List of 2-simplices as (node_a, node_b, node_c) tuples.
    
    Algorithm
    ---------
    1. Extract hydraulic head time series for all junctions
    2. Compute pairwise distance matrix using median absolute difference
    3. For each epsilon candidate:
       a. Build edges where distance <= epsilon (with kNN constraint)
       b. Find all 3-cliques (triangles) in the edge graph
       c. Check if triangle count is in target range
    4. Return first epsilon achieving target, or last tried
    
    Notes
    -----
    - Median absolute difference is robust to outliers from demand peaks
    - kNN constraint prevents hub nodes from dominating the complex
    - Epsilon is automatically tuned to achieve reasonable complex density
    
    Example
    -------
    >>> results = sim.run_sim()
    >>> junctions = wn.junction_name_list
    >>> eps, edges, tris = build_functional_complex(results, junctions)
    >>> print(f"ε* = {eps:.4f} m, {len(tris)} triangles")
    """
    
    # Step 1: Extract hydraulic head data
    head_data = simulation_results.node['head']
    head_junctions = head_data[junction_list].astype(float)
    n_junctions = len(junction_list)
    
    # Step 2: Compute pairwise distance matrix
    # Distance = median of absolute head differences over simulation period
    distance_matrix = np.zeros((n_junctions, n_junctions), dtype=float)
    
    for i in range(n_junctions):
        for j in range(i + 1, n_junctions):
            # Time series of absolute differences
            abs_diff = np.abs(head_junctions.iloc[:, i] - head_junctions.iloc[:, j])
            # Use median for robustness to demand-driven spikes
            median_diff = float(np.median(abs_diff))
            distance_matrix[i, j] = median_diff
            distance_matrix[j, i] = median_diff
    
    # Get non-zero upper triangle for percentile calculation
    upper_triangle = distance_matrix[np.triu_indices(n_junctions, k=1)]
    nonzero_distances = upper_triangle[upper_triangle > 0]
    if len(nonzero_distances) == 0:
        nonzero_distances = np.array([1.0])
    
    def build_rips_at_epsilon(epsilon):
        """
        Build Vietoris-Rips complex at given epsilon threshold.
        
        Parameters
        ----------
        epsilon : float
            Distance threshold for edge creation.
        
        Returns
        -------
        edges : set
            Set of edges as tuples.
        triangles : list
            List of triangle tuples.
        """
        edges = set()
        
        # For each node, connect to k nearest neighbors within epsilon
        for i in range(n_junctions):
            # Find neighbors within epsilon
            neighbors = [
                (distance_matrix[i, k], k) 
                for k in range(n_junctions) 
                if k != i and distance_matrix[i, k] <= epsilon
            ]
            # Sort by distance and take k nearest
            neighbors.sort()
            for _, k in neighbors[:k_nearest_neighbors]:
                edge = tuple(sorted((junction_list[i], junction_list[k])))
                edges.add(edge)
        
        # Find triangles (3-cliques)
        edge_set = set(edges)
        triangles = []
        
        for a, b, c in combinations(junction_list, 3):
            edge_ab = tuple(sorted((a, b)))
            edge_ac = tuple(sorted((a, c)))
            edge_bc = tuple(sorted((b, c)))
            
            if edge_ab in edge_set and edge_ac in edge_set and edge_bc in edge_set:
                triangles.append((a, b, c))
        
        return edges, triangles
    
    # Step 3: Search for optimal epsilon
    selected_result = None
    last_result = None
    
    for percentile in epsilon_percentiles:
        epsilon = np.percentile(nonzero_distances, percentile)
        edges, triangles = build_rips_at_epsilon(epsilon)
        last_result = (epsilon, edges, triangles)
        
        # Check if triangle count is in target range
        if target_triangle_range[0] <= len(triangles) <= target_triangle_range[1]:
            selected_result = last_result
            break
    
    # Return selected or last tried
    result = selected_result if selected_result else last_result
    return result[0], result[1], result[2]


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def extract_edges_from_triangles(triangles):
    """
    Extract the set of edges from a list of triangles.
    
    Parameters
    ----------
    triangles : list of tuples
        List of triangles as (a, b, c) tuples.
    
    Returns
    -------
    edges : list of tuples
        Sorted list of unique edges as (node1, node2) tuples.
    
    Example
    -------
    >>> tris = [(1, 2, 3), (2, 3, 4)]
    >>> edges = extract_edges_from_triangles(tris)
    >>> print(edges)
    [(1, 2), (1, 3), (2, 3), (2, 4), (3, 4)]
    """
    edge_set = set()
    
    for a, b, c in triangles:
        edge_set.add(tuple(sorted((a, b))))
        edge_set.add(tuple(sorted((b, c))))
        edge_set.add(tuple(sorted((a, c))))
    
    return sorted(edge_set)


def get_vertices_from_triangles(triangles):
    """
    Extract the set of vertices from a list of triangles.
    
    Parameters
    ----------
    triangles : list of tuples
        List of triangles as (a, b, c) tuples.
    
    Returns
    -------
    vertices : list
        Sorted list of unique vertex IDs.
    """
    vertex_set = set()
    
    for a, b, c in triangles:
        vertex_set.add(a)
        vertex_set.add(b)
        vertex_set.add(c)
    
    return sorted(vertex_set)
