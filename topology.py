"""
================================================================================
topology.py - Topological Invariants for Simplicial Complexes
================================================================================

This module computes topological invariants (Betti numbers) and derived metrics
(Functional Cohesion Index) for simplicial complexes built from water networks.

BETTI NUMBERS:
    β₀ : Number of connected components
    β₁ : Number of independent loops (1-dimensional holes)
    β₂ : Number of voids (2-dimensional holes, cavities)

These are computed using the rank-nullity theorem applied to boundary matrices
over the field Z/2Z (binary arithmetic).

FUNCTIONAL COHESION INDEX (FCI):
    A composite metric combining triangle count, node participation, and
    connectivity. Higher FCI indicates more robust topological redundancy.

    FCI = |Δ₂| × mean_overlap / (β₀ + 1)

    where:
    - |Δ₂| = number of 2-simplices (triangles)
    - mean_overlap = average number of triangles per participating node
    - β₀ = number of connected components

References:
    - Munkres, J. R. (1984). Elements of Algebraic Topology.
    - Hatcher, A. (2002). Algebraic Topology.

Author: WDSA Research Team
Date: 2024
"""

import numpy as np


# =============================================================================
# BETTI NUMBER COMPUTATION
# =============================================================================

def _row_reduce_mod2(matrix):
    """
    Perform Gaussian elimination over Z/2Z (binary field).
    
    This is used to compute the rank of boundary matrices for Betti number
    calculation. All arithmetic is modulo 2.
    
    Parameters
    ----------
    matrix : numpy.ndarray
        Binary matrix (dtype uint8 recommended).
    
    Returns
    -------
    rank : int
        The rank of the matrix over Z/2Z.
    
    Algorithm
    ---------
    Standard Gaussian elimination with:
    - Pivot selection from current row downward
    - Row swaps when needed
    - XOR operations for elimination (addition mod 2)
    
    Notes
    -----
    - The input matrix is copied and not modified
    - Time complexity: O(min(m,n) * m * n) for m×n matrix
    """
    # Work with a copy to avoid modifying original
    A = (matrix % 2).copy()
    n_rows, n_cols = A.shape
    
    rank = 0  # Current pivot row
    
    for col in range(n_cols):
        # Find pivot in current column
        pivot_row = None
        for row in range(rank, n_rows):
            if A[row, col] == 1:
                pivot_row = row
                break
        
        if pivot_row is None:
            # No pivot in this column, move to next
            continue
        
        # Swap pivot row into position
        if pivot_row != rank:
            A[[rank, pivot_row]] = A[[pivot_row, rank]]
        
        # Eliminate all other 1s in this column
        for row in range(n_rows):
            if row != rank and A[row, col] == 1:
                A[row, :] ^= A[rank, :]  # XOR = addition mod 2
        
        rank += 1
        
        if rank == n_rows:
            break
    
    return rank


def compute_betti_numbers(vertices, edges, triangles):
    """
    Compute Betti numbers β₀, β₁, β₂ for a 2-dimensional simplicial complex.
    
    Uses the boundary matrix method with rank computation over Z/2Z.
    
    Parameters
    ----------
    vertices : list
        List of vertex IDs (0-simplices).
    edges : list of tuples
        List of edges as (v1, v2) tuples (1-simplices).
    triangles : list of tuples
        List of triangles as (v1, v2, v3) tuples (2-simplices).
    
    Returns
    -------
    betti_numbers : tuple of int
        (β₀, β₁, β₂) - the Betti numbers.
    simplex_counts : tuple of int
        (n_vertices, n_edges, n_triangles) - simplex counts.
    
    Theory
    ------
    For a simplicial complex K with chain groups C_i and boundary maps ∂_i:
    
        C_2 --∂₂--> C_1 --∂₁--> C_0
    
    The Betti numbers are:
        β_i = dim(ker(∂_i)) - dim(im(∂_{i+1}))
            = dim(C_i) - rank(∂_i) - rank(∂_{i+1})
    
    Specifically:
        β₀ = |V| - rank(∂₁)           [connected components]
        β₁ = |E| - rank(∂₁) - rank(∂₂) [independent loops]
        β₂ = |T| - rank(∂₂)           [voids/cavities]
    
    Example
    -------
    >>> V = [1, 2, 3, 4]
    >>> E = [(1,2), (2,3), (3,1), (3,4)]
    >>> T = [(1, 2, 3)]
    >>> betti, counts = compute_betti_numbers(V, E, T)
    >>> print(f"β₀={betti[0]}, β₁={betti[1]}, β₂={betti[2]}")
    """
    # Normalize and deduplicate simplices
    V = sorted(set(vertices))
    E = sorted({tuple(sorted(e)) for e in edges})
    T = sorted({tuple(sorted(t)) for t in triangles})
    
    # Create index mappings
    vertex_to_idx = {v: i for i, v in enumerate(V)}
    edge_to_idx = {e: i for i, e in enumerate(E)}
    
    n_vertices = len(V)
    n_edges = len(E)
    n_triangles = len(T)
    
    # Build boundary matrix ∂₁: C_1 → C_0
    # Entry (v, e) = 1 if vertex v is in edge e
    boundary_1 = np.zeros((n_vertices, n_edges), dtype=np.uint8)
    
    for j, (u, v) in enumerate(E):
        boundary_1[vertex_to_idx[u], j] = 1
        boundary_1[vertex_to_idx[v], j] = 1
    
    # Build boundary matrix ∂₂: C_2 → C_1
    # Entry (e, t) = 1 if edge e is in triangle t
    boundary_2 = np.zeros((n_edges, n_triangles), dtype=np.uint8)
    
    for j, (a, b, c) in enumerate(T):
        # Each triangle has 3 edges
        for edge in [tuple(sorted((a, b))), 
                     tuple(sorted((b, c))), 
                     tuple(sorted((a, c)))]:
            if edge in edge_to_idx:
                boundary_2[edge_to_idx[edge], j] = 1
    
    # Compute ranks
    rank_d1 = _row_reduce_mod2(boundary_1) if n_edges > 0 else 0
    rank_d2 = _row_reduce_mod2(boundary_2) if n_triangles > 0 else 0
    
    # Compute Betti numbers
    beta_0 = n_vertices - rank_d1      # Connected components
    beta_1 = n_edges - rank_d1 - rank_d2  # Independent loops
    beta_2 = n_triangles - rank_d2     # Voids
    
    return (beta_0, beta_1, beta_2), (n_vertices, n_edges, n_triangles)


# =============================================================================
# DERIVED METRICS
# =============================================================================

def compute_node_simplex_overlap(triangles):
    """
    Count how many triangles each node participates in.
    
    This measures the local redundancy at each node. Nodes participating
    in many triangles have more alternative flow paths available.
    
    Parameters
    ----------
    triangles : list of tuples
        List of triangles as (a, b, c) tuples.
    
    Returns
    -------
    overlap : dict
        Dictionary mapping node ID to triangle count.
    
    Example
    -------
    >>> tris = [(1,2,3), (1,3,4), (1,4,5)]
    >>> overlap = compute_node_simplex_overlap(tris)
    >>> print(overlap[1])  # Node 1 is in all 3 triangles
    3
    """
    overlap = {}
    
    for a, b, c in triangles:
        for node in (a, b, c):
            overlap[node] = overlap.get(node, 0) + 1
    
    return overlap


def compute_fci(triangles, overlap_dict, beta_0):
    """
    Compute the Functional Cohesion Index (FCI).
    
    FCI is a composite metric that captures:
    - Topological redundancy (number of triangles)
    - Local connectivity (node participation in triangles)
    - Global connectivity (inversely related to fragmentation)
    
    Formula:
        FCI = |Δ₂| × mean_overlap / (β₀ + 1)
    
    Parameters
    ----------
    triangles : list of tuples
        List of triangles in the complex.
    overlap_dict : dict
        Node-to-triangle-count mapping from compute_node_simplex_overlap().
    beta_0 : int
        Number of connected components (Betti-0).
    
    Returns
    -------
    fci : float
        The Functional Cohesion Index. Higher is better.
    
    Interpretation
    --------------
    - FCI increases with more triangles (more redundancy)
    - FCI increases when nodes participate in many triangles (tighter mesh)
    - FCI decreases as network fragments (β₀ increases)
    
    The +1 in denominator prevents division by zero and provides
    smooth behavior as the network approaches full connectivity.
    
    Example
    -------
    >>> triangles = [(1,2,3), (2,3,4), (3,4,5)]
    >>> overlap = compute_node_simplex_overlap(triangles)
    >>> beta0 = 1
    >>> fci = compute_fci(triangles, overlap, beta0)
    >>> print(f"FCI = {fci:.2f}")
    """
    if not triangles or not overlap_dict:
        return 0.0
    
    n_triangles = len(triangles)
    mean_overlap = np.mean(list(overlap_dict.values()))
    
    # FCI formula
    fci = n_triangles * mean_overlap / (beta_0 + 1)
    
    return fci


def compute_complex_metrics(vertices, edges, triangles):
    """
    Compute all topological metrics for a simplicial complex.
    
    This is a convenience function that computes Betti numbers, overlap,
    and FCI in one call.
    
    Parameters
    ----------
    vertices : list
        List of vertex IDs.
    edges : list of tuples
        List of edges as (v1, v2) tuples.
    triangles : list of tuples
        List of triangles as (v1, v2, v3) tuples.
    
    Returns
    -------
    metrics : dict
        Dictionary containing:
        - 'n_vertices': Number of vertices
        - 'n_edges': Number of edges
        - 'n_triangles': Number of triangles
        - 'beta0': Connected components (β₀)
        - 'beta1': Independent loops (β₁)
        - 'beta2': Voids (β₂)
        - 'fci': Functional Cohesion Index
        - 'overlap': Node overlap dictionary
    
    Example
    -------
    >>> metrics = compute_complex_metrics(V, E, T)
    >>> print(f"|Δ₂| = {metrics['n_triangles']}, FCI = {metrics['fci']:.2f}")
    """
    # Compute Betti numbers
    betti, counts = compute_betti_numbers(vertices, edges, triangles)
    beta0, beta1, beta2 = betti
    n_vertices, n_edges, n_triangles = counts
    
    # Compute overlap and FCI
    overlap = compute_node_simplex_overlap(triangles)
    fci = compute_fci(triangles, overlap, beta0)
    
    return {
        'n_vertices': n_vertices,
        'n_edges': n_edges,
        'n_triangles': n_triangles,
        'beta0': beta0,
        'beta1': beta1,
        'beta2': beta2,
        'fci': fci,
        'overlap': overlap
    }
