"""
================================================================================
visualization.py - Publication-Quality Figure Generation
================================================================================

This module provides functions for creating publication-quality figures
for water distribution network resilience analysis.

Figure types:
    - Network topology comparison (graph, structural, functional complex)
    - Scenario ensemble plots (multiple scenarios, separate panels)
    - Individual scenario detail plots
    - Correlation analysis plots

All figures are designed without titles (for LaTeX caption placement)
and with increased font sizes for readability.

Author: WDSA Research Team
Date: 2024
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import networkx as nx


# =============================================================================
# PLOT SETTINGS
# =============================================================================

def setup_publication_style():
    """
    Configure matplotlib for publication-quality figures.
    
    Sets serif fonts, increased font sizes, and clean styling.
    Call this at the start of your script.
    """
    plt.rcParams.update({
        'font.family': 'serif',
        'font.size': 12,
        'axes.labelsize': 14,
        'axes.titlesize': 14,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 11,
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'axes.grid': True,
        'grid.alpha': 0.3,
        'grid.linewidth': 0.5,
        'lines.linewidth': 2
    })


# Color scheme for scenarios
COLORS = {
    'S1a': '#1B9E77',  # Teal
    'S1b': '#D95F02',  # Orange
    'S1c': '#7570B3',  # Purple
    'S2a': '#1B9E77',
    'S2b': '#D95F02',
    'S2c': '#7570B3',
    'baseline': '#888888',
    'structural': '#377EB8',
    'functional': '#E6550D'
}


# =============================================================================
# BASELINE VISUALIZATION
# =============================================================================

def plot_network_comparison(graph, positions, structural_simplices, 
                            functional_simplices, save_path=None):
    """
    Create 3-panel figure comparing network graph and simplicial complexes.
    
    Panels:
        (a) Network graph - nodes and edges only
        (b) Structural complex - with 2-simplices shaded blue
        (c) Functional complex - with 2-simplices shaded orange
    
    Parameters
    ----------
    graph : networkx.Graph
        Network graph.
    positions : dict
        Node coordinates {node_id: (x, y)}.
    structural_simplices : list of tuples
        Triangles from structural complex.
    functional_simplices : list of tuples
        Triangles from functional complex.
    save_path : str, optional
        If provided, save figure to this path.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object.
    
    Example
    -------
    >>> fig = plot_network_comparison(G, pos, struct_tris, func_tris, 
    ...                               save_path='figure1.pdf')
    """
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    
    node_size = 15
    edge_width = 0.5
    
    # Panel (a): Network graph
    ax = axes[0]
    nx.draw_networkx_edges(graph, pos=positions, ax=ax, 
                           width=edge_width, alpha=0.6, edge_color='gray')
    nx.draw_networkx_nodes(graph, pos=positions, ax=ax, 
                           node_size=node_size, node_color='black')
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Panel (b): Structural complex
    ax = axes[1]
    # Draw triangles first (underneath)
    polygons = [
        [positions[a], positions[b], positions[c]] 
        for (a, b, c) in structural_simplices
        if a in positions and b in positions and c in positions
    ]
    nx.draw_networkx_edges(graph, pos=positions, ax=ax, 
                           width=edge_width, alpha=0.4, edge_color='gray')
    nx.draw_networkx_nodes(graph, pos=positions, ax=ax, 
                           node_size=node_size, node_color='black', alpha=0.7)
    if polygons:
        collection = PolyCollection(polygons, closed=True,
                                    facecolors=COLORS['structural'],
                                    edgecolors='none', alpha=0.3)
        ax.add_collection(collection)
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Panel (c): Functional complex
    ax = axes[2]
    polygons = [
        [positions[a], positions[b], positions[c]] 
        for (a, b, c) in functional_simplices
        if a in positions and b in positions and c in positions
    ]
    nx.draw_networkx_edges(graph, pos=positions, ax=ax, 
                           width=edge_width, alpha=0.4, edge_color='gray')
    nx.draw_networkx_nodes(graph, pos=positions, ax=ax, 
                           node_size=node_size, node_color='black', alpha=0.7)
    if polygons:
        collection = PolyCollection(polygons, closed=True,
                                    facecolors=COLORS['functional'],
                                    edgecolors='none', alpha=0.35)
        ax.add_collection(collection)
    ax.set_aspect('equal')
    ax.axis('off')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"  Saved: {save_path}")
    
    return fig


# =============================================================================
# SCENARIO ENSEMBLE VISUALIZATION
# =============================================================================

def plot_ensemble_separate_rows(all_results, scenario_prefix, 
                                 baseline_data, save_path=None):
    """
    Create ensemble figure with separate row for each scenario variant.
    
    Layout: 3 rows × 3 columns
        Rows: scenario variants (a, b, c for different start times)
        Columns: Triangles, FCI, Service coverage
    
    Each panel shows the full 24-hour timeline with:
        - Gray dashed line: baseline for reference
        - Colored line: scenario data
        - Red vertical line: disruption time
    
    Parameters
    ----------
    all_results : dict
        Dictionary of scenario results from run_all_scenarios().
    scenario_prefix : str
        'S1' or 'S2' to select which scenarios to plot.
    baseline_data : dict
        Results for baseline scenario (S0).
    save_path : str, optional
        If provided, save figure to this path.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object.
    
    Example
    -------
    >>> fig = plot_ensemble_separate_rows(results, 'S1', results['S0'],
    ...                                   save_path='figure2.pdf')
    """
    scenarios = [f'{scenario_prefix}a', f'{scenario_prefix}b', f'{scenario_prefix}c']
    row_labels = ['02h (min demand)', '08h (peak demand)', '14h (mid demand)']
    start_times = [2, 8, 14]
    
    # Baseline data for reference
    bl_timeline = baseline_data['timeline']
    bl_hydraulic = baseline_data['hydraulic']
    
    fig, axes = plt.subplots(3, 3, figsize=(14, 10), sharex=True)
    
    for row_idx, (scenario, row_label, start_h) in enumerate(
            zip(scenarios, row_labels, start_times)):
        
        if scenario not in all_results:
            continue
        
        data = all_results[scenario]
        timeline = data['timeline']
        hydraulic = data['hydraulic']
        color = COLORS[scenario]
        r_value = data['metrics']['r_tri_svc']
        
        # Column 0: Operational 2-simplices
        ax = axes[row_idx, 0]
        ax.plot(bl_timeline['time_hours'], bl_timeline['n_triangles_mesh'],
                color=COLORS['baseline'], linewidth=1.5, linestyle='--',
                alpha=0.7, label='Baseline')
        ax.plot(timeline['time_hours'], timeline['n_triangles_mesh'],
                color=color, linewidth=2, label=scenario)
        ax.axvline(x=start_h, color='red', linestyle='-', linewidth=2, alpha=0.7)
        ax.set_ylabel(f'{scenario}\n({row_label})\n\n$|\\Delta_2|$')
        ax.set_ylim(0, 300)
        ax.legend(loc='upper right', fontsize=10)
        
        # Column 1: FCI
        ax = axes[row_idx, 1]
        ax.plot(bl_timeline['time_hours'], bl_timeline['fci'],
                color=COLORS['baseline'], linewidth=1.5, linestyle='--', alpha=0.7)
        ax.plot(timeline['time_hours'], timeline['fci'],
                color=color, linewidth=2)
        ax.axvline(x=start_h, color='red', linestyle='-', linewidth=2, alpha=0.7)
        ax.set_ylim(0, 700)
        
        # Column 2: Service coverage
        ax = axes[row_idx, 2]
        ax.plot(bl_hydraulic['time_hours'], bl_hydraulic['service_pct'],
                color=COLORS['baseline'], linewidth=1.5, linestyle='--', alpha=0.7)
        ax.plot(hydraulic['time_hours'], hydraulic['service_pct'],
                color=color, linewidth=2)
        ax.axvline(x=start_h, color='red', linestyle='-', linewidth=2, alpha=0.7)
        ax.set_ylim(0, 105)
        
        # Correlation annotation
        ax.text(0.98, 0.05, f'r = {r_value:.3f}',
                transform=ax.transAxes, ha='right', va='bottom',
                fontsize=12, fontweight='bold',
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
    
    # X-axis labels on bottom row only
    for col in range(3):
        axes[2, col].set_xlabel('Time (hours)')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"  Saved: {save_path}")
    
    return fig


# =============================================================================
# INDIVIDUAL SCENARIO DETAIL VISUALIZATION
# =============================================================================

def plot_scenario_detail(scenario_data, scenario_name, baseline_data, 
                         save_path=None):
    """
    Create detailed 2×2 figure for a single scenario.
    
    Panels:
        (a) Triangles over time with baseline shaded
        (b) Service coverage over time with baseline shaded
        (c) FCI over time with baseline shaded
        (d) Scatter plot showing topology-service correlation
    
    Parameters
    ----------
    scenario_data : dict
        Results for this scenario from run_scenario().
    scenario_name : str
        Scenario identifier (e.g., 'S1a').
    baseline_data : dict
        Results for baseline scenario (S0).
    save_path : str, optional
        If provided, save figure to this path.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object.
    
    Example
    -------
    >>> fig = plot_scenario_detail(results['S1a'], 'S1a', results['S0'],
    ...                            save_path='figure_S1a.pdf')
    """
    config = scenario_data['config']
    start_h = config['start_time']
    metrics = scenario_data['metrics']
    
    timeline = scenario_data['timeline']
    hydraulic = scenario_data['hydraulic']
    merged = scenario_data['merged']
    
    bl_timeline = baseline_data['timeline']
    bl_hydraulic = baseline_data['hydraulic']
    
    color = COLORS.get(scenario_name, '#1B9E77')
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    
    # Panel (a): Triangles
    ax = axes[0, 0]
    ax.fill_between(bl_timeline['time_hours'], 0, bl_timeline['n_triangles_mesh'],
                    alpha=0.15, color='gray', label='Baseline range')
    ax.plot(timeline['time_hours'], timeline['n_triangles_mesh'],
            color=color, linewidth=2.5, label=scenario_name)
    ax.axvline(x=start_h, color='red', linestyle='-', linewidth=2,
               alpha=0.8, label='Disruption')
    ax.set_ylabel('Operational 2-simplices ($|\\Delta_2|$)')
    ax.legend(loc='upper right')
    ax.set_xlim(0, 24)
    ax.set_ylim(0, 300)
    
    # Panel (b): Service coverage
    ax = axes[0, 1]
    ax.fill_between(bl_hydraulic['time_hours'], 0, bl_hydraulic['service_pct'],
                    alpha=0.15, color='gray')
    ax.plot(hydraulic['time_hours'], hydraulic['service_pct'],
            color=color, linewidth=2.5)
    ax.axvline(x=start_h, color='red', linestyle='-', linewidth=2, alpha=0.8)
    ax.set_ylabel('Service coverage (%)')
    ax.set_xlim(0, 24)
    ax.set_ylim(0, 105)
    
    # Panel (c): FCI
    ax = axes[1, 0]
    ax.fill_between(bl_timeline['time_hours'], 0, bl_timeline['fci'],
                    alpha=0.15, color='gray')
    ax.plot(timeline['time_hours'], timeline['fci'],
            color=color, linewidth=2.5)
    ax.axvline(x=start_h, color='red', linestyle='-', linewidth=2, alpha=0.8)
    ax.set_xlabel('Time (hours)')
    ax.set_ylabel('Functional Cohesion Index')
    ax.set_xlim(0, 24)
    ax.set_ylim(0, 700)
    
    # Panel (d): Correlation scatter
    ax = axes[1, 1]
    
    pre_mask = merged['time_hours'] < start_h
    post_mask = merged['time_hours'] >= start_h
    
    ax.scatter(merged.loc[pre_mask, 'n_triangles_mesh'],
               merged.loc[pre_mask, 'service_pct'],
               c='gray', alpha=0.5, s=40, label='Pre-disruption', edgecolors='none')
    ax.scatter(merged.loc[post_mask, 'n_triangles_mesh'],
               merged.loc[post_mask, 'service_pct'],
               c=color, alpha=0.7, s=40, label='Post-disruption', edgecolors='none')
    
    ax.set_xlabel('Operational 2-simplices ($|\\Delta_2|$)')
    ax.set_ylabel('Service coverage (%)')
    ax.legend(loc='lower right')
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 105)
    
    # Add correlation text
    ax.text(0.05, 0.95, f'r = {metrics["r_tri_svc"]:.3f}',
            transform=ax.transAxes, ha='left', va='top',
            fontsize=14, fontweight='bold',
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"  Saved: {save_path}")
    
    return fig


def close_all_figures():
    """Close all open matplotlib figures to free memory."""
    plt.close('all')
