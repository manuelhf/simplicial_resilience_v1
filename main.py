#!/usr/bin/env python3
"""
================================================================================
main.py - WDSA-TDA Analysis Pipeline
================================================================================

Main entry point for running the complete TDA-based resilience analysis
on a water distribution network.

This script demonstrates the full workflow:
1. Load network and run baseline simulation
2. Build structural and functional complexes
3. Define and run disruption scenarios
4. Generate publication-quality figures and tables

Usage:
    python main.py                          # Uses default Parete.inp
    python main.py --input network.inp      # Specify input file
    python main.py --output results/        # Specify output directory
    python main.py --pressure 15.0          # Change pressure threshold

Author: WDSA Research Team
Date: 2024
"""

import os
import argparse
import numpy as np
import pandas as pd
import networkx as nx
import wntr

# Import from our package
from wdsa_tda_package.simplicial import (
    build_structural_complex, 
    build_functional_complex,
    extract_edges_from_triangles
)
from wdsa_tda_package.topology import compute_complex_metrics
from wdsa_tda_package.scenarios import (
    identify_critical_components,
    define_standard_scenarios,
    run_all_scenarios
)
from wdsa_tda_package.visualization import (
    setup_publication_style,
    plot_network_comparison,
    plot_ensemble_separate_rows,
    plot_scenario_detail,
    close_all_figures
)
from wdsa_tda_package.tables import (
    generate_baseline_table,
    generate_scenario_table
)


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='TDA-based resilience analysis for water distribution networks',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python main.py
    python main.py --input MyNetwork.inp --output my_results/
    python main.py --pressure 25.0
        """
    )
    
    parser.add_argument(
        '--input', '-i',
        type=str,
        default='Parete.inp',
        help='Path to EPANET .inp file (default: Parete.inp)'
    )
    
    parser.add_argument(
        '--output', '-o',
        type=str,
        default='paper_outputs',
        help='Output directory for figures and tables (default: paper_outputs)'
    )
    
    parser.add_argument(
        '--pressure', '-p',
        type=float,
        default=20.0,
        help='Minimum pressure threshold in meters (default: 20.0)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Print detailed progress information'
    )
    
    return parser.parse_args()


def run_analysis(inp_filename, output_dir='paper_outputs', pressure_threshold=20.0, 
                 verbose=True):
    """
    Run complete TDA resilience analysis pipeline.
    
    Parameters
    ----------
    inp_filename : str
        Path to EPANET .inp file.
    output_dir : str
        Directory for output files (created if doesn't exist).
    pressure_threshold : float
        Minimum pressure for adequate service (meters).
    verbose : bool
        If True, print progress information.
    
    Returns
    -------
    results : dict
        Dictionary containing all computed metrics and results.
    """
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Setup plotting style
    setup_publication_style()
    
    if verbose:
        print("="*70)
        print("WDSA-TDA: Topological Resilience Analysis")
        print("="*70)
        print(f"\nInput file: {inp_filename}")
        print(f"Output directory: {output_dir}")
        print(f"Pressure threshold: {pressure_threshold} m")
    
    # =========================================================================
    # PHASE 1: BASELINE CHARACTERIZATION
    # =========================================================================
    
    if verbose:
        print("\n" + "="*70)
        print("PHASE 1: BASELINE CHARACTERIZATION")
        print("="*70)
    
    # Load network model
    if verbose:
        print("\n[1.1] Loading network model...")
    
    wn = wntr.network.WaterNetworkModel(inp_filename)
    
    # Configure for 24-hour baseline simulation
    wn.options.time.duration = 24 * 3600
    wn.options.time.hydraulic_timestep = 3600
    wn.options.time.report_timestep = 3600
    
    # Run baseline simulation
    if verbose:
        print("[1.2] Running baseline hydraulic simulation...")
    
    sim = wntr.sim.EpanetSimulator(wn)
    baseline_results = sim.run_sim()
    
    # Extract network information
    n_junctions = len(wn.junction_name_list)
    n_pipes = len(wn.pipe_name_list)
    
    if verbose:
        print(f"      Network: {n_junctions} junctions, {n_pipes} pipes")
    
    # Build network graph and get node positions
    G = nx.Graph(wntr.network.to_graph(wn))
    pos = {name: obj.coordinates for name, obj in wn.nodes()}
    junctions = [n for n, obj in wn.nodes() if obj.node_type == 'Junction']
    
    # -------------------------------------------------------------------------
    # Build Structural Complex
    # -------------------------------------------------------------------------
    
    if verbose:
        print("\n[1.3] Building structural complex...")
    
    structural_triangles, _, n_cycles = build_structural_complex(G, pos, max_loop_length=14)
    structural_vertices = list(G.nodes())
    structural_edges = extract_edges_from_triangles(structural_triangles)
    structural_metrics = compute_complex_metrics(
        structural_vertices, structural_edges, structural_triangles
    )
    
    if verbose:
        print(f"      Triangles: {structural_metrics['n_triangles']}")
        print(f"      β₀={structural_metrics['beta0']}, β₁={structural_metrics['beta1']}")
        print(f"      FCI: {structural_metrics['fci']:.2f}")
    
    # -------------------------------------------------------------------------
    # Build Functional Complex
    # -------------------------------------------------------------------------
    
    if verbose:
        print("\n[1.4] Building functional complex...")
    
    epsilon, functional_edges, functional_triangles = build_functional_complex(
        baseline_results, junctions
    )
    functional_metrics = compute_complex_metrics(
        junctions, list(functional_edges), functional_triangles
    )
    
    if verbose:
        print(f"      ε* = {epsilon:.4f} m")
        print(f"      Triangles: {functional_metrics['n_triangles']}")
        print(f"      β₀={functional_metrics['beta0']}, β₁={functional_metrics['beta1']}")
        print(f"      FCI: {functional_metrics['fci']:.2f}")
    
    # Calculate redundancy gap
    n_str = structural_metrics['n_triangles']
    n_func = functional_metrics['n_triangles']
    redundancy_gap = (n_str - n_func) / n_str * 100 if n_str > 0 else 0
    
    if verbose:
        print(f"\n      Redundancy gap: {redundancy_gap:.1f}%")
    
    # -------------------------------------------------------------------------
    # Generate Figure 1: Network Comparison
    # -------------------------------------------------------------------------
    
    if verbose:
        print("\n[1.5] Generating Figure 1 (baseline comparison)...")
    
    plot_network_comparison(
        G, pos, structural_triangles, functional_triangles,
        save_path=os.path.join(output_dir, 'Figure1_baseline.pdf')
    )
    
    # -------------------------------------------------------------------------
    # Generate Table 1: Baseline Metrics
    # -------------------------------------------------------------------------
    
    if verbose:
        print("[1.6] Generating Table 1 (baseline metrics)...")
    
    generate_baseline_table(
        structural_metrics, functional_metrics, epsilon,
        save_path=os.path.join(output_dir, 'Table1_baseline.tex')
    )
    
    # =========================================================================
    # PHASE 2: SCENARIO ANALYSIS
    # =========================================================================
    
    if verbose:
        print("\n" + "="*70)
        print("PHASE 2: SCENARIO ANALYSIS")
        print("="*70)
    
    # -------------------------------------------------------------------------
    # Identify Critical Components
    # -------------------------------------------------------------------------
    
    if verbose:
        print("\n[2.1] Identifying critical network components...")
    
    components = identify_critical_components(wn, baseline_results)
    
    if verbose:
        print(f"      Trunk mains: {components['trunk_mains'][:3]}...")
        print(f"      Secondary mains: {components['secondary_mains'][:3]}...")
    
    # -------------------------------------------------------------------------
    # Define Scenarios
    # -------------------------------------------------------------------------
    
    if verbose:
        print("\n[2.2] Defining disruption scenarios...")
    
    scenarios = define_standard_scenarios(components)
    
    if verbose:
        print(f"      Defined {len(scenarios)} scenarios: {list(scenarios.keys())}")
    
    # -------------------------------------------------------------------------
    # Run All Scenarios
    # -------------------------------------------------------------------------
    
    if verbose:
        print("\n[2.3] Running scenario simulations...")
    
    all_results = run_all_scenarios(
        inp_filename, scenarios, G, pos, pressure_threshold
    )
    
    # Get baseline data for reference
    baseline_data = all_results.get('S0')
    
    if baseline_data is None:
        print("ERROR: Baseline scenario (S0) failed!")
        return None
    
    # =========================================================================
    # PHASE 3: FIGURE GENERATION
    # =========================================================================
    
    if verbose:
        print("\n" + "="*70)
        print("PHASE 3: FIGURE GENERATION")
        print("="*70)
    
    # -------------------------------------------------------------------------
    # Figure 2: S1 Ensemble
    # -------------------------------------------------------------------------
    
    if verbose:
        print("\n[3.1] Generating Figure 2 (S1 ensemble)...")
    
    plot_ensemble_separate_rows(
        all_results, 'S1', baseline_data,
        save_path=os.path.join(output_dir, 'Figure2_S1_ensemble.pdf')
    )
    
    # -------------------------------------------------------------------------
    # Figure 3: S2 Ensemble
    # -------------------------------------------------------------------------
    
    if verbose:
        print("[3.2] Generating Figure 3 (S2 ensemble)...")
    
    plot_ensemble_separate_rows(
        all_results, 'S2', baseline_data,
        save_path=os.path.join(output_dir, 'Figure3_S2_ensemble.pdf')
    )
    
    # -------------------------------------------------------------------------
    # Individual Scenario Detail Figures
    # -------------------------------------------------------------------------
    
    if verbose:
        print("[3.3] Generating individual scenario detail figures...")
    
    for scenario_name in ['S1a', 'S1b', 'S1c', 'S2a', 'S2b', 'S2c']:
        if scenario_name in all_results:
            plot_scenario_detail(
                all_results[scenario_name], scenario_name, baseline_data,
                save_path=os.path.join(output_dir, f'Figure_{scenario_name}_detail.pdf')
            )
    
    # =========================================================================
    # PHASE 4: TABLE GENERATION
    # =========================================================================
    
    if verbose:
        print("\n" + "="*70)
        print("PHASE 4: TABLE GENERATION")
        print("="*70)
    
    # -------------------------------------------------------------------------
    # Table 2: S1 Results
    # -------------------------------------------------------------------------
    
    if verbose:
        print("\n[4.1] Generating Table 2 (S1 results)...")
    
    generate_scenario_table(
        all_results, 'S1', baseline_data,
        save_path=os.path.join(output_dir, 'Table2_S1_results.tex')
    )
    
    # -------------------------------------------------------------------------
    # Table 3: S2 Results
    # -------------------------------------------------------------------------
    
    if verbose:
        print("[4.2] Generating Table 3 (S2 results)...")
    
    generate_scenario_table(
        all_results, 'S2', baseline_data,
        save_path=os.path.join(output_dir, 'Table3_S2_results.tex')
    )
    
    # -------------------------------------------------------------------------
    # Save Results CSV
    # -------------------------------------------------------------------------
    
    if verbose:
        print("[4.3] Saving results summary CSV...")
    
    results_data = []
    for name, data in all_results.items():
        if name == 'S0':
            continue
        m = data['metrics']
        c = data['config']
        results_data.append({
            'scenario': name,
            'type': c['type'],
            'start_time': c['start_time'],
            'tri_retention_pct': m['tri_retention'],
            'min_service_pct': m['min_service'],
            'correlation_r': m['r_tri_svc']
        })
    
    results_df = pd.DataFrame(results_data)
    csv_path = os.path.join(output_dir, 'results_summary.csv')
    results_df.to_csv(csv_path, index=False)
    
    if verbose:
        print(f"  Saved: {csv_path}")
    
    # =========================================================================
    # PHASE 5: SUMMARY
    # =========================================================================
    
    if verbose:
        print("\n" + "="*70)
        print("SUMMARY")
        print("="*70)
        
        # Baseline summary
        bl_tri_min = baseline_data['timeline']['n_triangles_mesh'].min()
        bl_tri_max = baseline_data['timeline']['n_triangles_mesh'].max()
        bl_svc_min = baseline_data['hydraulic']['service_pct'].min()
        bl_svc_max = baseline_data['hydraulic']['service_pct'].max()
        
        print(f"\n  BASELINE TOPOLOGY:")
        print(f"    Structural: |Δ₂| = {structural_metrics['n_triangles']}, FCI = {structural_metrics['fci']:.2f}")
        print(f"    Functional: |Δ₂| = {functional_metrics['n_triangles']}, FCI = {functional_metrics['fci']:.2f}")
        print(f"    Redundancy gap: {redundancy_gap:.1f}%")
        
        print(f"\n  BASELINE OPERATIONAL RANGE (24h demand cycle):")
        print(f"    |Δ₂|: {bl_tri_min} – {bl_tri_max}")
        print(f"    Service: {bl_svc_min:.1f}% – {bl_svc_max:.1f}%")
        
        # Scenario summary
        print(f"\n  SCENARIO RESULTS:")
        print(f"  {'Scenario':<10} {'Type':<10} {'Tri_ret%':<10} {'Svc_min%':<10} {'r':<8}")
        print("  " + "-"*50)
        
        correlations = []
        for name in ['S1a', 'S1b', 'S1c', 'S2a', 'S2b', 'S2c']:
            if name in all_results:
                m = all_results[name]['metrics']
                c = all_results[name]['config']
                print(f"  {name:<10} {c['type']:<10} {m['tri_retention']:<10.1f} {m['min_service']:<10.1f} {m['r_tri_svc']:<8.3f}")
                correlations.append(m['r_tri_svc'])
        
        if correlations:
            print(f"\n  CORRELATION SUMMARY:")
            print(f"    Mean r = {np.mean(correlations):.3f} ± {np.std(correlations):.3f}")
            print(f"    Range: [{min(correlations):.3f}, {max(correlations):.3f}]")
        
        # Output files
        print(f"\n  OUTPUT FILES in {output_dir}/:")
        print("    Figures:")
        print("      • Figure1_baseline.pdf")
        print("      • Figure2_S1_ensemble.pdf")
        print("      • Figure3_S2_ensemble.pdf")
        print("      • Figure_S1a_detail.pdf ... Figure_S2c_detail.pdf")
        print("    Tables:")
        print("      • Table1_baseline.tex")
        print("      • Table2_S1_results.tex")
        print("      • Table3_S2_results.tex")
        print("    Data:")
        print("      • results_summary.csv")
        
        print("\n" + "="*70)
        print("ANALYSIS COMPLETE")
        print("="*70)
    
    # Clean up
    close_all_figures()
    
    # Return all results for programmatic access
    return {
        'structural_metrics': structural_metrics,
        'functional_metrics': functional_metrics,
        'epsilon': epsilon,
        'redundancy_gap': redundancy_gap,
        'scenarios': all_results,
        'baseline': baseline_data,
        'components': components
    }


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    import sys
    
    # Parse command-line arguments
    args = parse_arguments()
    
    # Check if input file exists
    if not os.path.exists(args.input):
        print(f"ERROR: Input file not found: {args.input}")
        sys.exit(1)
    
    # Run the analysis
    results = run_analysis(
        inp_filename=args.input,
        output_dir=args.output,
        pressure_threshold=args.pressure,
        verbose=args.verbose or True  # Always verbose unless changed
    )
    
    if results is None:
        print("Analysis failed!")
        sys.exit(1)
    
    sys.exit(0)
