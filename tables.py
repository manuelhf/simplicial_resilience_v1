"""
================================================================================
tables.py - LaTeX Table Generation for Publication
================================================================================

This module provides functions for generating publication-ready LaTeX tables
summarizing baseline topology and scenario results.

Tables are formatted for direct inclusion in LaTeX documents using the
booktabs package for professional styling.

Author: WDSA Research Team
Date: 2024
"""


def generate_baseline_table(structural_metrics, functional_metrics, 
                            epsilon, save_path=None):
    """
    Generate LaTeX table comparing structural and functional complexes.
    
    Parameters
    ----------
    structural_metrics : dict
        Metrics from compute_complex_metrics() for structural complex.
    functional_metrics : dict
        Metrics from compute_complex_metrics() for functional complex.
    epsilon : float
        The optimal epsilon threshold used for functional complex.
    save_path : str, optional
        If provided, save LaTeX to this file.
    
    Returns
    -------
    latex : str
        LaTeX table code.
    
    Table Contents
    --------------
    - 2-simplices count
    - β₀ (connected components)
    - β₁ (independent loops)
    - FCI (Functional Cohesion Index)
    - ε* (optimal threshold)
    - Redundancy gap percentage
    
    Example
    -------
    >>> latex = generate_baseline_table(struct_m, func_m, eps, 
    ...                                 save_path='table1.tex')
    """
    # Calculate redundancy gap
    n_str = structural_metrics['n_triangles']
    n_func = functional_metrics['n_triangles']
    gap = (n_str - n_func) / n_str * 100 if n_str > 0 else 0
    
    latex = f"""
\\begin{{table}}[htbp]
\\centering
\\caption{{Baseline topological characterisation. The redundancy gap indicates 
the fraction of structural triangles that do not contribute to functional 
hydraulic connectivity.}}
\\label{{tab:baseline}}
\\begin{{tabular}}{{lcc}}
\\toprule
\\textbf{{Metric}} & \\textbf{{Structural}} & \\textbf{{Functional}} \\\\
\\midrule
2-simplices ($|\\Delta_2|$) & {structural_metrics['n_triangles']} & {functional_metrics['n_triangles']} \\\\
$\\beta_0$ (components) & {structural_metrics['beta0']} & {functional_metrics['beta0']} \\\\
$\\beta_1$ (loops) & {structural_metrics['beta1']} & {functional_metrics['beta1']} \\\\
FCI & {structural_metrics['fci']:.2f} & {functional_metrics['fci']:.2f} \\\\
\\midrule
$\\varepsilon^*$ (m) & -- & {epsilon:.4f} \\\\
Redundancy gap & \\multicolumn{{2}}{{c}}{{{gap:.1f}\\%}} \\\\
\\bottomrule
\\end{{tabular}}
\\end{{table}}
"""
    
    if save_path:
        with open(save_path, 'w') as f:
            f.write(latex)
        print(f"  Saved: {save_path}")
    
    return latex


def generate_scenario_table(all_results, scenario_prefix, baseline_data,
                            save_path=None):
    """
    Generate LaTeX table summarizing scenario ensemble results.
    
    Parameters
    ----------
    all_results : dict
        Dictionary of scenario results from run_all_scenarios().
    scenario_prefix : str
        'S1' or 'S2' to select which scenarios to include.
    baseline_data : dict
        Results for baseline scenario (S0).
    save_path : str, optional
        If provided, save LaTeX to this file.
    
    Returns
    -------
    latex : str
        LaTeX table code.
    
    Table Contents
    --------------
    For each scenario variant:
    - Start time
    - Mean post-disruption triangles (absolute and % of baseline)
    - Mean post-disruption service (absolute and % of baseline)
    - Correlation coefficient r
    
    Example
    -------
    >>> latex = generate_scenario_table(results, 'S1', results['S0'],
    ...                                 save_path='table2.tex')
    """
    scenarios = [f'{scenario_prefix}a', f'{scenario_prefix}b', f'{scenario_prefix}c']
    
    # Baseline reference values (minimum over 24h cycle)
    bl_tri_min = baseline_data['timeline']['n_triangles_mesh'].min()
    bl_svc_min = baseline_data['hydraulic']['service_pct'].min()
    
    # Build table rows
    rows = ""
    for scenario in scenarios:
        if scenario not in all_results:
            continue
        
        data = all_results[scenario]
        config = data['config']
        metrics = data['metrics']
        timeline = data['timeline']
        hydraulic = data['hydraulic']
        
        start_h = config['start_time']
        
        # Post-disruption statistics
        post_timeline = timeline[timeline['time_hours'] >= start_h]
        post_hydraulic = hydraulic[hydraulic['time_hours'] >= start_h]
        
        tri_mean = post_timeline['n_triangles_mesh'].mean()
        svc_mean = post_hydraulic['service_pct'].mean()
        
        # Normalized to baseline minimum
        tri_pct = (tri_mean / bl_tri_min * 100) if bl_tri_min > 0 else 0
        svc_pct = (svc_mean / bl_svc_min * 100) if bl_svc_min > 0 else 0
        
        rows += f"    {scenario} & {start_h:02d}:00 & {tri_mean:.0f} & {tri_pct:.0f} & {svc_mean:.1f} & {svc_pct:.0f} & {metrics['r_tri_svc']:.3f} \\\\\n"
    
    # Scenario type description
    if scenario_prefix == 'S1':
        scenario_type = "secondary main break"
    else:
        scenario_type = "trunk main break"
    
    latex = f"""
\\begin{{table}}[htbp]
\\centering
\\caption{{{scenario_prefix} ensemble results: {scenario_type} at three demand levels.
Mean post-disruption values shown with normalization to baseline minimum
($|\\Delta_2|_{{\\min}}$ = {bl_tri_min}, Service$_{{\\min}}$ = {bl_svc_min:.1f}\\%).}}
\\label{{tab:{scenario_prefix.lower()}_results}}
\\begin{{tabular}}{{lcccccc}}
\\toprule
\\textbf{{ID}} & \\textbf{{Start}} & \\textbf{{$|\\Delta_2|$}} & \\textbf{{(\\%)}} & \\textbf{{Svc (\\%)}} & \\textbf{{(\\%)}} & \\textbf{{$r$}} \\\\
\\midrule
{rows}\\bottomrule
\\end{{tabular}}
\\end{{table}}
"""
    
    if save_path:
        with open(save_path, 'w') as f:
            f.write(latex)
        print(f"  Saved: {save_path}")
    
    return latex


def generate_summary_table(all_results, save_path=None):
    """
    Generate LaTeX table with summary statistics for all scenarios.
    
    Parameters
    ----------
    all_results : dict
        Dictionary of all scenario results.
    save_path : str, optional
        If provided, save LaTeX to this file.
    
    Returns
    -------
    latex : str
        LaTeX table code.
    """
    rows = ""
    
    for name in ['S1a', 'S1b', 'S1c', 'S2a', 'S2b', 'S2c']:
        if name not in all_results:
            continue
        
        data = all_results[name]
        config = data['config']
        metrics = data['metrics']
        
        rows += f"    {name} & {config['type']} & {config['start_time']:02d}:00 & "
        rows += f"{metrics['tri_retention']:.1f} & {metrics['min_service']:.1f} & "
        rows += f"{metrics['r_tri_svc']:.3f} \\\\\n"
    
    latex = f"""
\\begin{{table}}[htbp]
\\centering
\\caption{{Summary of all disruption scenarios.}}
\\label{{tab:summary}}
\\begin{{tabular}}{{llcccc}}
\\toprule
\\textbf{{ID}} & \\textbf{{Type}} & \\textbf{{Start}} & \\textbf{{$|\\Delta_2|$ ret. (\\%)}} & \\textbf{{Svc$_{{\\min}}$ (\\%)}} & \\textbf{{$r$}} \\\\
\\midrule
{rows}\\bottomrule
\\end{{tabular}}
\\end{{table}}
"""
    
    if save_path:
        with open(save_path, 'w') as f:
            f.write(latex)
        print(f"  Saved: {save_path}")
    
    return latex
