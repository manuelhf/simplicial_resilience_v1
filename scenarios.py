"""
================================================================================
scenarios.py - Disruption Scenario Definition and Simulation
================================================================================

This module provides functions for defining and simulating disruption scenarios
in water distribution networks. Scenarios represent various failure modes that
test network resilience:

SCENARIO TYPES:
    - Pipe breaks (sudden closure)
    - Pump failures
    - Valve operations
    - Cascading failures (sequential breaks)
    - Combined stress (multiple simultaneous failures)

The module automatically identifies critical network components (trunk mains,
secondary mains) based on flow analysis, enabling hydraulically-justified
scenario selection.

Author: WDSA Research Team
Date: 2024
"""

import wntr
import networkx as nx
import pandas as pd
from scipy.stats import pearsonr

from .operational import track_operational_complex, compute_hydraulic_metrics


# =============================================================================
# NETWORK COMPONENT IDENTIFICATION
# =============================================================================

def identify_critical_components(water_network, simulation_results):
    """
    Identify critical network components based on hydraulic importance.
    
    This function ranks pipes by their maximum flow rate to identify:
    - Trunk mains: Top 10 pipes by flow (primary supply arteries)
    - Secondary mains: Rank 11-20 (distribution feeders)
    - Peripheral pipes: Lower-ranked pipes
    
    Parameters
    ----------
    water_network : wntr.network.WaterNetworkModel
        The WNTR water network model.
    simulation_results : wntr.sim.SimulationResults
        Results from a baseline hydraulic simulation.
    
    Returns
    -------
    components : dict
        Dictionary containing:
        - 'sources': List of reservoir and tank names
        - 'pumps': List of pump names
        - 'trunk_mains': Top 10 pipes by max flow
        - 'secondary_mains': Pipes ranked 11-20 by max flow
        - 'peripheral': Pipes ranked 31-50
        - 'pipe_ranking': Full DataFrame with pipe flow rankings
    
    Notes
    -----
    The flow-based ranking provides a hydraulically-meaningful way to
    prioritize components for failure analysis. Trunk mains carry the
    most flow and their failure typically has the largest impact.
    
    Example
    -------
    >>> results = sim.run_sim()
    >>> components = identify_critical_components(wn, results)
    >>> print(f"Trunk mains: {components['trunk_mains'][:3]}")
    """
    # Identify sources (reservoirs and tanks)
    sources = list(water_network.reservoir_name_list) + list(water_network.tank_name_list)
    
    # Identify pumps
    pumps = list(water_network.pump_name_list)
    
    # Rank pipes by maximum absolute flow
    flow_data = simulation_results.link['flowrate']
    pipe_flows = []
    
    for pipe_name, pipe in water_network.pipes():
        if pipe_name in flow_data.columns:
            max_flow = flow_data[pipe_name].abs().max()
            pipe_flows.append({
                'pipe': pipe_name,
                'max_flow': max_flow,
                'diameter': pipe.diameter
            })
    
    # Create ranking DataFrame
    pipe_df = pd.DataFrame(pipe_flows).sort_values('max_flow', ascending=False)
    pipe_df['rank'] = range(1, len(pipe_df) + 1)
    
    # Extract component groups
    trunk_mains = pipe_df.head(10)['pipe'].tolist()
    secondary_mains = pipe_df.iloc[10:20]['pipe'].tolist()
    peripheral = pipe_df.iloc[30:50]['pipe'].tolist()
    
    return {
        'sources': sources,
        'pumps': pumps,
        'trunk_mains': trunk_mains,
        'secondary_mains': secondary_mains,
        'peripheral': peripheral,
        'pipe_ranking': pipe_df
    }


# =============================================================================
# SCENARIO DEFINITIONS
# =============================================================================

def define_standard_scenarios(components):
    """
    Define a standard set of disruption scenarios for resilience testing.
    
    This function creates scenario configurations for:
    - S0: Baseline (no disruption)
    - S1a/b/c: Secondary main break at different demand levels
    - S2a/b/c: Trunk main break at different demand levels
    
    Parameters
    ----------
    components : dict
        Output from identify_critical_components().
    
    Returns
    -------
    scenarios : dict
        Dictionary of scenario configurations, each containing:
        - 'type': Scenario category (baseline, gradual, abrupt)
        - 'description': Human-readable description
        - 'start_time': Hour when disruption begins
        - 'duration_h': Simulation duration in hours
        - 'actions': List of disruption actions to apply
    
    Scenario Design Rationale
    -------------------------
    The three timing variants (02h, 08h, 14h) represent:
    - 02:00 - Minimum demand period (low stress)
    - 08:00 - Peak demand period (high stress)
    - 14:00 - Mid-day demand (moderate stress)
    
    This tests whether topological response patterns are consistent
    across different operational states.
    
    Example
    -------
    >>> components = identify_critical_components(wn, results)
    >>> scenarios = define_standard_scenarios(components)
    >>> print(list(scenarios.keys()))
    ['S0', 'S1a', 'S1b', 'S1c', 'S2a', 'S2b', 'S2c']
    """
    trunk = components['trunk_mains']
    secondary = components['secondary_mains']
    
    scenarios = {
        # Baseline: No disruption
        'S0': {
            'type': 'baseline',
            'description': 'Baseline (no failure)',
            'start_time': 0,
            'duration_h': 24,
            'actions': []
        },
        
        # S1: Gradual depressurisation (secondary main break)
        'S1a': {
            'type': 'gradual',
            'description': 'Secondary main break at 02h (minimum demand)',
            'start_time': 2,
            'duration_h': 24,
            'actions': [{'type': 'pipe_break', 'element': secondary[0], 'time_h': 2}]
        },
        'S1b': {
            'type': 'gradual',
            'description': 'Secondary main break at 08h (peak demand)',
            'start_time': 8,
            'duration_h': 24,
            'actions': [{'type': 'pipe_break', 'element': secondary[0], 'time_h': 8}]
        },
        'S1c': {
            'type': 'gradual',
            'description': 'Secondary main break at 14h (mid demand)',
            'start_time': 14,
            'duration_h': 24,
            'actions': [{'type': 'pipe_break', 'element': secondary[0], 'time_h': 14}]
        },
        
        # S2: Abrupt isolation (trunk main break)
        'S2a': {
            'type': 'abrupt',
            'description': 'Trunk main break at 02h (minimum demand)',
            'start_time': 2,
            'duration_h': 24,
            'actions': [{'type': 'pipe_break', 'element': trunk[0], 'time_h': 2}]
        },
        'S2b': {
            'type': 'abrupt',
            'description': 'Trunk main break at 08h (peak demand)',
            'start_time': 8,
            'duration_h': 24,
            'actions': [{'type': 'pipe_break', 'element': trunk[0], 'time_h': 8}]
        },
        'S2c': {
            'type': 'abrupt',
            'description': 'Trunk main break at 14h (mid demand)',
            'start_time': 14,
            'duration_h': 24,
            'actions': [{'type': 'pipe_break', 'element': trunk[0], 'time_h': 14}]
        },
    }
    
    return scenarios


def create_custom_scenario(name, description, start_time, duration_h, actions):
    """
    Create a custom scenario configuration.
    
    Parameters
    ----------
    name : str
        Scenario identifier.
    description : str
        Human-readable description.
    start_time : int
        Hour when first disruption begins.
    duration_h : int
        Total simulation duration in hours.
    actions : list of dict
        List of disruption actions. Each action is a dict with:
        - 'type': 'pipe_break', 'pump_off', or 'valve_close'
        - 'element': Name of the element to affect
        - 'time_h': Hour when action occurs
    
    Returns
    -------
    scenario : dict
        Scenario configuration dictionary.
    
    Example
    -------
    >>> scenario = create_custom_scenario(
    ...     name='cascading',
    ...     description='Sequential pipe failures',
    ...     start_time=8,
    ...     duration_h=24,
    ...     actions=[
    ...         {'type': 'pipe_break', 'element': 'pipe_1', 'time_h': 8},
    ...         {'type': 'pipe_break', 'element': 'pipe_2', 'time_h': 10},
    ...     ]
    ... )
    """
    return {
        'type': 'custom',
        'description': description,
        'start_time': start_time,
        'duration_h': duration_h,
        'actions': actions
    }


# =============================================================================
# SCENARIO EXECUTION
# =============================================================================

def run_scenario(inp_filename, scenario_name, scenario_config, 
                 network_graph, node_positions, pressure_threshold=20.0):
    """
    Execute a single disruption scenario and collect results.
    
    This function:
    1. Loads a fresh copy of the network
    2. Applies disruption actions as time-based controls
    3. Runs the hydraulic simulation
    4. Tracks topological and hydraulic metrics
    5. Computes summary statistics and correlations
    
    Parameters
    ----------
    inp_filename : str
        Path to the EPANET .inp file.
    scenario_name : str
        Identifier for this scenario.
    scenario_config : dict
        Configuration from define_standard_scenarios() or create_custom_scenario().
    network_graph : networkx.Graph
        Full network graph for complex construction.
    node_positions : dict
        Node coordinates {node_id: (x, y)}.
    pressure_threshold : float, optional
        Minimum pressure for service (default 20.0 m).
    
    Returns
    -------
    results : dict or None
        Dictionary containing:
        - 'timeline': Topological metrics over time (DataFrame)
        - 'hydraulic': Hydraulic metrics over time (DataFrame)
        - 'merged': Combined DataFrame for correlation
        - 'config': Original scenario configuration
        - 'metrics': Summary statistics dict
        
        Returns None if simulation fails.
    
    Applied Actions
    ---------------
    PIPE_BREAK: Sets pipe status to closed (0) at specified time.
    PUMP_OFF: Sets pump status to closed (0) at specified time.
    VALVE_CLOSE: Sets valve status to closed (0) at specified time.
    
    Example
    -------
    >>> results = run_scenario('network.inp', 'S1a', scenario_config, G, pos)
    >>> if results:
    ...     print(f"Correlation: {results['metrics']['r_tri_svc']:.3f}")
    """
    print(f"  Running {scenario_name}: {scenario_config['description'][:50]}...")
    
    # Load fresh network model
    wn = wntr.network.WaterNetworkModel(inp_filename)
    
    # Configure simulation timing
    wn.options.time.duration = scenario_config['duration_h'] * 3600
    wn.options.time.hydraulic_timestep = 300  # 5-minute intervals
    wn.options.time.report_timestep = 300
    
    # Apply disruption actions
    for i, action in enumerate(scenario_config['actions']):
        action_time_seconds = action['time_h'] * 3600
        element_name = action['element']
        
        try:
            if action['type'] == 'pipe_break':
                element = wn.get_link(element_name)
                control_action = wntr.network.controls.ControlAction(
                    element, 'status', 0  # 0 = closed
                )
            elif action['type'] == 'pump_off':
                element = wn.get_link(element_name)
                control_action = wntr.network.controls.ControlAction(
                    element, 'status', 0
                )
            elif action['type'] == 'valve_close':
                element = wn.get_link(element_name)
                control_action = wntr.network.controls.ControlAction(
                    element, 'status', 0
                )
            else:
                print(f"    Warning: Unknown action type '{action['type']}'")
                continue
            
            # Create time-based control
            condition = wntr.network.controls.SimTimeCondition(
                wn, '>=', action_time_seconds
            )
            control = wntr.network.controls.Control(condition, control_action)
            wn.add_control(f'{scenario_name}_action_{i}', control)
            
        except Exception as e:
            print(f"    Warning: Could not apply action {i}: {e}")
    
    # Run simulation
    try:
        sim = wntr.sim.EpanetSimulator(wn)
        sim_results = sim.run_sim()
    except Exception as e:
        print(f"    Simulation failed: {e}")
        return None
    
    # Track metrics
    timeline = track_operational_complex(
        wn, sim_results, network_graph, node_positions, pressure_threshold
    )
    hydraulic = compute_hydraulic_metrics(wn, sim_results, pressure_threshold)
    
    # Merge for correlation analysis
    merged = pd.merge(timeline, hydraulic, on='time_hours')
    
    # Compute correlation between topology and service
    if merged['n_triangles_mesh'].std() > 0:
        r_tri_svc, _ = pearsonr(merged['n_triangles_mesh'], merged['service_pct'])
    else:
        r_tri_svc = 0.0
    
    # Summary statistics
    tri_init = timeline['n_triangles_mesh'].iloc[0]
    tri_min = timeline['n_triangles_mesh'].min()
    tri_retention = (tri_min / tri_init * 100) if tri_init > 0 else 0
    min_service = hydraulic['service_pct'].min()
    
    print(f"    |Δ₂|: {tri_init} → {tri_min} ({tri_retention:.1f}%)")
    print(f"    Service min: {min_service:.1f}%, r = {r_tri_svc:.3f}")
    
    return {
        'timeline': timeline,
        'hydraulic': hydraulic,
        'merged': merged,
        'config': scenario_config,
        'metrics': {
            'tri_init': tri_init,
            'tri_min': tri_min,
            'tri_retention': tri_retention,
            'min_service': min_service,
            'r_tri_svc': r_tri_svc
        }
    }


def run_all_scenarios(inp_filename, scenarios, network_graph, 
                      node_positions, pressure_threshold=20.0):
    """
    Execute all scenarios and collect results.
    
    Parameters
    ----------
    inp_filename : str
        Path to the EPANET .inp file.
    scenarios : dict
        Dictionary of scenario configurations.
    network_graph : networkx.Graph
        Full network graph.
    node_positions : dict
        Node coordinates.
    pressure_threshold : float, optional
        Minimum pressure for service.
    
    Returns
    -------
    all_results : dict
        Dictionary mapping scenario names to their results.
    
    Example
    -------
    >>> scenarios = define_standard_scenarios(components)
    >>> all_results = run_all_scenarios('network.inp', scenarios, G, pos)
    >>> for name, data in all_results.items():
    ...     print(f"{name}: r = {data['metrics']['r_tri_svc']:.3f}")
    """
    all_results = {}
    
    for name, config in scenarios.items():
        result = run_scenario(
            inp_filename, name, config, 
            network_graph, node_positions, pressure_threshold
        )
        if result is not None:
            all_results[name] = result
    
    return all_results
