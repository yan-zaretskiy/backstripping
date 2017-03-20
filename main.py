from backstripping import prepare_events, compute_subsidence, plot_results

# setup layers
rock_properties = [
    {'surface_porosity': 0.63, 'compaction_rate': 0.51e-3, 'sediment_density': 2720.0},
    {'surface_porosity': 0.49, 'compaction_rate': 0.27e-3, 'sediment_density': 2650.0},
    {'surface_porosity': 0.70, 'compaction_rate': 0.71e-3, 'sediment_density': 2710.0},
    {'surface_porosity': 0.40, 'compaction_rate': 0.60e-3, 'sediment_density': 2720.0},
    {'surface_porosity': 0.20, 'compaction_rate': 0.60e-3, 'sediment_density': 2870.0},
    {'surface_porosity': 0.05, 'compaction_rate': 0.20e-3, 'sediment_density': 2960.0},
]

ages = [260, 245, 210, 160, 145, 125, 100, 80, 55, 45, 0]
sea_levels = [10, 0, 0, -20, -40, 70, 80, 100, 50, 40, 0]
bathymetries = [-20, 0, 20, 10, 20, 20, 200, 300, 350, 325, 300]
rock_types = [4, 5, 1, 4, 3, 1, 2, 0, 1, 0]
thicknesses = [400, 750, 250, 400, 200, 900, 1300, 750, 250, 200]

event_manager = prepare_events(ages, bathymetries, sea_levels, thicknesses, rock_types, rock_properties)
subsidence, thickness_evolution = compute_subsidence(event_manager)
plot_results(ages, subsidence, thickness_evolution, sea_levels, bathymetries)
