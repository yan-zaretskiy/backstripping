from itertools import accumulate, chain
import matplotlib.pyplot as plt
from matplotlib.pyplot import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

from traits.api import HasStrictTraits, Float

from layer import Layer, athy_porosity
from events import Deposition, EventManager


class PhysicalConstants(HasStrictTraits):
    gravity = Float(9.81)
    water_density = Float(1000)
    mantle_density = Float(3300)


def prepare_events(ages, bathymetries, sea_levels, thicknesses, rock_types, rock_properties):
    """ Package inputs into an EventManager instance. """
    events = []
    for age, bathymetry, sea_level, thickness, rock_type in zip(
          ages[1:], bathymetries[1:], sea_levels[1:], thicknesses, rock_types):
        layer = Layer(present_thickness=thickness, porosity_function=athy_porosity)
        layer.set_rock_properties(rock_properties[rock_type])
        event = Deposition(age=age, bathymetry=bathymetry, sea_level=sea_level, layer=layer)
        events.append(event)

    event_manager = EventManager(initial_age=ages[0],
                                 initial_sea_level=sea_levels[0],
                                 initial_bathymetry=bathymetries[0])
    event_manager.add_events(events)
    event_manager.reconstruct_burial_history()
    return event_manager


def compute_deflection(sediment_weight, sea_level_change, constants):
    """ helper function for Airy isostasy. """
    total_weight = sediment_weight + constants.gravity * constants.water_density * sea_level_change
    return total_weight / (constants.gravity * (constants.mantle_density - constants.water_density))


def compute_subsidence(event_manager, constants=PhysicalConstants()):
    """ Actual backstripping is performed here. """
    subsidence = []
    thickness_evolution = []
    for event_id in range(len(event_manager.events)):
        thickness, weight = event_manager.decompact_layers(event_id, constants)
        total_thickness = sum(thickness)
        total_weight = sum(weight)
        sea_level_change = event_manager.sea_level_change(event_id)
        bathymetry = event_manager.bathymetry(event_id)
        deflection = compute_deflection(total_weight, sea_level_change, constants)
        s = (bathymetry + total_thickness - deflection - sea_level_change
             - event_manager.initial_bathymetry)
        subsidence.append(s)
        thickness_evolution.append(thickness)
    return subsidence[::-1], thickness_evolution[::-1]


def plot_results(ages, subsidence, thickness_list, sea_levels, bathymetries):
    """ Plots tectonic subsidence and sediment thickness change over time. """
    # plot setup
    fig = plt.figure(figsize=(12, 9), facecolor='white')
    axes = plt.gca()
    plt.grid()
    axes.invert_xaxis()
    axes.set_xlabel("Time [Ma]", labelpad=15)
    axes.set_ylabel("Depth [m]", labelpad=15)
    axes.tick_params(axis='both', which='major', pad=10, direction='out', size=5)

    axes.xaxis.label.set_fontsize(18)
    axes.yaxis.label.set_fontsize(18)
    for item in axes.get_xticklabels() + axes.get_yticklabels():
        item.set_fontsize(14)

    # actual data
    subs, = plt.plot(ages, [0] + subsidence, '--', color='#0077B8', lw=5, label='Subsidence')
    plt.legend(handles=[subs], loc=3)
    horizon_offset = [w - bathymetries[0] - (s - sea_levels[0])
                      for w, s in zip(bathymetries, sea_levels)]
    horizons = [list(accumulate(chain([ho], t))) for ho, t in zip(horizon_offset, thickness_list)]

    axes.set_xlim([max(ages), min(ages)])
    axes.set_ylim([max(max(h) for h in horizons), min(horizon_offset)])

    patches = []
    n_patches = len(ages) - 1
    x_indices = [n_patches]
    for i, j in enumerate(reversed(range(n_patches))):
        x_indices = [j] + x_indices + [j+1]
        y_indices = list(range(i+1)) + list(range(i+1, -1, -1))
        points = [[ages[x_id], horizons[x_id-1][y_id]] for x_id, y_id in zip(x_indices, y_indices)]
        patches.append(Polygon(points))
    # we have to fix the last polygon
    xy = patches[-1].get_xy()
    xy[0, 1] = xy[-1, 1] = 0.0
    patches[-1].set_xy(xy)

    p = PatchCollection(patches, cmap=discrete_cmap(n_patches, 'terrain'), alpha=0.7)
    p.set_array(np.arange(n_patches))
    axes.add_collection(p)
    plt.show()


def discrete_cmap(N, base_cmap=None):
    """ Create an N-bin discrete colormap from the specified input map. """
    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return LinearSegmentedColormap.from_list(cmap_name, color_list, N)
