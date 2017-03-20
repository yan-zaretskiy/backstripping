from sortedcontainers import SortedListWithKey
from traits.api import HasStrictTraits, Float, Instance

from layer import Layer


class Deposition(HasStrictTraits):
    age = Float
    bathymetry = Float
    sea_level = Float
    layer = Instance(Layer)


class EventManager(HasStrictTraits):
    events = Instance(SortedListWithKey, kw={'key': lambda e: e.age})
    initial_age = Float
    initial_sea_level = Float
    initial_bathymetry = Float

    def add_events(self, events):
        self.events.update(events)

    def reconstruct_burial_history(self):
        """ Compute maximum burial depths for all the deposited layers. """
        current_burial = 0.0
        for event in self.events:
            event.layer.maximum_burial = current_burial
            current_burial += event.layer.present_thickness

    def decompact_layers(self, starting_event_id, constants):
        """ Decompaction of a sediment column. """
        current_burial = 0.0
        thickness_list = []
        weight_list = []
        for event in self.events[starting_event_id:]:
            thickness = event.layer.thickness_at_depth(current_burial)
            weight = event.layer.sediment_weight(constants)
            current_burial += thickness
            thickness_list.append(thickness)
            weight_list.append(weight)

        return thickness_list, weight_list

    def sea_level_change(self, event_id):
        """ Sea level change for a given event ID. """
        return self.events[event_id].sea_level - self.initial_sea_level

    def bathymetry(self, event_id):
        """ Water depth value for a given event ID. """
        return self.events[event_id].bathymetry
