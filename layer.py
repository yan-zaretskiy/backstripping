import math

from traits.api import (HasStrictTraits, Function, Float,
                        Property, cached_property)


class Layer(HasStrictTraits):
    maximum_burial = Float
    present_thickness = Float
    porosity_function = Function
    compaction_rate = Float
    surface_porosity = Float
    sediment_density = Float
    sediment_thickness = Property(depends_on=['present_thickness,'
                                              'maximum_burial'])

    def set_rock_properties(self, properties_dict):
        self.surface_porosity = properties_dict['surface_porosity']
        self.compaction_rate = properties_dict['compaction_rate']
        self.sediment_density = properties_dict['sediment_density']

    @cached_property
    def _get_sediment_thickness(self):
        """ Compute sediment thickness. """
        z0 = self.maximum_burial
        z1 = z0 + self.present_thickness
        water_thickness = self.integrate_porosity_function(z0, z1)
        return self.present_thickness - water_thickness

    def integrate_porosity_function(self, z0, z1):
        """ Numerically integrate porosity function over the given interval. """
        w = 0.5773502691896257  # sqrt(3)/3
        halflength = 0.5 * (z1 - z0)
        midpoint = 0.5 * (z0 + z1)

        porosity_0 = self.porosity_function(self, midpoint + halflength * w)
        porosity_1 = self.porosity_function(self, midpoint - halflength * w)
        return halflength * (porosity_0 + porosity_1)

    def thickness_at_depth(self, depth, eps=1e-6):
        """ Computes layer's thickness if buried at a given depth. """
        thickness = self.present_thickness  # initial guess
        # Newton iteration
        carry_on = True
        while carry_on:
            water_thickness = self.integrate_porosity_function(depth, depth + thickness)
            function_value = thickness - self.sediment_thickness - water_thickness
            derivative_value = 1.0 - self.porosity_function(self, depth + thickness)
            thickness -= function_value / derivative_value
            carry_on = abs(function_value) > eps
        return thickness

    def sediment_weight(self, constants):
        """ Layer weight above that of water. """
        return (self.sediment_density - constants.water_density) \
            * constants.gravity * self.sediment_thickness


def athy_porosity(layer, z):
    """ Athy's porosity-depth relationship. """
    return layer.surface_porosity * math.exp(-layer.compaction_rate * z)
