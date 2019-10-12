"""
.. module:: freeform
    :platform: Unix, Windows
    :synopsis: Provides freeform geometry classes

.. moduleauthor:: Onur R. Bingol <contact@onurbingol.net>

"""

from geomdl import abstract


class Freeform(abstract.Geometry):
    """ n-dimensional freeform geometry """
    def __init__(self, *args, **kwargs):
        super(Freeform, self).__init__(*args, **kwargs)
        self._idt['geometry_type'] = "freeform"
        self.name = "freeform geometry"

    @property
    def data(self):
        """ Returns a dict which contains the geometry data.

        Please refer to the `wiki <https://github.com/orbingol/NURBS-Python/wiki/Using-Python-Properties>`_ for details
        on using this class member.
        """
        return dict(
            type=self.type,
            points=tuple(self.evalpts)
        )

    def evaluate(self, **kwargs):
        """ Sets points that form the geometry.

        Keyword Arguments:
            * ``points``: sets the points
        """
        self._eval_points = kwargs.get('points', self._init_array())
        self._idt['dimension'] = len(self._eval_points[0])
