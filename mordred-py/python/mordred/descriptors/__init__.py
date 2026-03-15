"""All available molecular descriptors.

Pass this module to ``Calculator()`` to use all descriptors::

    calc = Calculator(descriptors, ignore_3D=True)

This module serves as a sentinel -- when passed to ``Calculator()``, it
signals "use all available descriptors." Individual descriptor classes will
be added as they are implemented.
"""
