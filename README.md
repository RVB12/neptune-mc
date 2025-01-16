# EIRON2

So this package is just me (Thijs) playing around with some monte carlo stuff. I've called it eiron2, because there is some hope that if it goes well enough, it could become something of an expansion of the capabilities of eiron. But that is a long way off and not currently an actual goal.

### Running the example

To run the example, you'll need:

- Python
- CMake
- C++ compilers
- pybind11

If you have everything installed, running the example is as simple as:

- `pip install .`
- `python3 examples/pythontest_2d.py`

### General design

I really want to have a strong separation between the geometry and the physics. The idea is that the geometry should only contain geometry things, like whether a point is inside a cell, the size of a cell, information about where boundaries are, when a particle with given position and direction will change cells, ... The physics side should be able to do all the physics without direction doing any geometry.

A downside of this generality is that there can be no optimization for special cases. For example, the boundary conditions use the surface normal for reflections, but this could be simplified for square domains.

