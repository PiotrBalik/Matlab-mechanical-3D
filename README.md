# Matlab-mechanical-3D
Collaboration with a PhD researcher in structural engineering at Nottingham Trent University, UK. Project investigated the optimum design of three-dimensional steel frames in terms of real-world cost function (weight, embodied carbon). Imposed behavioral constraints were based on EC3 specifications.

Since the problem is _hard_ to formulate in [YALMIP](https://yalmip.github.io/), metaheuristics were considered. Fully functional optimization algorithms interfaced:
* [Genetic Algorithm](https://www.springer.com/gp/book/9783319521558) (Matlab built-in)
* [Particle Swarm Optimization](https://www.springer.com/gp/book/9783642378454) (Matlab built-in)
* [Ant Colony Optimization](https://www.springer.com/gp/book/9783030673796) (custom)
* [Multi-Objective GA](https://www.springer.com/gp/book/9789811314704) (Matlab built-in)
* [Multi-Objective PSO](https://www.springer.com/gp/book/9783642051647) (custom)

Elements of speeding up the code (precomputing, inlining) and metaprogramming - dynamic creation in _buildGroup_. Possible extension of program editing it's code and compiling it - for hopefully faster execution instead of _evaluating_ the anonymous function.
