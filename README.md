# Matlab-mechanical-3D
Collaboration with a PhD researcher in structural engineering at Nottingham Trent University, UK. Project investigated the optimum design of three-dimensional steel frames in terms of real-world cost function (weight, embodied carbon). Imposed behavioral constraints were based on EC3 specifications.

Fully functionnal optimization algorithms interfaced:
* Genetic Algorithm (Matlab built-in)
* Particle Swarm Optimization (Matlab built-in)
* Ant Colony Optimization (custom)
* Multi-Objective GA (Matlab built-in)
* Multi-Objective PSO (custom)

Elements of speeding up the code (precomputing, inlining) and elements of metaprogramming - dynamic creation in _buildGroup_. Possible extension of program editing it's code and compiling it - for hopefully faster execution instead of _evaluating_ the anonymous function.
