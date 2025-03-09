# Hull Membership
A feasibility solver that checks if a point $b$ in Euclidean space can be represented as a convex or conic combination (ie. $b \in \text{conv}S, b \in \text{cone}S$ . Implements Phase-I Simplex to find such coefficients.
In $\mathbb{R}^n$, he number of coefficients needed is at most $n+1$ and $n$, respectively. If interested, see Caratheodory's theorem (https://en.wikipedia.org/wiki/Carath%C3%A9odory%27s_theorem_(convex_hull)). For educational purposes, this was done "by scratch" and without using MATLAB's built-in `linprog` function.


To-DO:
1. Implement in C++ or Rust
2. Write test suite
