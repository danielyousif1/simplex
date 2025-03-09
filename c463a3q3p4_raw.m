%% CO 463/663 ASSIGNMENT 3
% We write an algorithm to find if a point $b$ is in the convex hull and/or
% conic hull of a sequence of $k$ points (stored in a matrix $S$). We avoid
% using MATLAB's built-in linprog function, and instead implement Phase-I
% of the Simplex algorithm where we encode the appropriate constraints for
% each problem. Because we are only looking for feasibility (ie. whether a
% convex or conical combination exists) there is nothing to optimize and as
% we know from Phase-I there is no objective function. As expected by
% Carathéodory's Theorem, our convex hull algorithm uses at most $n + 1$
% points and our conical hull membership algorithm uses at most $n$ points.
% For illustration, we use $n = 3, k = 4$.
%% Setup, Seeds
rng(34)

n = 3;
k = 4;

S = rand(n, k);
b = rand(n, 1);
%% 1. Convex Hull Membership -- Using Phase-1 Simplex
% We implement Phase 1 Simplex. 0 objective function with the constraints:
% (i) Sx = b, (ii) x >= 0, (iii) sum(x) = 1

A = [S; ones(1, k)];
b_ext = [b; 1];

% Enforce b >= 0
for i = 1:length(b_ext)
    if b_ext(i) < 0
        A(i, :) = -A(i, :);
        b_ext(i) = -b_ext(i);
    end
end

% Introduce artificial variables to obtain basic feasible solution
m = size(A, 1);
I_art = eye(m);
A_ext = [A, I_art];
c = [zeros(1, k), ones(1, m)];

% Basic Feasible Solution
B_idx = k+1 : k+m; % Basic Variables
N_idx = 1 : k; % Non Basic Variables
x_B = b_ext;

% Phase-1 Simplex iterations for convex hull:
maxIter = 1000;
feasible_convex = false;
for iter = 1:maxIter
    B = A_ext(:, B_idx);
    N = A_ext(:, N_idx);
    c_B = c(B_idx);
    c_N = c(N_idx);
    y = (c_B / B);
    r_N = c_N - y * N;
    % Check optimality (for Phase-1, optimality means no negative reduced cost, since we minimize sum of artificials).
    if all(r_N >= 0)
        current_obj = c_B * x_B;
        if abs(current_obj) < 1e-9 % For float
            feasible_convex = true;
        end
        break;
    end
    [~, enter_idx_rel] = min(r_N);
    entering = N_idx(enter_idx_rel);
    d = B \ A_ext(:, entering);
    if all(d <= 1e-12)
        disp('Phase-1 objective unbounded');
        break;
    end
    theta = inf;
    leave_idx_rel = -1;
    for i = 1:m
        if d(i) > 1e-12
            if x_B(i) / d(i) < theta
                theta = x_B(i) / d(i);
                leave_idx_rel = i;
            end
        end
    end
    leaving = B_idx(leave_idx_rel);
    B_idx(leave_idx_rel) = entering;
    N_idx(enter_idx_rel) = leaving;
    B = A_ext(:, B_idx);
    x_B = B \ b_ext;
end

% Check feasibility
if feasible_convex
    % Extract first k coefficients
    full_solution = zeros(k + m, 1);
    full_solution(B_idx) = x_B;
    x_conv_sol = full_solution(1:k);
    disp('Point b IS inside the convex hull of S.');
    disp('Convex combination coefficients:');
    disp(x_conv_sol');
else
    disp('b IS NOT in the convex hull of S.');
end

%% 2. Conical Hull Membership -- Using Phase-1 Simplex
% We implement Phase 1 Simplex. 0 objective function with the constraints:
% (i) Sx = b, (ii) x >= 0

A_cone = S;
b_cone = b;
m_cone = n;

% Enforce b >= 0
for i = 1:m_cone
    if b_cone(i) < 0
        A_cone(i, :) = -A_cone(i, :);
        b_cone(i) = -b_cone(i);
    end
end

% Introduce artificial variables to obtain basic feasible solution
I_art_cone = eye(m_cone);
A_ext_cone = [A_cone, I_art_cone];
c_cone = [zeros(1, k), ones(1, m_cone)];

% Basic Feasible Solution
B_idx = k+1 : k+m_cone; % Basic Variables
N_idx = 1 : k; % Non Basic Variables
x_B = b_cone;

% Phase-1 Simplex iterations for conic hull:
feasible_conical = false;
for iter = 1:maxIter
    B = A_ext_cone(:, B_idx);
    N = A_ext_cone(:, N_idx);
    c_B = c_cone(B_idx);
    c_N = c_cone(N_idx);
    y = (c_B / B);
    r_N = c_N - y * N;
    if all(r_N >= 0)
        if abs(c_B * x_B) < 1e-9 % For float
            feasible_conical = true;
        end
        break;
    end
    [~, enter_idx_rel] = min(r_N);
    entering = N_idx(enter_idx_rel);
    d = B \ A_ext_cone(:, entering);
    if all(d <= 1e-12)
        disp('Phase-1 unbounded');
        break;
    end
    theta = inf;
    leave_idx_rel = -1;
    for i = 1:length(x_B)
        if d(i) > 1e-12
            if x_B(i) / d(i) < theta
                theta = x_B(i) / d(i);
                leave_idx_rel = i;
            end
        end
    end
    leaving = B_idx(leave_idx_rel);
    B_idx(leave_idx_rel) = entering;
    N_idx(enter_idx_rel) = leaving;
    B = A_ext_cone(:, B_idx);
    x_B = B \ b_cone;
end

% Check feasibility
if feasible_conical
    full_solution = zeros(k + m_cone, 1);
    full_solution(B_idx) = x_B;
    x_conical_sol = full_solution(1:k);
    disp('Point b IS inside the conical hull of S.');
    disp('Conical combination coefficients:');
    disp(x_conical_sol');
else
    disp('Point b IS NOT in the conical hull of S.');
end