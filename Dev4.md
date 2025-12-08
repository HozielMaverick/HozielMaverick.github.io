```matlab
function Dev4_472_MH
% ECSE 472 – Project Deliverable 4
% Harmonic Balance (HB) using Newton–Raphson, Fourier series and Γ matrix.
%
% Circuit: single node with C || R || diode driven by:
%   i_s(t) = I0 * cos(2*pi*f0*t)
%
% KCL in time domain:
%   i_R(t) + i_C(t) + i_D(t) = i_s(t)
%
% For HB we:
%   - Represent v(t) by a truncated Fourier series with harmonics -H..H
%   - Use a Γ matrix to move between time samples and Fourier coefficients
%   - Apply the linear admittance Y(jkω0) in the harmonic domain
%   - Transform back to time to build:
%         Phi(x) = G_HB * x + f(x) - b_HB  =  0
%   - Solve with Newton using the pattern:
%         1) f(x), 2) f'(x), 3) Phi, 4) J, 5) Δx, 6) x:=x+Δx, 7) ||Δx||_inf

    clc;
    warning('off', 'MATLAB:nearlySingularMatrix');
    warning('off', 'MATLAB:singularMatrix');
    warning('off', 'MATLAB:rankDeficientMatrix');


    %------------------------
    % Circuit parameters 
    %------------------------
    f0  = 1e3;          % [Hz] fundamental
    w0  = 2*pi*f0; %omega = 2pif
    T   = 1/f0;

    I0  = 10e-3;        % [A] current source amplitude
    R   = 1e3;          % [Ω]
    C   = 1e-6;         % [F]

    Is  = 1e-14;        % [A] diode saturation current
    Vt  = 0.025;        % [V] thermal voltage

    %------------------------
    % HB setup: harmonics and sampling (Γ)
    %------------------------
    H = 10;                      % number of harmonics on each side (−H..H)
    N = 2*H + 1;                % number of time samples = number of coeffs

    % Time samples over one period:
    % t_k = k*T/N,  k = 0..N-1
    t_HB = (0:N-1).' * (T/N); %col vector of T/N, one full period sampled at N points

    % Γ matrix and its inverse (Γ⁻¹ = (1/N)*Γᴴ for this uniform sampling)
    [Gamma, invGamma, Yvec] = build_gamma_and_Y(R, C, f0, H, t_HB);

    % Linear HB operator:
    %   i_lin(t) = Γ * diag(Y) * Γ⁻¹ * v(t)
    G_HB = Gamma * diag(Yvec) * invGamma;
    G_HB = real(G_HB); % should be real; strip tiny imag noise

    % Current source in time at the HB sampling instants
    i_src_HB = I0 * cos(w0 * t_HB);    % b_HB

    %------------------------
    % Nonlinear diode term (time domain)
    %------------------------
    f_nl  = @(v) diode_current_vector(v, Is, Vt);   % f(x)
    df_nl = @(v) diode_jacobian_diag(v, Is, Vt);    % f'(x) (diagonal)

    %------------------------
    % Initial guess for HB solution (one period)
    %------------------------
    v0_HB = zeros(N,1);

    %------------------------
    % Newton solve for Harmonic Balance (using Γ/Fourier-based G_HB)
    %------------------------
    [v_HB, converged, iters] = newton_NR_BE(G_HB, i_src_HB, v0_HB, f_nl, df_nl);

    if ~converged
        warning('HB Newton did NOT converge in %d iterations.', iters);
    else
        fprintf('HB Newton converged in %d iterations.\n', iters);
    end

    %------------------------
    % Time-domain Backward Euler simulation (Deliverable 3 style),
    % used here only to verify the HB steady-state.
    %------------------------
    nPeriods_BE = 10;
    [t_BE, v_BE] = simulate_BE_time_domain(f0, I0, R, C, Is, Vt, nPeriods_BE);

    %------------------------
    % Plots
    %------------------------
    figure;

    % For plotting, append v(0) at t = T so the curve clearly shows full period
    t_plot = [t_HB; T];
    v_plot = [v_HB; v_HB(1)];

    subplot(2,1,1);
    plot(t_plot*1e3, v_plot, 'LineWidth', 1.5);
    xlabel('t [ms]');
    ylabel('v_{HB}(t) [V]');
    title('Harmonic Balance steady-state (one period)');
    grid on;

    subplot(2,1,2);
    % show last period of BE
    T = 1/f0;
    t_start = t_BE(end) - T;
    idx = t_BE >= t_start;
    plot((t_BE(idx)-t_start)*1e3, v_BE(idx), 'LineWidth', 1.5);
    xlabel('t [ms]');
    ylabel('v_{BE}(t) [V]');
    title('Backward Euler steady-state (last period)');
    grid on;

    sgtitle('Deliverable 4 – HB vs BE');

end



%====================================================================
% Local functions
%====================================================================

function [Gamma, invGamma, Yvec] = build_gamma_and_Y(R, C, f0, H, t)
% Build Γ and Γ⁻¹ matrices and the harmonic admittances Y(jkω0).
%
%   Harmonic indices: n = -H, ..., H   (total N = 2H+1)
%   Time samples    : t_k = given in vector t (N x 1), one period.
%
% Γ(k,m) = exp(j * n(m) * ω0 * t_k)
% Γ⁻¹    = (1/N) * Γᴴ    (for this uniform sampling)
%
% Linear admittance at each harmonic:
%   for n = 0:  Y = 1/R
%   for n ≠ 0:  Y = 1/R + j * n * ω0 * C

    w0 = 2*pi*f0;
    N  = numel(t);
    n_vec = (-H:H).'; % column vector of harmonic indices (N x 1)

    % Build Γ
    Gamma = zeros(N, N);
    for k = 1:N %time samples
        Gamma(k,:) = exp(1j * n_vec.' * w0 * t(k)); %row vector e^(j(-H)wotk)..e^(j(+H)wotk)
    end

    % Inverse Γ for this DFT-like matrix
    invGamma = (1/N) * Gamma'; %invgamma is 1/N * conjugate transpose of gamma

    % Harmonic admittances Y(j n ω0)
    Yvec = zeros(N,1); %1 col vector
    for m = 1:N
        n = n_vec(m); %fill n_vec with time samples
        if n == 0
            Yvec(m) = 1/R; % DC: capacitor open
        else
            Yvec(m) = 1/R + 1j*n*w0*C;  % R || C admittance at harmonic n
        end
    end
end


function i_d = diode_current_vector(v, Is, Vt)
% Vector of diode currents in time:
%   I_d(v) = Is (exp(v/Vt) - 1)

    i_d = Is * (exp(v./Vt) - 1);
end


function J_diag = diode_jacobian_diag(v, Is, Vt)
% Diagonal entries of f'(x) for the diode nonlinearity:
%   d/dv [Is (exp(v/Vt) - 1)] = (Is/Vt) exp(v/Vt)

    J_diag = (Is / Vt) * exp(v./Vt);
end


function [x, converged, it] = newton_NR_BE(G_BE, b_BE, x0, f_handle, df_handle)
% Generic Newton–Raphson for:
%   Phi(x) = G_BE * x + f(x) - b_BE = 0
%
% Pattern (as required in the assignment):
%   1) f(x)
%   2) f'(x)
%   3) Phi(x) = G_BE * x + f(x) - b_BE
%   4) J(x)   = G_BE + f'(x)
%   5) Δx     = - J \ Phi
%   6) x      := x + Δx
%   7) converge if ||Δx||_inf <= 1e-6 (max 100 iterations)

    max_iters = 100;
    tol       = 1e-6;

    x = x0;
    converged = false;

    for it = 1:max_iters
        % 1) f(x)
        f_val = f_handle(x);

        % 2) f'(x)
        df_val = df_handle(x);

        % if df_handle returns a vector, treat as diagonal Jacobian
        if isvector(df_val)
            J_nl = spdiags(df_val, 0, numel(df_val), numel(df_val));
        else
            J_nl = df_val;
        end

        % 3) Phi(x)
        Phi = G_BE * x + f_val - b_BE;

        % 4) J(x)
        J = G_BE + J_nl;

        % 5) Δx = -J \ Phi
        dx = - J \ Phi;

        % 6) update x
        x = x + dx;

        % 7) convergence
        if norm(dx, inf) <= tol
            converged = true;
            break;
        end
    end
end


function [t_vec, v_vec] = simulate_BE_time_domain(f0, I0, R, C, Is, Vt, nPeriods)
% Nonlinear Backward Euler transient simulation
%   C dv/dt + v/R + I_d(v) = I0 cos(2*pi*f0 t)

    w0 = 2*pi*f0;
    T  = 1/f0;

    % Time step (you can change if you want)
    N_per_period = 128;
    dt = T / N_per_period;

    t_end = nPeriods * T;
    t_vec = (0:dt:t_end).';
    Nt    = length(t_vec);

    v_vec = zeros(Nt,1);
    v_prev = 0;     % initial condition

    for n = 2:Nt
        t_n   = t_vec(n);
        i_src = I0 * cos(w0 * t_n);

        % Scalar BE equation:
        %   C (v_n - v_prev)/dt + v_n/R + I_d(v_n) = i_src
        %
        % => G_BE * v_n + f(v_n) = b_BE
        % with
        %   G_BE = C/dt + 1/R
        %   f(v) = I_d(v)
        %   b_BE = C/dt * v_prev + i_src

        G_BE_scalar = C/dt + 1/R;
        b_BE_scalar = (C/dt)*v_prev + i_src;

        G_BE = G_BE_scalar;
        b_BE = b_BE_scalar;

        f_nl  = @(v) Is*(exp(v./Vt) - 1);
        df_nl = @(v) (Is/Vt)*exp(v./Vt);

        % Newton with previous value as initial guess
        [v_new, ~, ~] = newton_NR_BE(G_BE, b_BE, v_prev, f_nl, df_nl);

        v_prev   = v_new;
        v_vec(n) = v_new;
    end
end
```
