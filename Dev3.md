```matlab
function [t, Xhist, meta] = Dev3_472_MH(spice_file)
format shortE
% DEV3_472_MH — Deliverable 3: Transient (linear + nonlinear) with Backward Euler
%
% We reuse the same MNA structure and diode Newton-Raphson style
% as in Dev2_472_MH, but now:
%   - Parse .tran Tstop dt
%   - Parse COS sources (Vin n1 0 COS A F)
%   - Do Backward Euler: (C/h)(x^k - x^{k-1}) + G x^k + f(x^k) = b(t_k)
%     => (G + C/h) x^k + f(x^k) = b(t_k) + (C/h) x^{k-1}
%   - If there are diodes, solve each time step with Newton in the
%     same style as Deliverable 2 (Phi, J, dx, etc.).
%

%% ========= PASS 1: read / strip / collect (with .tran) =========
lines_raw = read_text_lines(spice_file);
[lines, has_op] = strip_comments_end_and_op(lines_raw); %#ok<NASGU>
[elems, node_map, counts, tran] = pass_collect_tran(lines);

if ~tran.enabled
    error('No .tran line found in netlist. Need ".tran Tstop dt".');
end

% Dimensions 
N   = counts.N;
nV  = counts.nV;   % independent V + 0-ohm R + L-as-short in .op
nL  = counts.nL;   % time-domain inductors
nE  = counts.nE;
nH  = counts.nH;
X   = N + nV + nL + nE + nH;

% Indices
idx.v  = 1:N;
idx.iV = N + (1:nV);
idx.iL = N + nV + (1:nL);
idx.iE = N + nV + nL + (1:nE);
idx.iH = N + nV + nL + nE + (1:nH);

% Allocate
G = zeros(X,X);
C = zeros(X,X);
b_dc = zeros(X,1);       % DC/time-independent part of b

% Branch incidence helpers
Bv = zeros(N,nV);  % independent V + 0-ohm R + L@.op
Bl = zeros(N,nL);  % inductors (time-domain)
Be = zeros(N,nE);  % VCVS
Bh = zeros(N,nH);  % CCVS

% Map "element name" -> absolute column index in x for branch-current unknowns
name2branch = containers.Map('KeyType','char','ValueType','int32');

%% ========= Pre-assign branch columns (same as Dev2) =========
vcolcount = 0; ecolcount = 0; hcolcount = 0; lcol_seen = 0;
for k = 1:numel(elems)
    e = elems{k}; T = upper(e.type);
    switch T
        case 'V'
            vcolcount = vcolcount + 1; e.branch_col = vcolcount;
            name2branch(upper(e.name)) = idx.iV(vcolcount);

        case 'R'
            if isfield(e,'is_short') && e.is_short
                vcolcount = vcolcount + 1; e.branch_col = vcolcount;
                name2branch(upper(e.name)) = idx.iV(vcolcount);
            else
                e.branch_col = 0;
            end

        case 'L'
            if isfield(e,'dc_short') && e.dc_short
                vcolcount = vcolcount + 1; e.branch_col = vcolcount;
                name2branch(upper(e.name)) = idx.iV(vcolcount);
            else
                lcol_seen = lcol_seen + 1; e.branch_col = lcol_seen;
            end

        case 'E'
            ecolcount = ecolcount + 1; e.branch_col = ecolcount;
            name2branch(upper(e.name)) = idx.iE(ecolcount);

        case 'H'
            hcolcount = hcolcount + 1; e.branch_col = hcolcount;
            name2branch(upper(e.name)) = idx.iH(hcolcount);

        otherwise
            e.branch_col = 0; % D, C, I, G, F
    end
    elems{k} = e;
end

%% ========= PASS 2: stamp linear devices =========
lcol_seen = 0;
for k = 1:numel(elems)
    e = elems{k}; T = upper(e.type);
    switch T
        case 'R'
            [n1,n2] = nn(e.np, e.nn, node_map);
            if isfield(e,'is_short') && e.is_short
                col = e.branch_col; Bv = stamp_branch(Bv, n1, n2, col);
            else
                g = 1/numval(e.value);
                G = stamp_passive(G, n1, n2, g);
            end

        case 'C'
            [n1,n2] = nn(e.np, e.nn, node_map);
            Cval = numval(e.value);
            C(idx.v, idx.v) = stamp_passive(C(idx.v,idx.v), n1, n2, Cval);

        case 'L'
            [n1,n2] = nn(e.np, e.nn, node_map);
            if isfield(e,'dc_short') && e.dc_short
                % For .op, L is a short; for transient we actually
                % treat it as inductor: here dc_short should be false.
                col = e.branch_col; Bv = stamp_branch(Bv, n1, n2, col);
            else
                lcol_seen = lcol_seen + 1;
                Bl = stamp_branch(Bl, n1, n2, lcol_seen);
                Lval = numval(e.value);
                il = idx.iL(lcol_seen);
                C(il, il) = C(il, il) + Lval; % L * diL/dt
            end

        case 'I'
            % SPICE: positive current flows from n+ to n-
            [nplus,nminus] = nn(e.np, e.nn, node_map);
            Ival = source_value(e);
            if nplus>0,  b_dc(nplus)  = b_dc(nplus)  - Ival; end
            if nminus>0, b_dc(nminus) = b_dc(nminus) + Ival; end

        case 'V'
            [n1,n2] = nn(e.np, e.nn, node_map);
            col = e.branch_col; Bv = stamp_branch(Bv, n1, n2, col);

            % DC or AC sources contribute to b_dc here.
            % COS sources are time-varying: we *skip* them here and
            % handle them inside the transient loop.
            if isfield(e,'src_kind') && strcmpi(e.src_kind,'COS')
                % nothing here; time-dependent part later
            else
                Vval = source_value(e);
                b_dc(idx.iV(col)) = b_dc(idx.iV(col)) + Vval;
            end

        case 'E'
            [n1,n2]   = nn(e.np,  e.nn,  node_map);
            [nc1,nc2] = nn(e.nc1, e.nc2, node_map);
            col = e.branch_col; Be = stamp_branch(Be, n1, n2, col);
            A = numval(e.value); ie = idx.iE(col);
            if nc1>0, G(ie, nc1) = G(ie, nc1) - A; end
            if nc2>0, G(ie, nc2) = G(ie, nc2) + A; end

        case 'G'
            [n1,n2]   = nn(e.np,  e.nn,  node_map);
            [nc1,nc2] = nn(e.nc1, e.nc2, node_map);
            gm = numval(e.value);
            if n1>0 && nc1>0, G(n1,nc1) = G(n1,nc1) + gm; end
            if n1>0 && nc2>0, G(n1,nc2) = G(n1,nc2) - gm; end
            if n2>0 && nc1>0, G(n2,nc1) = G(n2,nc1) - gm; end
            if n2>0 && nc2>0, G(n2,n2) = G(n2,n2) + gm; end

        case 'F'
            [n1,n2] = nn(e.np, e.nn, node_map);
            beta = numval(e.value);
            ctrl = upper(char(e.ctrl));
            assert(isKey(name2branch, ctrl), ...
                'CCCS %s controls unknown branch "%s".', e.name, ctrl);
            ictrl = name2branch(ctrl);
            if n1>0, G(n1, ictrl) = G(n1, ictrl) + beta; end
            if n2>0, G(n2, ictrl) = G(n2, ictrl) - beta; end

        case 'H'
            [n1,n2] = nn(e.np, e.nn, node_map);
            col = e.branch_col; Bh = stamp_branch(Bh, n1, n2, col);
            Rm = numval(e.value); ih = idx.iH(col);
            ctrl = upper(char(e.ctrl));
            assert(isKey(name2branch, ctrl), ...
                'CCVS %s controls unknown branch "%s".', e.name, ctrl);
            ictrl = name2branch(ctrl);
            G(ih, ictrl) = G(ih, ictrl) - Rm;

        case 'D'
            % handled via f(x), f'(x) in Newton

        otherwise
            % .tran already stripped inside pass_collect_tran
            error('Unsupported element type "%s" on line %d (%s).', ...
                  T, e.line_no, e.raw);
    end
end

% Couple branch equations into G
if nV>0
    G(idx.v,  idx.iV) = G(idx.v,  idx.iV) + Bv;
    G(idx.iV, idx.v ) = G(idx.iV, idx.v ) + Bv.';
end
if nL>0
    G(idx.v,  idx.iL) = G(idx.v,  idx.iL) + Bl;
    G(idx.iL, idx.v ) = G(idx.iL, idx.v ) + Bl.';
end
if nE>0
    G(idx.v,  idx.iE) = G(idx.v,  idx.iE) + Be;
    G(idx.iE, idx.v ) = G(idx.iE, idx.v ) + Be.';
end
if nH>0
    G(idx.v,  idx.iH) = G(idx.v,  idx.iH) + Bh;
    G(idx.iH, idx.v ) = G(idx.iH, idx.v ) + Bh.';
end

%% ========= Labels & meta =========
var_names = strings(1,X);
inv_node = invert_map(node_map, N);
for i = 1:N,  var_names(i) = "v(" + inv_node(i) + ")"; end
for i = 1:nV, var_names(N+i) = "i(V" + i + ")"; end
for i = 1:nL, var_names(N+nV+i) = "i(L" + i + ")"; end
for i = 1:nE, var_names(N+nV+nL+i) = "i(E" + i + ")"; end
for i = 1:nH, var_names(N+nV+nL+nE+i) = "i(H" + i + ")"; end

meta.node_map  = node_map;
meta.var_names = var_names;
meta.idx       = idx;
meta.counts    = counts;
meta.elems     = elems;
meta.tran      = tran;

% Summary
fprintf('\n=== MNA Summary for "%s" ===\n', spice_file);
fprintf('Nodes (non-ground): %d\n', N);
fprintf('Voltage sources:   %d\n', nV);
fprintf('Inductors:         %d\n', nL);
fprintf('VCVS:              %d\n', nE);
fprintf('CCVS (H):          %d\n', nH);
fprintf('Diodes:            %d\n', counts.nD);
disp(' '); disp('Variable ordering (x):'); disp(meta.var_names.');
disp(' '); disp('G matrix :'); disp(G);
disp('C matrix:'); disp(C);
disp('b_dc    :'); disp(b_dc);

%% ========= Transient Backward Euler solve =========
%computes time domain responses
nonlinear = counts.nD > 0;
h   = tran.dt;
t_end = tran.tstop;
nsteps = floor(t_end/h); %floor matlab function
t = (0:nsteps).' * h;   % column vector, time from 0 to t end

Xdim  = X;
Xhist = zeros(Xdim, nsteps+1); %nsteps +1 since IC + all steps
x_prev = zeros(Xdim,1);   % initial condition x(0) = 0 


G_BE = G + (C / h); %circuit equation LHS

for k = 1:nsteps
    tk = t(k+1);  % current time t_k

    % Build b(t_k) = b_dc + contributions from time varying sources cos
    b_t = b_dc;
    for eei = 1:numel(elems)
        ee = elems{eei};
        if upper(ee.type) == 'V'
            if isfield(ee, 'src_kind') && strcmpi(ee.src_kind,'COS')
                % Vin(t) = A * cos(2*pi*f*t)
                A = numval(ee.A);
                f = numval(ee.f);
                vinst = A * cos(2*pi*f*tk);
                iVidx = idx.iV(ee.branch_col); %index of this voltage source branch current
                b_t(iVidx) = b_t(iVidx) + vinst;
            end
        end
    end

    % Backward Euler RHS: b_BE = b(t_k) + (C/h)*x_prev
    b_BE = b_t + (C / h) * x_prev;

    if ~nonlinear
        % Linear case no f(x)
        x_new = G_BE \ b_BE;
        info = struct('iters', 1, 'converged', true); %#ok<NASGU>
    else
        % Phi(x) = G_BE * x + f(x) - b_BE
        % J(x)   = G_BE + f'(x)
        x_init = x_prev;   % initial guess, initially all zeros

        %Xdim=length(x_prev);
        %x_init=zeros(Xdim, 1);

        %assign values


        [x_new, info] = nr_solve_BE(G_BE, b_BE, elems, node_map, x_init);
        if ~info.converged
            warning('Newton did not converge at t=%.4g (iters=%d).', tk, info.iters);
        end
    end

    Xhist(:,k+1) = x_new; %k+1 since first step is t=0
    x_prev = x_new;
end

fprintf('\nTransient simulation finished: t in [0, %.4g] s, h=%.4g s, steps=%d\n', ...
    t_end, h, nsteps);

%% ========= Example plotting for this deliverable =========
% Plot v(n1) and v(n2) if those nodes exist.
try
    n1 = node_lookup("N1", node_map);
    n2 = node_lookup("N2", node_map);

    if n1>0
        v1 = Xhist(n1,:).'; %v1 for all time steps
    else
        v1 = [];
    end
    if n2>0
        v2 = Xhist(n2,:).';
    else
        v2 = [];
    end

    figure;
    hold on;
    if ~isempty(v1), plot(t, v1, 'b', 'LineWidth', 2); end
    if ~isempty(v2), plot(t, v2, 'k', 'LineWidth', 2); end
    grid on;
    xlabel('Time (s)');
    ylabel('Voltage (V)');
    legend_entries = {};
    if ~isempty(v1), legend_entries{end+1} = 'v(n1)'; end
    if ~isempty(v2), legend_entries{end+1} = 'v(n2)'; end
    if ~isempty(legend_entries), legend(legend_entries, 'Location','best'); end
    title('Transient response');
    hold off;
catch
    % If something goes wrong with node names, just skip plotting
end

fprintf('===========================================\n\n');
end  % ===== end main function =====


%% ========================= Newton solver for BE =========================
function [x, did_it_converged] = nr_solve_BE(G_BE, b_BE, elems0, node_map0, x_init)
% Same Newton pattern as Deliverable 2:
%
%   1) f(x)
%   2) f'(x)
%   3) Phi(x) = G_BE * x + f(x) - b_BE
%   4) J(x)   = G_BE + f'(x)
%   5) Solve Δx = -J^(-1) * Phi
%   6) Update x := x + Δx
%   7) Converge if ||Δx||_inf <= 1e-6 (max 100 iters)

    maxit = 100;
    tol   = 1e-6;

    if nargin >= 5 && ~isempty(x_init) %nargin matlab function
        x = x_init;
    else
        x = zeros(size(b_BE));
    end

    for it = 1:maxit
        % (a) Compute f(x) and f'(x) from diodes
        [fvector, Jacobian] = diode_fx_Jacobian(x, elems0, node_map0, length(G_BE));

        % (b) Phi and J
        Phi = G_BE*x + fvector - b_BE;  % Phi(x)
        J   = G_BE + Jacobian;          % J(x)

        % (c) Newton step: J*dx = -Phi
        dx  = - (J \ Phi);

        % (d) Update
        x = x + dx;

        % (e) Convergence
        if norm(dx,inf) <= tol
            did_it_converged = struct('iters', it, 'converged', true);
            return;
        end
    end
    did_it_converged = struct('iters', maxit, 'converged', false);
end

%% ========================= f(x) and f'(x) for diodes =========================
function [fvector, Jacobian] = diode_fx_Jacobian(xv, elems0, node_map0, Xdim)
% same as DEV 2
% Diode: I = Is * (exp(Vd/Vt) - 1), Vd = v(n+) - v(n-)
% Defaults: Is = 1e-15 A, Vt = 25e-3 V

    Vt_default = 25e-3;
    fvector = zeros(Xdim,1);
    Jacobian = zeros(Xdim,Xdim);

    for kk = 1:numel(elems0)
        ee = elems0{kk}; if upper(ee.type) ~= 'D', continue; end
        [n1,n2] = nn(ee.np, ee.nn, node_map0);

        Is = ee.Is; if isempty(Is), Is = 1e-15; end
        Vt = ee.Vt; if isempty(Vt), Vt = Vt_default; end

        Vd = node_v(xv, n1) - node_v(xv, n2);
        arg = max(min(Vd / Vt, 40), -40);
        ev  = exp(arg);

        Id = Is*(ev - 1);      % current n1 -> n2
        didv = (Is / Vt)*ev;   % dI/dVd

        % f(x): diode current KCL contributions
        if n1>0, fvector(n1) = fvector(n1) + Id; end
        if n2>0, fvector(n2) = fvector(n2) - Id; end

        % f'(x): Jacobian contributions
        if n1>0
            Jacobian(n1,n1) = Jacobian(n1,n1) + didv;
            if n2>0, Jacobian(n1,n2) = Jacobian(n1,n2) - didv; end
        end
        if n2>0
            Jacobian(n2,n2) = Jacobian(n2,n2) + didv;
            if n1>0, Jacobian(n2,n1) = Jacobian(n2,n1) - didv; end
        end
    end
end

%% ========================= PASS COLLECT with .tran & COS =========================
%goal to scan list and get all necessary info for stamping
function [elems, node_map, counts, tran] = pass_collect_tran(lines)
    node_map = containers.Map('KeyType','char','ValueType','int32');
    elems = {};  nV=0; nL=0; nE=0; nH=0; nD=0; N=0;

    tran.enabled = false;
    tran.tstop   = 0;
    tran.dt      = 0;

    for i = 1:numel(lines)
        raw = char(lines(i));
        tok = tokenize(lines(i)); if isempty(tok), continue; end
        name = char(tok(1));

        % Handle '.tran' here and no stamps
        if name(1) == '.'
            cmd = upper(strrep(name,'.',''));
            switch cmd
                case 'TRAN'
                    if numel(tok) < 3
                        error('.tran line must be ".tran Tstop dt" on line %d: %s', i, raw);
                    end
                    tran.tstop = time_numval(tok(2));
                    tran.dt    = time_numval(tok(3));
                    tran.enabled = true;
                otherwise
                    % Ignore any other .xxx here (we already handled .op and .end)
            end
            continue;
        end

        type = upper(name(1));
        e = struct('raw',raw,'type',type,'name',name,'line_no',i);

        switch type
            case {'R','C'}
                [e.np, e.nn, e.value] = deal(tok(2), tok(3), tok(4));
                if type=='R'
                    vnum = safe_numval(e.value);
                    if ~isnan(vnum) && vnum==0
                        e.is_short = true; nV = nV + 1;
                    else
                        e.is_short = false;
                    end
                end

            case 'L'
                [e.np, e.nn, e.value] = deal(tok(2), tok(3), tok(4));
                % For transient only, we treat L as inductor
                e.dc_short = false;
                nL = nL + 1;

            case {'I','V'}
                e.np = tok(2); e.nn = tok(3);

                if upper(type) == 'V' && numel(tok)>=4 && strcmpi(tok(4),'COS')
                    % V name n+ n- COS A f
                    e.src_kind = 'COS';
                    if numel(tok) < 6
                        error('Voltage COS source must be "V name n+ n- COS A f" on line %d: %s', i, raw);
                    end
                    e.A = tok(5);
                    e.f = tok(6);
                    e.value = e.A;    % amplitude, used if needed
                    nV = nV + 1;
                else
                    % No cos case
                    if numel(tok) >= 5 && any(strcmpi(tok(4), {'DC','AC'}))
                        e.src_kind = upper(tok(4)); e.value = tok(5);
                    elseif numel(tok) >= 4 && ~any(strcmpi(tok(4), {'DC','AC'}))
                        e.src_kind = 'DC';         e.value = tok(4);
                    else
                        error('Missing numeric source value on line %d: %s', i, raw);
                    end
                    if type=='V', nV = nV + 1; end
                end

            case 'E'
                [e.np, e.nn, e.nc1, e.nc2, e.value] = deal(tok(2), tok(3), tok(4), tok(5), tok(6));
                nE = nE + 1;

            case 'G'
                [e.np, e.nn, e.nc1, e.nc2, e.value] = deal(tok(2), tok(3), tok(4), tok(5), tok(6));

            case 'F'
                [e.np, e.nn, e.ctrl, e.value] = deal(tok(2), tok(3), tok(4), tok(5));

            case 'H'
                [e.np, e.nn, e.ctrl, e.value] = deal(tok(2), tok(3), tok(4), tok(5));
                nH = nH + 1;

            case 'D'
                % Dname n+ n- [IS=...] [VT=...]
                e.np = tok(2); e.nn = tok(3);
                e.Is = 1e-15;
                e.Vt = 25e-3;
                for t = 4:numel(tok)
                    s = upper(strtrim(char(tok(t))));
                    if startsWith(s,'IS='), e.Is = numval(extractAfter(s,'IS='));
                    elseif startsWith(s,'VT='), e.Vt = numval(extractAfter(s,'VT='));
                    end
                end
                nD = nD + 1;

            otherwise
                error('Unsupported element "%s" on line %d: %s', type, i, raw);
        end

        elems{end+1} = e; %#ok<AGROW>

        % Node map for any referenced node (non-ground)
        for field = {'np','nn','nc1','nc2'}
            if isfield(e,field{1})
                nodename = upper(strtrim(string(e.(field{1})))); 
                if nodename=="" || nodename=="0" || nodename=="GND" || nodename=="GROUND", continue; end
                key = char(nodename);
                if ~isKey(node_map, key)
                    N = N + 1; node_map(key) = N;
                end
            end
        end
    end

    counts = struct('N',N,'nV',nV,'nL',nL,'nE',nE,'nH',nH,'nD',nD);
end


% Shared helpers ---------------------------------------------------
function lines = read_text_lines(fname)
    fid = fopen(fname,'r'); assert(fid>0, 'Cannot open %s', fname);
    cell_array = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    fclose(fid); lines = string(cell_array{1});
end

function [lines, has_op] = strip_comments_end_and_op(lines)
    result = strings(0); has_op = false;
    for i = 1:numel(lines)
        L = string(lines(i)); Lt = strtrim(L);
        if Lt == "" || startsWith(Lt,{'*',';'}), continue; end
        if startsWith(Lt, '.end', 'IgnoreCase', true), break; end
        if startsWith(Lt, '.op',  'IgnoreCase', true), has_op = true; continue; end
        cpos = min([char_pos(Lt,';'), char_pos(Lt,'*')]);
        if isfinite(cpos), Lt = strtrim(extractBefore(Lt, cpos)); end
        if Lt ~= "", result(end+1) = Lt; end %#ok<AGROW>
    end
    lines = result;
end

function p = char_pos(str, ch)
    rowvect = strfind(str, ch);
    if isempty(rowvect), p = inf; else, p = rowvect(1); end
end

function tok = tokenize(L)
    parts = regexp(char(L), '\s+', 'split');
    parts = parts(~cellfun(@isempty,parts));
    tok = string(parts);
end

function [n1, n2] = nn(ns1, ns2, node_map)
    n1 = node_lookup(ns1, node_map);
    n2 = node_lookup(ns2, node_map);
end

function idx = node_lookup(ns, node_map)
    ns = upper(strtrim(string(ns)));
    if ns=="0" || ns=="GND" || ns=="GROUND", idx = 0; return; end
    key = char(ns);
    if ~isKey(node_map, key), node_map(key) = node_map.Count + 1; end
    idx = node_map(key);
end

function M = stamp_passive(M, n1, n2, val)
    if n1>0, M(n1,n1) = M(n1,n1) + val; end
    if n2>0, M(n2,n2) = M(n2,n2) + val; end
    if n1>0 && n2>0
        M(n1,n2) = M(n1,n2) - val;
        M(n2,n1) = M(n2,n1) - val;
    end
end

function B = stamp_branch(B, n1, n2, k)
    if n1>0, B(n1,k) = B(n1,k) + 1; end
    if n2>0, B(n2,k) = B(n2,k) - 1; end
end

function v = numval(tok)
    s = strtrim(char(tok));
    [num, mult] = parse_eng_suffix(s);
    v = num * mult;
    if isnan(v), error('Bad numeric value: %s', s); end
end

function v = time_numval(tok)
    % Parse time values like 0.05s, 0.1ms, 1us using the same
    % engineering suffix rules as numval, but allowing a trailing 's'.
    s = strtrim(char(tok));

    % If it ends with 's' or 'S' (seconds), strip that unit.
    if ~isempty(s) && (s(end) == 's' || s(end) == 'S')
        s(end) = [];
    end

    % Now parse engineering suffix on the remaining part (k, m, u, n, p...)
    [num, mult] = parse_eng_suffix(s);
    v = num * mult;
    if isnan(v)
        error('Bad time numeric value: %s', s);
    end
end


function v = safe_numval(tok)
    try v = numval(tok); catch, v = NaN; end
end

function [num, mult] = parse_eng_suffix(s)
    mult = 1;
    if isempty(s), num = NaN; return; end
    s = upper(strtrim(s));
    if endsWith(s, "MEG"), mult = 1e6; s = extractBefore(s, strlength(s) - 2);
    elseif endsWith(s, "T"), mult = 1e12; s(end) = [];
    elseif endsWith(s, "G"), mult = 1e9;  s(end) = [];
    elseif endsWith(s, "K"), mult = 1e3;  s(end) = [];
    elseif endsWith(s, "M"), mult = 1e-3; s(end) = [];
    elseif endsWith(s, "U"), mult = 1e-6; s(end) = [];
    elseif endsWith(s, "N"), mult = 1e-9; s(end) = [];
    elseif endsWith(s, "P"), mult = 1e-12; s(end) = [];
    elseif endsWith(s, "F"), mult = 1e-15; s(end) = [];
    end
    num = str2double(s);
end

function val = source_value(e)
    val = numval(e.value);
end

function inv = invert_map(node_map, N)
    inv = strings(1,N);
    ks = keys(node_map); vs = values(node_map);
    for i=1:numel(vs), inv(vs{i}) = string(ks{i}); end
end

function v = node_v(x, n)
    if n>0, v = x(n); else, v = 0; end
end
```
