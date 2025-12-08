```matlab
function [G, C, b, meta] = Dev2_472_MH(spice_file)
format shortE
% DEV2_472_MH — A2 DC .op Newton
%
% Newton steps:
%   1) Compute f(x)
%   2) Compute f'(x)
%   3) Phi(x) = G*x + f(x) - b
%   4) J(x)   = G + f'(x)
%   5) Δx = -J^{-1} * Phi
%   6) x := x + Δx
%   7) Convergence if ||Δx||_inf <= 1e-6 (max 100 iters)
%

%% ========= PASS 1: read / strip / collect =========
lines = read_text_lines(spice_file);
[lines, has_op] = strip_comments_end_and_op(lines);
[elems, node_map, counts] = pass_collect(lines, has_op);   % parses D with Is, Vt only

% Dimensions 
N   = counts.N;
nV  = counts.nV;   % independent V + 0-ohm R + L-as-short in .op
nL  = counts.nL;   % time-domain inductors (not used in .op)
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
b = zeros(X,1);

% Branch incidence helpers
Bv = zeros(N,nV);  % independent V + 0-ohm R + L@.op
Bl = zeros(N,nL);  % inductors (time-domain)
Be = zeros(N,nE);  % VCVS
Bh = zeros(N,nH);  % CCVS

% Map "element name" -> absolute column index in x for branch-current unknowns
name2branch = containers.Map('KeyType','char','ValueType','int32');

%% Pre-assign branch columns
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
            if has_op && isfield(e,'dc_short') && e.dc_short
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

%% ========= PASS 2: stamp linear devices (diodes handled in Newton) =========
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
            if has_op && isfield(e,'dc_short') && e.dc_short
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
            if nplus>0,  b(nplus)  = b(nplus)  - Ival; end % leaves n+
            if nminus>0, b(nminus) = b(nminus) + Ival; end % enters n-

        case 'V'
            [n1,n2] = nn(e.np, e.nn, node_map);
            col = e.branch_col; Bv = stamp_branch(Bv, n1, n2, col);
            Vval = source_value(e);
            b(idx.iV(col)) = b(idx.iV(col)) + Vval;

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
            % handled in Newton

        otherwise
            error('Unsupported element type "%s" on line %d (%s).', T, e.line_no, e.raw);
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

% Labels & meta
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
meta.has_op    = has_op;

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
disp('b vector :'); disp(b);

%% ========= DC /.op SOLVE==============================
nonlinear = counts.nD > 0;
if has_op || nonlinear
    fprintf('\nNonlinear DC .op Newton\n');
    [x, did_it_converged] = nr_solve(G, b, elems, node_map);  % defined below
    if ~did_it_converged.converged
        warning('Newton did not converge in %d iterations; showing last iterate.', did_it_converged.iters);
    end
    disp(table(meta.var_names(:), x, 'VariableNames', {'Variable','Value'}));

    % report diode currents (optional)
    for k = 1:numel(elems)
        e = elems{k}; if upper(e.type)~='D', continue; end %skip anythig thats not diode
        [n1,n2] = nn(e.np, e.nn, node_map); %n1 positive node and n2 negative node
        Vd = node_v(x, n1) - node_v(x, n2); %voltage difference
        Is = e.Is; if isempty(Is), Is = 1e-15; end %VT and IS for diodes if not specialed
        Vt = e.Vt; if isempty(Vt), Vt = 25e-3; end
        arg = max(min(Vd/Vt, 40), -40); %ensure we have non infinite value arg is what is in the exponent
        Id = Is*(exp(arg) - 1);
        fprintf('Id(%s) = %.6g A  (Vd=%.4f V)\n', e.name, Id, Vd);
    end
else
   %just classic x=G/b
    x = G \ b;   % throws if singular, as desired
    disp(' '); disp('DC / operating-point solution (x = G\b):'); %\ is matlab built in system divider
    disp(table(meta.var_names(:), x, 'VariableNames', {'Variable','Value'}));
end

fprintf('===========================================\n\n');
end  % ===== end main function =====


%% ========================= Newton solver =========================
function [x, did_it_converged] = nr_solve(G0, b0, elems0, node_map0)
%   1) f(x)
%   2) f'(x)
%   3) Phi(x) = Gx + f(x) - b
%   4) J(x)   = G + f'(x)
%   5) Solve Δx = -J^(-1)* Phi
%   6) Update x := x + Δx
%   7) Converge if ||Δx||_inf <= 1e-6 (max 100 iters)

    maxit = 100;
    tol   = 1e-6;
    %x     = zeros(size(b0));  % Initial guess of zeros
    x=[50; 49; 0];

    for it = 1:maxit
        % (a) Compute f(x) and f'(x) 
        [fvector, Jacobian] = diode_fx_Jacobian(x, elems0, node_map0, length(G0)); %local function that computes f(x) and J (f'(x)) for diodes

        % (b) Phi and J
        Phi = G0*x + fvector - b0;    % Phi(x) = Gx + f(x) - b
        J   = G0 + Jacobian;           % J(x)   = G + f'(x)

        % (c) Newton step: J*dx = -Phi
        dx  = - (J \ Phi);         % \ matlab system solver
        %solves the Jdx=-Phi

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
% Diode: I = Is * (exp(Vd/Vt) - 1), Vd = v(n+) - v(n-)
% Defaults: Is = 1e-15 A, Vt = 25e-3 V

    Vt_default = 25e-3;
    fvector = zeros(Xdim,1);
    Jacobian = zeros(Xdim,Xdim);

    for kk = 1:numel(elems0)
        ee = elems0{kk}; if upper(ee.type) ~= 'D', continue; end %not a diode get out
        [n1,n2] = nn(ee.np, ee.nn, node_map0); %n1 pn and n2 nn (negative node)

        Is = ee.Is; if isempty(Is), Is = 1e-15; end %no specs to diode give default specs
        Vt = ee.Vt; if isempty(Vt), Vt = Vt_default; end

        Vd = node_v(xv, n1) - node_v(xv, n2); %node_v returns voltage value of node so that vd can be a number
        arg = max(min(Vd / Vt, 40), -40);  % prevent infinity e^40+
        ev  = exp(arg); %e^...

        Id = Is*(ev - 1);      % current n1 to n2
        didv = (Is / Vt)*ev;     % dI/dVd goes into jacobian

        % f(x)
        %DIODE Stamp======================================================
        if n1>0, fvector(n1) = fvector(n1) + Id; end
        if n2>0, fvector(n2) = fvector(n2) - Id; end

        % f'(x) whihc is the acobian
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

%% ========================= Shared helpers =========================
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

function [elems, node_map, counts] = pass_collect(lines, has_op)
    node_map = containers.Map('KeyType','char','ValueType','int32');
    elems = {};  nV=0; nL=0; nE=0; nH=0; nD=0; N=0;

    for i = 1:numel(lines)
        raw = char(lines(i));
        tok = tokenize(lines(i)); if isempty(tok), continue; end
        name = char(tok(1)); type = upper(name(1));
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
                if has_op
                    e.dc_short = true; nV = nV + 1;  % L shorted in .op
                else
                    e.dc_short = false; nL = nL + 1;
                end

            case {'I','V'}
                e.np = tok(2); e.nn = tok(3);
                if numel(tok) >= 5 && any(strcmpi(tok(4), {'DC','AC'}))
                    e.src_kind = upper(tok(4)); e.value = tok(5);
                elseif numel(tok) >= 4 && ~any(strcmpi(tok(4), {'DC','AC'}))
                    e.src_kind = 'DC';         e.value = tok(4);
                else
                    error('Missing numeric source value on line %d: %s', i, raw);
                end
                if type=='V', nV = nV + 1; end

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
                % Dname n+ n- [IS=...] [VT=...] (defaults below)
                e.np = tok(2); e.nn = tok(3);
                e.Is = 1e-15;   % default Is (10e-16)
                e.Vt = 25e-3;   % default Vt (25 mV)
                for t = 4:numel(tok)
                    s = upper(strtrim(char(tok(t))));
                    if startsWith(s,'IS='), e.Is = numval(extractAfter(s,'IS='));
                    elseif startsWith(s,'VT='), e.Vt = numval(extractAfter(s,'VT='));
                    % no N=...
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
