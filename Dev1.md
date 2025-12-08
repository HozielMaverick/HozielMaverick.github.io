```matlab
function [G, C, b, meta] = Dev1_472_MH(spice_file)
format shortE 
% DEV1_472_MH  Parse a SPICE-like netlist and assemble MNA matrices.
%   [G, C, b, meta] = Dev1_472_MH('circuit.sp' or 'circuit.txt')
%
% MNA form:  G*x + C*xdot = b
% Supported: R, C, L, I, V, E(VCVS), G(VCCS), F(CCCS), H(CCVS), 0-ohm R→short
% Ground aliases: 0, GND, GROUND
% Values: SPICE suffixes (case-insensitive): T,G,MEG,K,M,U,N,P,F


% Two-pass flow:
%   PASS 1 collect
%   PASS 2 stamp


% Pass 1:
% Read text, remove comments/.end, and collect elements + bookkeeping.
lines = read_text_lines(spice_file); %load file content into string array 1 per element
lines = strip_comments_and_end(lines);
[elems, node_map, counts] = pass_collect(lines); %tok each lines non ground nodes

% Dimensions
N  = counts.N;    % non-ground nodes first N unknowns are node voltages
nV = counts.nV;   % independent Voltage sources 
nL = counts.nL;   % inductors 
nE = counts.nE;   % VCVS E
nH = counts.nH;   % CCVS H
X  = N + nV + nL + nE + nH; % length of unknown vector x

% Indexing
% Ordering (fixed): [ v(1..N), iV(1..nV), iL(1..nL), iE(1..nE), iH(1..nH) ]
%v voltage source iv branche currents iL inductor branche current Ie VCVS
%branche currents iH CCVS branche currents
idx.v  = 1:N;
idx.iV = N + (1:nV);
idx.iL = N + nV + (1:nL);
idx.iE = N + nV + nL + (1:nE);
idx.iH = N + nV + nL + nE + (1:nH);

% Allocate
% Final matrix sizes are now known; initialize dense matrices (could be sparse).
G = zeros(X,X);   % Conductance / topology 
C = zeros(X,X);   % Dynamic terms
b = zeros(X,1);   % RHS vector

% Helper accumulators (branch incidence)
Bv = zeros(N,nV);  % for V and 0-ohm R shorts
Bl = zeros(N,nL);  % for L
Be = zeros(N,nE);  % for VCVS
Bh = zeros(N,nH);  % for CCVS
 
% map element NAMEs to absolute column indices in x for any element
% creates a branch-current unknown
% name2branch: v1 to idx in x for i(V1)
name2branch = containers.Map('KeyType','char','ValueType','int32');

vcolcount = 0; ecolcount = 0; hcolcount = 0;
for k = 1:numel(elems)
    elements_of_struct = elems{k}; T = upper(elements_of_struct.type);
    switch T
        case 'V'
            % Independent V-source: allocate next iV slot and remember it by name
            vcolcount = vcolcount + 1;
            elements_of_struct.branch_col = vcolcount;
            name2branch(upper(elements_of_struct.name)) = idx.iV(vcolcount);

        case 'R'
            % 0-ohm resistor is treated as an ideal short acts like a V-source branch
            if isfield(elements_of_struct,'is_short') && elements_of_struct.is_short
                vcolcount = vcolcount + 1;
                elements_of_struct.branch_col = vcolcount;
                name2branch(upper(elements_of_struct.name)) = idx.iV(vcolcount);
            else
                elements_of_struct.branch_col = 0; % ordinary resistor has no branch unknown
            end

        case 'E'
            % VCVS introduces a branch current (iE); allocate and map by name
            ecolcount = ecolcount + 1;
            elements_of_struct.branch_col = ecolcount;
            name2branch(upper(elements_of_struct.name)) = idx.iE(ecolcount);

        case 'H'
            % CCVS introduces a branch current (iH); allocate and map by name
            hcolcount = hcolcount + 1;
            elements_of_struct.branch_col = hcolcount;
            name2branch(upper(elements_of_struct.name)) = idx.iH(hcolcount);

        otherwise
            % Other devices do not introduce a branch-current unknown
            elements_of_struct.branch_col = 0;
    end
    elems{k} = elements_of_struct; % write back the annotated element no structural changes
end

%Pass 2:
% Now that indices are fixed and branch columns are known, time to stamp G
% C b
lcol_seen = 0; % running index for inductors

for k = 1:numel(elems)
    elements_of_struct = elems{k}; T = upper(elements_of_struct.type);

    switch T
        case 'R' % Resistor (0-ohm treated as ideal short branch like V)
            % Ordinary R, G nodal stamp; 0-ohm, treat as a V-branch with KVL=0
            [n1,n2] = nn(elements_of_struct.np, elements_of_struct.nn, node_map);
            if isfield(elements_of_struct,'is_short') && elements_of_struct.is_short
                col = elements_of_struct.branch_col;                    
                Bv = stamp_branch(Bv, n1, n2, col);    
                
            else
                g = 1/numval(elements_of_struct.value);                 % conductance
                G = stamp_passive(G, n1, n2, g);       % standard R stamp
            end

        case 'C' % Capacitor like g in G but in C
            [n1,n2] = nn(elements_of_struct.np, elements_of_struct.nn, node_map);
            Cval = numval(elements_of_struct.value);
            C(idx.v, idx.v) = stamp_passive(C(idx.v,idx.v), n1, n2, Cval);

        case 'L' % Inductor add iL
            [n1,n2] = nn(elements_of_struct.np, elements_of_struct.nn, node_map);
            lcol_seen = lcol_seen + 1;                 % which inductor is this?
            Bl = stamp_branch(Bl, n1, n2, lcol_seen);  % KCL incidence for iL
            Lval = numval(elements_of_struct.value);
            il = idx.iL(lcol_seen);                    % absolute index of iL in x
            C(il, il) = C(il, il) + Lval;              % L * diL/dt (state equation)

        case 'I' % Independent current source (DC value)
            
            [n1,n2] = nn(elements_of_struct.np, elements_of_struct.nn, node_map);
            Ival = source_value(elements_of_struct);
            if n1>0, b(n1) = b(n1) + Ival; end
            if n2>0, b(n2) = b(n2) - Ival; end
            if isfield(elements_of_struct,'src_kind') && strcmpi(elements_of_struct.src_kind,'AC')
                % For this deliverable, AC is treated as DC magnitude; warn once per source.
                warning('AC source found; ignored (treated as DC magnitude).');
            end

        case 'V' % Independent voltage source: add iV and KVL RHS
            % Creates branch current iV and a KVL row with b = V.
            [n1,n2] = nn(elements_of_struct.np, elements_of_struct.nn, node_map);
            col = elements_of_struct.branch_col;                        % known from pre-assign
            Bv = stamp_branch(Bv, n1, n2, col);        % KCL incidence for iV
            Vval = source_value(elements_of_struct);
            b(idx.iV(col)) = b(idx.iV(col)) + Vval;    % KVL RHS
            if isfield(elements_of_struct,'src_kind') && strcmpi(elements_of_struct.src_kind,'AC')
                warning('AC source found; ignored (treated as DC magnitude).');
            end

        case 'E' % VCVS: v(n1)-v(n2) = A*(v(nc1)-v(nc2))
            % Introduces branch current iE  + controlled KVL row.
            [n1,n2]   = nn(elements_of_struct.np,  elements_of_struct.nn,  node_map);
            [nc1,nc2] = nn(elements_of_struct.nc1, elements_of_struct.nc2, node_map);
            col = elements_of_struct.branch_col;                        % known from pre-assign
            Be = stamp_branch(Be, n1, n2, col);        % KCL incidence for iE
            A = numval(elements_of_struct.value);
            ie = idx.iE(col);                          % KVL row index
            % KVL row: v(n1)-v(n2) - A*(v(nc1)-v(nc2)) = 0
            if nc1>0, G(ie, nc1) = G(ie, nc1) - A; end
            if nc2>0, G(ie, nc2) = G(ie, nc2) + A; end

        case 'G' % VCCS: i = g*(v(nc1)-v(nc2)) injected from n+ to n-
            
            [n1,n2]   = nn(elements_of_struct.np,  elements_of_struct.nn,  node_map);
            [nc1,nc2] = nn(elements_of_struct.nc1, elements_of_struct.nc2, node_map);
            gm = numval(elements_of_struct.value);
            % Four entries: +g at (n+,nc1), -g at (n+,nc2), -g at (n-,nc1), +g at (n-,nc2)
            if n1>0 && nc1>0, G(n1,nc1) = G(n1,nc1) + gm; end
            if n1>0 && nc2>0, G(n1,nc2) = G(n1,nc2) - gm; end
            if n2>0 && nc1>0, G(n2,nc1) = G(n2,nc1) - gm; end
            if n2>0 && nc2>0, G(n2,nc2) = G(n2,nc2) + gm; end

        case 'F' % CCCS: Fname n+ n- Vctrl beta  (A/A)
            
            [n1,n2] = nn(elements_of_struct.np, elements_of_struct.nn, node_map);
            beta = numval(elements_of_struct.value);
            ctrl = upper(char(elements_of_struct.ctrl));                % control element name
            assert(isKey(name2branch, ctrl), ...
                'CCCS %s controls unknown branch "%s". Control must be a voltage-source-like element.', elements_of_struct.name, ctrl);
            ictrl = name2branch(ctrl);                 % absolute index of i_ctrl in x
            if n1>0, G(n1, ictrl) = G(n1, ictrl) + beta; end
            if n2>0, G(n2, ictrl) = G(n2, ictrl) - beta; end

        case 'H' % CCVS: Hname n+ n- Vctrl Rm  (V/A)
            
            [n1,n2] = nn(elements_of_struct.np, elements_of_struct.nn, node_map);
            col = elements_of_struct.branch_col;                        % known from pre-assign
            Bh = stamp_branch(Bh, n1, n2, col);        % KCL incidence for iH
            Rm = numval(elements_of_struct.value);
            ih = idx.iH(col);                          % KVL row index for iH
            ctrl = upper(char(elements_of_struct.ctrl));
            assert(isKey(name2branch, ctrl), ...
                'CCVS %s controls unknown branch "%s". Control must be a voltage-source-like element.', elements_of_struct.name, ctrl);
            ictrl = name2branch(ctrl);
            % place -Rm in (ih, ictrl)
            G(ih, ictrl) = G(ih, ictrl) - Rm;

        otherwise
            error('Unsupported element type "%s" on line %d (%s).', T, elements_of_struct.line_no, elements_of_struct.raw);
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


% Build human-readable labels for x 
var_names = strings(1,X);
inv_node = invert_map(node_map, N); %inverse map means node to name as opposed to name to node
for i = 1:N, var_names(i) = "v(" + inv_node(i) + ")"; end
for i = 1:nV, var_names(N+i) = "i(V" + i + ")"; end %current
for i = 1:nL, var_names(N+nV+i) = "i(L" + i + ")"; end %inductor current
for i = 1:nE, var_names(N+nV+nL+i) = "i(E" + i + ")"; end %VCVS Branch I
for i = 1:nH, var_names(N+nV+nL+nE+i) = "i(H" + i + ")"; end %CCVS branch I

%Struct meta
meta.node_map  = node_map; %name to index
meta.var_names = var_names; %lbel for x
meta.idx       = idx; %index range
meta.counts    = counts; %N, nV,...
meta.elems     = elems; %parsed elements 

% Summary who prints everytime
fprintf('\n=== MNA Summary for "%s" ===\n', spice_file);
fprintf('Nodes (non-ground): %d\n', meta.counts.N);
fprintf('Voltage sources:   %d\n', meta.counts.nV);
fprintf('Inductors:         %d\n', meta.counts.nL);
fprintf('VCVS:              %d\n', meta.counts.nE);
fprintf('CCVS (H):          %d\n', meta.counts.nH);

disp(' ');
disp('Variable ordering (x):');
disp(meta.var_names.');

disp(' '); %blank line 
disp('G matrix:'); disp(G);
disp('C matrix:'); disp(C);
disp('b vector:'); disp(b);

% Always attempt a DC/operating-point solve (steady state: xdot = 0)
% If G is singular, warn (common with floating nodes / insufficient references).
try
    x = G \ b;
    disp(' ');
    disp('DC / operating-point solution (x = G\\b):');
    disp(table(meta.var_names(:), x, 'VariableNames', {'Variable','Value'}));
catch
    disp(' '); %blank line
    disp('Matrix G is singular cannot compute');
end

fprintf('===========================================\n\n');

end

% Helpers

function lines = read_text_lines(fname)
    % Read entire text file into a string array, one line per element
    fid = fopen(fname,'r'); %file ID
    assert(fid>0, 'Cannot open %s', fname); %makes sure it opens
    cell_array = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', ''); %built in matlab function
    fclose(fid);
    lines = string(cell_array{1});
end

function lines = strip_comments_and_end(lines)
    %cleans raw line 
    result = strings(0); %empty string array
    for i = 1:numel(lines)
        L = string(lines(i));
        Lt = strtrim(L); %trimmed line
        if Lt == "" || startsWith(Lt,{'*',';'}), continue; end %skip if Lt is empty or starts with...
        if startsWith(Lt, '.end', 'IgnoreCase', true), break; end %if line starts with .end get out
        % remove inline comments
        cpos = min([char_pos(Lt,';'), char_pos(Lt,'*')]); %comment position
        %remove first one

        if isfinite(cpos), Lt = strtrim(extractBefore(Lt, cpos)); end %matlab functions
        %if lT is not empty append lt to end of result and telling matlab i
        %know my array can grow inside a loop
        if Lt ~= "", result(end+1) = Lt; end %#ok<AGROW>

    end
    lines = result;
end

function p = char_pos(str, ch) %char position
    % Return first position of character ch in str inf if not present
    rowvect = strfind(str, ch); %matlab function
    if isempty(rowvect), p = inf; else, p = rowvect(1); end
end

function [elems, node_map, counts] = pass_collect(lines)
    % used for pass 1
    % Tokenize non-comment line, classify element type and store it
    % fields (nodes, params, control), and:
    %   • count branch-creating devices (nV, nL, nE, nH)
    %   • build node_map (name → 1..N) for all non-ground nodes
    % No matrices are stamped here—only metadata is collected.
    node_map = containers.Map('KeyType','char','ValueType','int32');
    elems = {};  nV=0; nL=0; nE=0; nH=0; N=0;

    for i = 1:numel(lines)
        raw = char(lines(i));
        tok = tokenize(lines(i)); if isempty(tok), continue; end
        name = char(tok(1)); type = upper(name(1)); %name R1 type R
        e = struct('raw',raw,'type',type,'name',name,'line_no',i); %creating struct e

        switch type
            case {'R','C','L'}
                % Passive 2 terminal name n+ n- value
                [e.np, e.nn, e.value] = deal(tok(2), tok(3), tok(4));
                if type=='R'
                    % Treat 0-ohm R as ideal short (reserve iV)
                    vnum = safe_numval(e.value);
                    if ~isnan(vnum) && vnum==0 %parsing good and value is 0
                        e.is_short = true; nV = nV + 1;
                    else
                        e.is_short = false;
                    end
                else
                    e.is_short = false;
                end
                if type=='L', nL = nL + 1; end

            case {'I','V'}
                % Independent sources name n+ n- [DC/AC] value
                e.np = tok(2); e.nn = tok(3); %numel counts nb of elem matlab F
                if numel(tok)>=5 && any(strcmpi(tok(4),{'DC','AC'})) %strcmpi matlab function true if vlues = not case sensitive
                    e.src_kind = upper(tok(4)); e.value = tok(5);
                else
                    e.src_kind = 'DC'; e.value = tok(4); %if it doesnt say assume DC
                end
                if type=='V', nV = nV + 1; end %voltage source current
                e.is_short = false;

            case 'E' % VCVS: name n+ n- nc+ nc- gain
                [e.np, e.nn, e.nc1, e.nc2, e.value] = deal(tok(2), tok(3), tok(4), tok(5), tok(6));
                nE = nE + 1; e.is_short = false;

            case 'G' % VCCS: name n+ n- nc+ nc- gm
                [e.np, e.nn, e.nc1, e.nc2, e.value] = deal(tok(2), tok(3), tok(4), tok(5), tok(6));
                e.is_short = false;

            case 'F' % CCCS: name n+ n- Vctrl beta
                [e.np, e.nn, e.ctrl, e.value] = deal(tok(2), tok(3), tok(4), tok(5));
                e.is_short = false;

            case 'H' % CCVS: name n+ n- Vctrl Rm
                [e.np, e.nn, e.ctrl, e.value] = deal(tok(2), tok(3), tok(4), tok(5));
                nH = nH + 1; e.is_short = false;

            otherwise
                error('Unsupported element "%s" on line %d: %s', type, i, raw);
        end

        elems{end+1} = e; %#ok<AGROW>
        %tells matlab to not give me warningg
        %appends struct e to end of elems

        % Ground nodes handling
        for field = {'np','nn','nc1','nc2'} %+ - nodes % control nodes
            if isfield(e,field{1}) %check if it exist matlab function
                %strtrim matlab F removes white space
                nodename = upper(strtrim(string(e.(field{1})))); % node name
                %disregard ground
                if nodename=="" || nodename=="0" || nodename=="GND" || nodename=="GROUND", continue; end
                key = char(nodename);
                if ~isKey(node_map, key)
                    N = N + 1; node_map(key) = N; %add node if it wasnt already there
                end
            end
        end
    end

    counts = struct('N',N,'nV',nV,'nL',nL,'nE',nE,'nH',nH); %struct counts that counts
end

function tok = tokenize(L)
    % Split a line on whitespace into tokens, keeping order.
    parts = regexp(char(L), '\s+', 'split'); %reg expression tokenize all even white space
    parts = parts(~cellfun(@isempty,parts)); %only keep not empty tokens cellfun matlab function
    tok = string(parts);
end

function [n1, n2] = nn(ns1, ns2, node_map) %node number helper
    % Convenience wrapper look up two nodes (returns 0 for ground).
    n1 = node_lookup(ns1, node_map);
    n2 = node_lookup(ns2, node_map);
end

function idx = node_lookup(ns, node_map) %lookup helper

    ns = upper(strtrim(string(ns))); %node string
    if ns=="0" || ns=="GND" || ns=="GROUND", idx = 0; return; end
    key = char(ns);
    if ~isKey(node_map, key)
        node_map(key) = node_map.Count + 1;
    end
    idx = node_map(key); %read index
end

function M = stamp_passive(M, n1, n2, val) %passive elements 
    % stamps for g and C since same way to stamp in G and C matrix
    if n1>0, M(n1,n1) = M(n1,n1) + val; end %terminal not gnd? yes= no diaogonal stamping
    if n2>0, M(n2,n2) = M(n2,n2) + val; end
    if n1>0 && n2>0 %diagonal stamping now anly take care of -g/c
        M(n1,n2) = M(n1,n2) - val;
        M(n2,n1) = M(n2,n1) - val;
    end
end

function B = stamp_branch(B, n1, n2, k)
    % Fill row n at column k
    % +1 at node n1, -1 at node n2 no ground 
    if n1>0, B(n1,k) = B(n1,k) + 1; end
    if n2>0, B(n2,k) = B(n2,k) - 1; end
end

function v = numval(tok)
    % Convert a token with magnitude unit
    s = strtrim(char(tok));
    [num, mult] = parse_eng_suffix(s);
    v = num * mult;
    if isnan(v), error('Bad numeric value: %s', s); end
end

function v = safe_numval(tok)
    % Like numval, but returns NaN instead of throwing (used for is_short check).
    try v = numval(tok); catch, v = NaN; end
end

function [num, mult] = parse_eng_suffix(s)
    % SPICE magnitudes
    mult = 1;
    if isempty(s), num = NaN; return; end
    s = upper(strtrim(s));
    if endsWith(s, "MEG") %Since SpICE is not case sensitive
        mult = 1e6; s = extractBefore(s, strlength(s) - 2);
    elseif endsWith(s, "T")
        mult = 1e12; s(end) = [];
    elseif endsWith(s, "G")
        mult = 1e9; s(end) = [];
    elseif endsWith(s, "K")
        mult = 1e3; s(end) = [];
    elseif endsWith(s, "M")
        mult = 1e-3; s(end) = [];
    elseif endsWith(s, "U")
        mult = 1e-6; s(end) = [];
    elseif endsWith(s, "N")
        mult = 1e-9; s(end) = [];
    elseif endsWith(s, "P")
        mult = 1e-12; s(end) = [];
    elseif endsWith(s, "F")
        mult = 1e-15; s(end) = [];
    end
    num = str2double(s); %double type
end

function val = source_value(e) %makes numeric value
    val = numval(e.value); % For Deliverable 1 AC is DC magntiude 1
end

function inv = invert_map(node_map, N)
    % Build inverse mapping index to node name, no change order just builds
    % link. Used to name elements once circuit solved
    inv = strings(1,N);
    ks = keys(node_map); vs = values(node_map);
    for i=1:numel(vs), inv(vs{i}) = string(ks{i}); end
end
```
