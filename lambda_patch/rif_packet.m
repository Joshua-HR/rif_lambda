function P = rif_packet(R, varargin)
%RIF_PACKET  Packet-level decision of the RIF validity detector: combine the fragment statistics of the
%   RIF fragments that one receiver uses for one ranging measurement (NumRIF = 4 -> 8 RIF fragments).
%
%   The requirement is set per packet (FiRa "medium": false acceptance 1e-6 per ranging), and the
%   fragment statistic Lam_f of rif_cir_analyze is chi-square with 2 DOF under H0 for both chips:
%       new chip   ('Chip','new',   5-column dumps)  Lambda from the raw-sample moments Q, Pi
%       fixed chip ('Chip','fixed', any dump)        chi2(2) equivalent of Lambda-hat (Q, Pi from the CIR)
%   The STS changes every symbol, so the fragments of a packet are independent under H0.
%
%   P = rif_packet(R)                   R: rif_cir_analyze output of ONE job folder (needs Lam_W)
%   P = rif_packet(R, Name, Value, ...)
%       'Combine'       'soft'    'soft' | 'and' | 'kofn' | 'strict' | 'single'   (see below)
%       'FragPerPacket' 8         RIF fragments per packet. The fragments of one frame (sorted by sample) are
%                                 split into consecutive groups of this size; incomplete groups are dropped.
%                                 1 = one-fragment configuration (every rule reduces to 'single')
%       'PfaTarget'     1e-6      packet false-acceptance target (not used by 'strict')
%       'KofN'          []        'kofn': fragments that must pass ([] = FragPerPacket - 2, at least 1)
%       'PfaFragment'   1e-3      'strict': design Pfa of the fragment threshold
%       'Floor'         1         'soft': every fragment needs Lam >= Floor on the selected path, so a packet
%                                 with an empty (jammed / replaced) fragment is discarded; 0 = off
%       'Drift'         0         'soft': peak-drift hypotheses [taps per fragment], e.g. -0.5:1/7:0.5.
%                                 Fragment f of a path starting at tap j uses tap j + round(d (f - 1)).
%       'Margin'        1.5       threshold margin (design for target / Margin)
%
%   Rules (W = |W_s| taps, K = FragPerPacket, Lmax_f = max over W_s of Lam_f)
%       single : first fragment only, Lmax_1 >= T,     W exp(-T/2) = PfaTarget / Margin
%       and    : every fragment Lmax_f >= T_f,          per-fragment budget p = PfaTarget^(1/K)
%       kofn   : at least KofN fragments Lmax_f >= T_f, p with P(Binomial(K, p) >= KofN) = PfaTarget
%       strict : every fragment at the fragment threshold of PfaFragment (packet Pfa ~ PfaFragment^K)
%       soft   : S = max over start taps j and drifts d of sum_f Lam_f[j + round(d (f-1))] >= T with
%                W N_d P(chi2_2K > T) = PfaTarget / Margin, then the Floor check on the selected path
%       T_f = 2 ln(W / p) + 2 ln(Margin) for the fragment rules.
%
%   Returned struct
%       rule, K, KofN, T (packet or fragment threshold), pFrag (per-fragment budget; NaN for soft),
%       Tchk (fragment threshold used for the independence check: T_f, or per-fragment 0.2 for soft / single),
%       drift, nPkt, nDrop (fragments in incomplete groups), valid, stat, kSel (LLS tap of the selected
%       path in the first fragment, or of the strongest fragment), dSel, floorFail, frameOf,
%       fragHit / nFrag (fragments with Lmax >= Tchk), pairHit / nPair (disjoint pairs of fragments of one
%       packet both >= Tchk: about (fragHit/nFrag)^2 if independent), pairSum = [n sx sy sxx syy sxy] of the
%       Lmax pairs and rhoPair (their correlation, ~0 +- 1/sqrt(n) if independent), tailT / tailCnt / nPath
%       (soft: zero-drift paths with S >= t, to compare with P(chi2_2K >= t)), opt

p = inputParser;
p.addRequired('R', @isstruct);
p.addParameter('Combine', 'soft', @(x) ischar(x) && any(strcmpi(x, {'soft', 'and', 'kofn', 'strict', 'single'})));
p.addParameter('FragPerPacket', 8, @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == round(x));
p.addParameter('PfaTarget', 1e-6, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('KofN', [], @isnumeric);
p.addParameter('PfaFragment', 1e-3, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('Floor', 1, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('Drift', 0, @isnumeric);
p.addParameter('Margin', 1.5, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.parse(R, varargin{:});
o = rmfield(p.Results, 'R');                       % options only (R holds the CIR arrays)
if ~isfield(R, 'Lam_W') || isempty(R.Lam_W)
    error('rif_packet:nolam', ['No Lambda in the input. Run rif_cir_analyze with ''Chip'', ''fixed'' ' ...
          '(Lambda-hat from the CIR, any dump) or ''Chip'', ''new'' (5-column dumps).']);
end
K = o.FragPerPacket;
W = size(R.Lam_W, 1);
rule = lower(o.Combine);
if K == 1, rule = 'single'; end
kN = o.KofN;
if isempty(kN), kN = max(1, K - 2); end
if strcmp(rule, 'kofn') && (kN < 1 || kN > K || kN ~= round(kN))
    error('rif_packet:kofn', '''KofN'' must be an integer in 1..%d.', K);
end
d = o.Drift(:).';
if isempty(d), d = 0; end
nd = numel(d);
sig1 = 119;
if isfield(R, 'opt') && isfield(R.opt, 'Sig'), sig1 = R.opt.Sig(1); end

% ---------------------------------------------------------------- packets
frame = R.frame(:).';
frag = R.frag(:).';
[pid, nDrop] = packet_ids(frame, frag, K);
nPkt = max([pid 0]);
if nPkt == 0 && ~isempty(frame)
    warning('rif_packet:nopacket', ['No frame holds %d fragments: no packet formed (check ''FragPerPacket'' ' ...
            'against the fragments per frame, e.g. with check_packet_assumptions).'], K);
end
pos = mod(frag - 1, K) + 1;                        % position of the fragment in its packet
L3 = zeros(W, K, nPkt);
frameOf = zeros(1, nPkt);
for i = find(pid > 0)
    L3(:, pos(i), pid(i)) = R.Lam_W(:, i);
    frameOf(pid(i)) = frame(i);
end
Lmax3 = reshape(max(L3, [], 1), K, nPkt);          % K x nPkt
[~, im] = max(L3, [], 1);
jMax3 = reshape(im, K, nPkt);                      % window index of each fragment's maximum
Tf = @(pf) 2 * log(W / pf) + 2 * log(o.Margin);

valid = false(1, nPkt); stat = zeros(1, nPkt); jSel = ones(1, nPkt); dSel = zeros(1, nPkt);
floorFail = 0; tailT = []; tailCnt = []; nPath = 0;
switch rule
    case 'single'
        pFrag = o.PfaTarget; T = Tf(pFrag); Tchk = Tf(0.2);
        stat = Lmax3(1, :); jSel = jMax3(1, :);
    case 'and'
        pFrag = o.PfaTarget ^ (1 / K); T = Tf(pFrag); Tchk = T;
        stat = min(Lmax3, [], 1);
    case 'kofn'
        pFrag = kofn_budget(kN, K, o.PfaTarget); T = Tf(pFrag); Tchk = T;
        s = sort(Lmax3, 1, 'descend');
        stat = s(kN, :);
    case 'strict'
        pFrag = o.PfaFragment; T = Tf(pFrag); Tchk = T;
        stat = min(Lmax3, [], 1);
    case 'soft'
        pFrag = NaN; Tchk = Tf(0.2);
        T = soft_threshold(K, W * nd, o.PfaTarget / o.Margin);
        stat = -Inf(1, nPkt);
        for a = 1:nd
            sh = round(d(a) * (0:K - 1));          % tap shift of fragment f relative to fragment 1
            j0 = max(1, 1 - min(sh)); j1 = min(W, W - max(sh));
            if j0 > j1, continue; end              % no path of this drift stays inside W_s
            js = (j0:j1).';
            S = zeros(numel(js), nPkt);
            for f = 1:K
                S = S + reshape(L3(js + sh(f), f, :), numel(js), nPkt);
            end
            [m, ix] = max(S, [], 1);
            up = m > stat;
            stat(up) = m(up); jSel(up) = js(ix(up)).'; dSel(up) = d(a);
        end
        S0 = reshape(sum(L3, 2), W, nPkt);         % zero-drift paths, for the chi2_2K tail check
        tailT = 0:ceil(T);
        if nPkt > 0
            cnt = histcounts(S0(:), [tailT Inf]);
            tailCnt = fliplr(cumsum(fliplr(cnt)));  % paths with S >= tailT(i)
        else
            tailCnt = zeros(size(tailT));
        end
        nPath = W * nPkt;
end
valid = stat >= T;
if ~strcmp(rule, 'soft') && ~strcmp(rule, 'single') && nPkt > 0
    [~, fBest] = max(Lmax3, [], 1);                % position reported from the strongest fragment
    jSel = jMax3(sub2ind([K nPkt], fBest, 1:nPkt));
end
if strcmp(rule, 'soft') && o.Floor > 0
    ok = true(1, nPkt);
    for q = find(valid)
        jj = jSel(q) + round(dSel(q) * (0:K - 1));
        v = L3(sub2ind([W K nPkt], jj, 1:K, q * ones(1, K)));
        ok(q) = all(v >= o.Floor);
    end
    floorFail = sum(valid & ~ok);
    valid = valid & ok;
end

% ---------------------------------------------------------------- independence check (disjoint pairs)
hits = Lmax3 >= Tchk;
np = floor(K / 2);
pairHit = 0; x = zeros(1, 0); y = zeros(1, 0);
for a = 1:np
    pairHit = pairHit + sum(hits(2 * a - 1, :) & hits(2 * a, :));
    x = [x, Lmax3(2 * a - 1, :)]; y = [y, Lmax3(2 * a, :)]; %#ok<AGROW>
end
pairSum = [numel(x), sum(x), sum(y), sum(x .^ 2), sum(y .^ 2), sum(x .* y)];

P = struct('rule', rule, 'K', K, 'KofN', kN, 'T', T, 'pFrag', pFrag, 'Tchk', Tchk, 'drift', d, ...
           'nPkt', nPkt, 'nDrop', nDrop, 'valid', valid, 'stat', stat, 'kSel', sig1 + jSel - 1, 'dSel', dSel, ...
           'floorFail', floorFail, 'frameOf', frameOf, 'fragHit', sum(hits(:)), 'nFrag', K * nPkt, ...
           'pairHit', pairHit, 'nPair', np * nPkt, 'pairSum', pairSum, 'rhoPair', pair_rho(pairSum), ...
           'tailT', tailT, 'tailCnt', tailCnt, 'nPath', nPath, 'opt', o);
end


% ======================================================================
function r = pair_rho(s)
% correlation coefficient from the sums [n sx sy sxx syy sxy] (NaN for fewer than 3 pairs)
r = NaN;
if s(1) < 3, return; end
n = s(1);
cxy = s(6) - s(2) * s(3) / n;
vx = s(4) - s(2) ^ 2 / n;
vy = s(5) - s(3) ^ 2 / n;
if vx > 0 && vy > 0, r = cxy / sqrt(vx * vy); end
end


% ======================================================================
function [pid, nDrop] = packet_ids(frame, frag, K)
% packet number of every fragment (0 = in an incomplete group): consecutive groups of K fragments per frame
n = numel(frame);
pid = zeros(1, n); nDrop = 0;
if n == 0, return; end
key = [frame(:), floor((frag(:) - 1) / K)];
[~, ~, ic] = unique(key, 'rows', 'stable');
cnt = accumarray(ic, 1);
full = cnt == K;
newId = zeros(size(cnt));
newId(full) = 1:sum(full);
pid = newId(ic).';
nDrop = sum(cnt(~full));
end


% ======================================================================
function q = chi2tail(T, K)
% P(chi-square with 2K degrees of freedom > T)
h = T / 2; s = 0; term = 1;
for i = 0:K - 1
    if i > 0, term = term * h / i; end
    s = s + term;
end
q = exp(-h) * s;
end


% ======================================================================
function T = soft_threshold(K, nHyp, target)
% smallest T with nHyp * P(chi2_2K > T) <= target (union bound over taps and drift hypotheses)
lo = 0; hi = 2 * K + 400;
for it = 1:200
    mid = 0.5 * (lo + hi);
    if nHyp * chi2tail(mid, K) > target, lo = mid; else hi = mid; end
end
T = hi;
end


% ======================================================================
function pf = kofn_budget(k, K, target)
% per-fragment probability p with P(Binomial(K, p) >= k) = target
lo = 0; hi = 1;
for it = 1:100
    mid = 0.5 * (lo + hi);
    pk = 0;
    for i = k:K
        pk = pk + nchoosek(K, i) * mid ^ i * (1 - mid) ^ (K - i);
    end
    if pk > target, hi = mid; else lo = mid; end
end
pf = lo;
end
