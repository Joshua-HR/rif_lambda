function export_stageb_csv(root, outDir)
% EXPORT_STAGEB_CSV  Write the stage-B results (current LLS, 2-column dumps) as CSV files for the report plots.
%
%   export_stageb_csv                          % dumps under D:/rif, CSV files into lambda_patch/results/stage_b
%   export_stageb_csv('D:/rif', 'D:/rif/report/csv')
%
%   Re-runs h0_pfa_analyze and h1_pd_analyze with the options of TEST_PLAN.md 3.2 (no figures) on the stage-B
%   dump sets and writes the curves that the checklist notes only summarize. A missing set folder is skipped.
%   The h0 runs reuse h0_stats_cache.mat; the h1 runs re-read every dump and take the longest (8 runs).
%
%   Sets: awgn = h0_awgn / h1_awgn (AWGN, 0.25 ppm), mp = h0_mp / h1_mp (ChanType=RESIDENTIAL LoSPath=1),
%         cfo0 = h0_cfo0 (residual CFO 0). Powers are written in dBm (negative).
%   Detectors: Z = existing Z detector (MD rule), L = CIR Lambda-hat (fixed chip, per fragment).
%
%   Files (one header line each):
%     h1_frag.csv       set,len,dBm,n, Z_k,Z_Pd,Z_lo,Z_hi,Z_PosAcc, L_k,L_Pd,L_lo,L_hi,L_PosAcc
%                       (fragment detection, Clopper-Pearson interval, position accuracy around PeakTap 125)
%     h1_pkt.csv        set,rule,len,dBm,n,k,Pd,lo,hi,PosAcc
%                       rule: soft8 (Lambda-hat soft sum, 8 RIF, 1e-6), kofn6 (6 of 8), single1 / soft2 / soft4
%                       (RIF 1 / 2 / 4 configurations from the 8-fragment dumps)
%     h1_levels.csv     set,rule,len,L90,L99,L90loc   (rule fragZ, fragL or a packet rule above; NaN = not reached)
%     h0_frag.csv       set,len,dBm,n, Z_k,Z_pfa,Z_upper, D_med_dB,kappa_med, L_k,L_pfa,L_upper, L_tail_t10,_t14,_t18
%     h0_frag_tail.csv  set,len,dBm,t,count,nTap,P_meas,P_chi2_2
%                       (Lambda-hat window taps with Lambda >= t, t = 0:0.25:40, against exp(-t/2))
%     h0_pkt.csv        set,K,len,dBm,packets,accepts,rhoPair,pairRatio,t1,t2,t3,ratio1,ratio2,ratio3
%                       (soft packet decision at 1e-6; t1..t3 = grid points where P(chi2_2K >= t) first drops
%                       below 1e-2 / 1e-3 / 1e-4, NaN when fewer than 5 paths are expected there)
%     h0_pkt_tail.csv   set,K,len,dBm,t,count,paths,P_meas,P_chi2_2K
%                       (zero-drift soft-sum paths with S >= t on the integer grid 0:ceil(T))
%
%   See also h0_pfa_analyze, h1_pd_analyze, rif_packet.

if nargin < 1 || isempty(root), root = 'D:/rif'; end
if nargin < 2 || isempty(outDir), outDir = fullfile(fileparts(mfilename('fullpath')), 'results', 'stage_b'); end
if exist(outDir, 'dir') ~= 7, mkdir(outDir); end
pkt = {'Chip', 'fixed', 'Profile', 'fira-medium', 'Drift', 0};

% ---------------------------------------------------------------- H0: fragment (Z, Lambda-hat) and packet tails
H0 = {'awgn', 'h0_awgn'; 'mp', 'h0_mp'; 'cfo0', 'h0_cfo0'};
tg = 0:0.25:40;                                   % the tail grid of h0_pfa_analyze (tail_grid)
fF = open_csv(outDir, 'h0_frag.csv', ['set,len,dBm,n,Z_k,Z_pfa,Z_upper,D_med_dB,kappa_med,' ...
                                     'L_k,L_pfa,L_upper,L_tail_t10,L_tail_t14,L_tail_t18']);
fT = open_csv(outDir, 'h0_frag_tail.csv', 'set,len,dBm,t,count,nTap,P_meas,P_chi2_2');
fP = open_csv(outDir, 'h0_pkt.csv', 'set,K,len,dBm,packets,accepts,rhoPair,pairRatio,t1,t2,t3,ratio1,ratio2,ratio3');
fQ = open_csv(outDir, 'h0_pkt_tail.csv', 'set,K,len,dBm,t,count,paths,P_meas,P_chi2_2K');
for s = 1:size(H0, 1)
    d = fullfile(root, H0{s, 2});
    if exist(d, 'dir') ~= 7, fprintf('skip %s (no folder)\n', d); continue; end
    note(sprintf('%s: fragment Z and Lambda-hat', H0{s, 2}));
    Rc = h0_pfa_analyze(d, 'Suggest', false, 'Moments', 'cir', 'Plot', false);
    C = Rc.cond;
    for c = 1:numel(C)
        fprintf(fF, '%s,%d,%g,%d,%d,%.6g,%.6g,%.4f,%.4f,%d,%.6g,%.6g,%.4f,%.4f,%.4f\n', H0{s, 1}, C(c).len, ...
                -C(c).pow, C(c).n, Rc.E.k(c), Rc.E.pfa(c), Rc.E.upper(c), median(C(c).D), median(C(c).kappa), ...
                Rc.EL.k(c), Rc.EL.pfa(c), Rc.EL.upper(c), Rc.EL.tailRatio(c, :));
        for j = 1:numel(tg)
            fprintf(fT, '%s,%d,%g,%.2f,%d,%d,%.6g,%.6g\n', H0{s, 1}, C(c).len, -C(c).pow, tg(j), ...
                    C(c).tailCnt(j), C(c).nTap, C(c).tailCnt(j) / max(C(c).nTap, 1), exp(-tg(j) / 2));
        end
    end
    Ks = 8;
    if strcmp(H0{s, 1}, 'awgn'), Ks = [2 4 8]; end  % RIF 2 / 4 configurations (B-1R); RIF 1 has no soft sum
    for K = Ks
        note(sprintf('%s: soft packet, %d fragments', H0{s, 2}, K));
        Rp = h0_pfa_analyze(d, 'Suggest', false, pkt{:}, 'FragPerPacket', K, 'Plot', false);
        P = Rp.PK;  g = Rp.pkInfo.tailT;
        for c = 1:numel(P)
            fprintf(fP, '%s,%d,%d,%g,%d,%d,%.5f,%.4f,%g,%g,%g,%.4f,%.4f,%.4f\n', H0{s, 1}, K, P(c).len, -P(c).pow, ...
                    P(c).nPkt, P(c).k, P(c).rhoPair, P(c).pairRatio, P(c).tailT, P(c).tailRatio);
            for j = 1:numel(P(c).tailCnt)
                fprintf(fQ, '%s,%d,%d,%g,%g,%d,%d,%.6g,%.6g\n', H0{s, 1}, K, P(c).len, -P(c).pow, g(j), ...
                        P(c).tailCnt(j), P(c).nPath, P(c).tailCnt(j) / max(P(c).nPath, 1), chi2tail(g(j), K));
            end
        end
    end
end
fclose(fF); fclose(fT); fclose(fP); fclose(fQ);

% ---------------------------------------------------------------- H1: Pd curves and levels
H1 = {'awgn', 'h1_awgn', 3; 'mp', 'h1_mp', 2};      % set, folder, 'Tol' (as in B-H1 / B-H1-3)
fH = open_csv(outDir, 'h1_frag.csv', 'set,len,dBm,n,Z_k,Z_Pd,Z_lo,Z_hi,Z_PosAcc,L_k,L_Pd,L_lo,L_hi,L_PosAcc');
fK = open_csv(outDir, 'h1_pkt.csv', 'set,rule,len,dBm,n,k,Pd,lo,hi,PosAcc');
fL = open_csv(outDir, 'h1_levels.csv', 'set,rule,len,L90,L99,L90loc');
for s = 1:size(H1, 1)
    d = fullfile(root, H1{s, 2});
    if exist(d, 'dir') ~= 7, fprintf('skip %s (no folder)\n', d); continue; end
    pos = {'PeakTap', 125, 'Tol', H1{s, 3}, 'Plot', false};
    note(sprintf('%s: fragment Z and Lambda-hat', H1{s, 2}));
    Rc = h1_pd_analyze(d, 'Moments', 'cir', pos{:});
    C = Rc.cond;
    for c = 1:numel(C)
        fprintf(fH, '%s,%d,%g,%d,%d,%.5f,%.5f,%.5f,%.5f,%d,%.5f,%.5f,%.5f,%.5f\n', H1{s, 1}, C(c).len, C(c).dBm, ...
                C(c).n, C(c).k(1), C(c).Pd(1), C(c).PdLo(1), C(c).PdHi(1), C(c).PosAcc(1), ...
                C(c).k(2), C(c).Pd(2), C(c).PdLo(2), C(c).PdHi(2), C(c).PosAcc(2));
    end
    for q = 1:numel(Rc.len)
        fprintf(fL, '%s,fragZ,%d,%.3f,%.3f,%.3f\n', H1{s, 1}, Rc.len(q), Rc.L90(q, 1), Rc.L99(q, 1), Rc.L90loc(q, 1));
        fprintf(fL, '%s,fragL,%d,%.3f,%.3f,%.3f\n', H1{s, 1}, Rc.len(q), Rc.L90(q, 2), Rc.L99(q, 2), Rc.L90loc(q, 2));
    end
    runs = {'soft8', {}};
    if strcmp(H1{s, 1}, 'awgn')
        runs = {'soft8', {}; 'kofn6', {'Combine', 'kofn', 'KofN', 6}; 'single1', {'FragPerPacket', 1}; ...
                'soft2', {'FragPerPacket', 2}; 'soft4', {'FragPerPacket', 4}};
    end
    for u = 1:size(runs, 1)
        note(sprintf('%s: packet %s', H1{s, 2}, runs{u, 1}));
        Rp = h1_pd_analyze(d, pkt{:}, pos{:}, runs{u, 2}{:});
        PC = Rp.pkt.cond;
        for c = 1:numel(PC)
            fprintf(fK, '%s,%s,%d,%g,%d,%d,%.5f,%.5f,%.5f,%.5f\n', H1{s, 1}, runs{u, 1}, PC(c).len, PC(c).dBm, ...
                    PC(c).n, PC(c).k, PC(c).Pd, PC(c).PdLo, PC(c).PdHi, PC(c).PosAcc);
        end
        for q = 1:numel(Rp.len)
            fprintf(fL, '%s,%s,%d,%.3f,%.3f,%.3f\n', H1{s, 1}, runs{u, 1}, Rp.len(q), Rp.pkt.L90(q), ...
                    Rp.pkt.L99(q), Rp.pkt.L90loc(q));
        end
    end
end
fclose(fH); fclose(fK); fclose(fL);
fprintf('\nCSV files written to %s\n', outDir);
end


% ======================================================================
function fid = open_csv(d, name, header)
fid = fopen(fullfile(d, name), 'w');
if fid < 0, error('export_stageb_csv:open', 'Cannot write %s', fullfile(d, name)); end
fprintf(fid, '%s\n', header);
end


% ======================================================================
function note(msg)
fprintf('\n[%s] %s\n', datestr(now, 'HH:MM:SS'), msg);
end


% ======================================================================
function q = chi2tail(t, K)
% P(chi-square with 2K degrees of freedom >= t)
h = t / 2; s = 0; term = 1;
for i = 0:K - 1
    if i > 0, term = term * h / i; end
    s = s + term;
end
q = exp(-h) * s;
end
