function res = prism_fresnel_probes(case_dir, varargin)
%PRISM_FRESNEL_PROBES  Probe histories of a PRISM run against the analytic Fresnel solution of a plasma slab.
%
%   res = prism_fresnel_probes()                 % case in the current folder
%   res = prism_fresnel_probes(case_dir)
%   res = prism_fresnel_probes(case_dir, 'nu', 1e5, 'floor', 1e-3, 'window', 'auto', 'save', true)
%
%   Reads <case_dir>/input.ini and the probe files written by [probes]: <output_basename>-probe_NN.dat for each
%   point_N and <output_basename>-plane_NN.dat for each plane_N (and the former one-file-per-probe output of the
%   [probe_N] sections, read the same way); the number
%   and positions of the probes come from the input. It then
%     1. plots E and B in time at every probe, with the analytic solution when it exists;
%     2. splits the field into the forward and backward waves, F+- = (E +- s c B)/2, s = propagation verse;
%     3. Fourier transforms the probe signals and compares the spectra with the analytic ones;
%     4. measures R(w) (front probe) and T(w) (front and rear probes) and compares them with the Fresnel
%        (Fabry-Perot) coefficients of a cold plasma slab, eps(w) = 1 - wp^2/(w (w + i nu));
%     5. for a packet with a carrier (gaussian_derivative_pulse with k0 > 0) measures the physics at w0 = c k0 with
%        the same estimator on PRISM and on the linear reference, against the ideal cold plasma: amplitude, phase and
%        group delay of E(w0) against the depth in the slab (slopes -kappa, Re k, 1/v_g), R(w0), T(w0) with their
%        phases and group delays, and the energy fractions of the reflected and transmitted waves.
%   Probes that share their abscissa along the pulse axis form a plane (one probe per transverse cell centre in the
%   Fresnel cases): their mean is the exact k_perp = 0 mode of the grid and replaces them ('planes', true); a lone
%   probe stays as it is.
%
%   References (S. J. Orfanidis, Electromagnetic Waves and Antennas, rev. 2016, freely distributed by the author;
%   chapters 1, 5, 6 saved next to this script): plasma permittivity Eq. (1.15.1); single slab reflection and
%   transmission Eqs. (5.4.3)-(5.4.6), propagation and matching matrices Secs. 5.1-5.2. Orfanidis uses exp(+j w t):
%   his formulas map to the exp(-i w t) ones of this script with j -> -i.
%
%   Analytic solution (exists when [initial_conditions] type is gaussian_pulse or gaussian_derivative_pulse, with or
%   without its carrier k0, the
%   pulse travels along the slab normal and the run is EM, or PIC with plasma_domain = slab): the incident pulse
%   is the input profile, translated at c in vacuum; each frequency is propagated through the slab with the exact
%   normal-incidence solution (transfer matrix), and the time signals are synthesised by inverse FFT: in frequency
%   the solution is a closed formula evaluated exactly, in time it is the exact integral evaluated by quadrature
%   (FFT: sampling, padding and the nu regularisation set its error). The domain is
%   open in the analytic solution: wall reflections (PEC/PML) are not modelled. The dashed vertical line in the
%   time plots is a lower bound of the arrival time of anything reflected by the x walls, with the pulse front at
%   r0 + 5 sigma (the field is below 1e-4 of its peak beyond it); with 'window','auto' the spectra are computed on the
%   signal before that time.
%
%   Plasma: plasma_density is the TOTAL density, electrons and H+ ions n_e = n_i = (1 - neutral_fraction) n/2;
%   wp^2 = n_e e^2/eps0 (1/m_e + 1/m_i). The transverse cold-plasma response holds for a warm plasma too
%   (non-relativistic); particle noise, finite slab edges of the deposit and nonlinear effects are not modelled.
%
%   Options (name/value):
%     'ini'    input file name                                    (default 'input.ini')
%     'nu'     collision frequency [rad/s] in eps(w); a tiny one regularises w = wp   (default 1e-9 wp)
%     'floor'  relative spectral floor: R, T are not shown where |incident| < floor*max   (default 1e-3)
%     'window' 'auto' (cut at the wall-reflection bound), 'none', or [t0 t1] in seconds  (default 'auto')
%     'save'   save the figures as PNG in case_dir                                     (default false)
%     'planes' replace the probes of each plane by their mean                           (default true)
%     'selftest' check the Fresnel solution (energy, Born & Wolf matrix, continuity) and return, no case read
%
%   Returned struct res: probes (time, E, B, F+-, analytic; planes when merged), spectra, R, T, carrier observables
%   (packets), plasma and pulse parameters.

if nargin < 1 || isempty(case_dir), case_dir = pwd; end
opt = parse_options(varargin{:});

% PRISM constants (adam_prism_parameters), so that the analytic solution uses the same numbers as the code
K.EPS0 = 8.854187e-12; K.MU0 = 1.256637e-6; K.C0 = 1/sqrt(K.MU0*K.EPS0);
K.QE = 1.6e-19; K.ME = 9.11e-31; K.MI = 1.67262192595e-27; K.KB = 1.380649e-23;
if opt.selftest, res = selftest(K); return; end

ini = read_ini(fullfile(case_dir, opt.ini));

% ---------------------------------------------------------------- model, units
model = lower(ini_get(ini, 'physics', 'physical_model', 'EM'));
is_adim = any(strcmp(model, {'adimem', 'adimensionalizedelectromagnetic', 'adimensionalelectromagnetic', ...
                              'adimpic', 'adim_pic'}));
is_pic = any(strcmp(model, {'pic', 'particleincell', 'adimpic', 'adim_pic'}));
U.T = 1; U.D = 1; U.B = 1;
if is_adim
   L0 = ini_num(ini, 'physics', 'L0'); B0s = ini_num(ini, 'physics', 'B0');
   U.T = L0/K.C0; U.D = B0s*K.C0*K.EPS0; U.B = B0s;
end

% ---------------------------------------------------------------- pulse
pulse = read_pulse(ini, K);

% ---------------------------------------------------------------- plasma slab
plasma = read_plasma(ini, K, is_pic, pulse);
if plasma.wp > 0 && isempty(opt.nu), opt.nu = 1e-9*plasma.wp; end
if isempty(opt.nu), opt.nu = 0; end

% ---------------------------------------------------------------- probes
basename = strtrim(ini_get(ini, 'IO', 'output_basename', ''));
Pin = probe_inputs(ini, case_dir, basename);             % points, then the lattices of the planes
np = numel(Pin);
if np < 1, error('prism_fresnel_probes: no [probes] in %s', opt.ini); end
ax = pulse.axis; if isempty(ax), ax = 1; end              % sampled components follow the pulse polarization
comp = {{'Dy', 'Bz'}, {'Dz', 'Bx'}, {'Dx', 'By'}};
comp = comp{ax};
names_axis = 'xyz';
dom = [ini_num(ini, 'grid', ['emin_' names_axis(ax)]), ini_num(ini, 'grid', ['emax_' names_axis(ax)])];
s = pulse.verse;                                           % x' = s x: the pulse travels along +x'
P = struct([]);
plane_cache = containers.Map();                            % a plane file is read once
for p = 1:np
   P(p).xyz = Pin(p).xyz;
   P(p).file = Pin(p).file;
   if Pin(p).nq == 0
      [P(p).names, data] = read_probe(P(p).file);
   else
      if ~isKey(plane_cache, P(p).file)
         [nm, d] = read_probe(P(p).file);
         keep = ~ismember(nm, {'i', 'j'});
         nt = floor(size(d, 1)/Pin(p).nq);                   % complete saves only
         plane_cache(P(p).file) = struct('names', {nm(keep)}, 'data', d(1:nt*Pin(p).nq, keep));
      end
      c = plane_cache(P(p).file);
      P(p).names = c.names;
      data = c.data(Pin(p).q:Pin(p).nq:end, :);             % per save, n1*n2 rows, i fastest
   end
   cadence = data(min(2, end), 1) - data(1, 1);            % the last step is saved off cadence: drop it
   if cadence > 0, data = data(mod(data(:, 1) - data(1, 1), cadence) == 0, :); end
   P(p).it = data(:, 1);
   P(p).t = data(:, 2)*U.T;
   P(p).E = column(P(p).names, data, comp{1})*U.D/K.EPS0;
   P(p).B = column(P(p).names, data, comp{2})*U.B;
   P(p).Fp = 0.5*(P(p).E + s*K.C0*P(p).B);                % forward  (+x') wave
   P(p).Fm = 0.5*(P(p).E - s*K.C0*P(p).B);                % backward (-x') wave
   P(p).xa = P(p).xyz(ax);                                 % position along the pulse axis
   P(p).xp = s*P(p).xa;                                    % position along x'
   P(p).where = probe_region(P(p).xp, plasma);
   P(p).n = 1;
   P(p).label = sprintf('probe %d (%s, %.4g m)', p, P(p).where, P(p).xa);
end
dt = check_uniform_time(P);
if opt.planes, P = plane_means(P, s, K.C0); end
np = numel(P);

% ---------------------------------------------------------------- wall-reflection bound and windows
xw = sort(s*dom);                                          % domain walls along x'
for p = 1:np
   P(p).t_wall = wall_time(P(p).xp, pulse, plasma, xw, K.C0);
   P(p).win = time_window(opt.window, P(p).t, P(p).t_wall);
end

% ---------------------------------------------------------------- analytic solution at the probes
has_analytic = pulse.ok && plasma.ok;
if has_analytic
   for p = 1:np
      [P(p).Ea, P(p).Ba] = analytic_probe(P(p).t, dt, P(p).xp, pulse, plasma, opt.nu, K.C0, 4);
      % quadrature check: the same reference with twice the zero padding (wrap of slowly decaying responses)
      E8 = analytic_probe(P(p).t, dt, P(p).xp, pulse, plasma, opt.nu, K.C0, 8);
      P(p).ref_conv = max(abs(P(p).Ea - E8))/max(max(abs(E8)), realmin);
      P(p).Fpa = 0.5*(P(p).Ea + s*K.C0*P(p).Ba);
      P(p).Fma = 0.5*(P(p).Ea - s*K.C0*P(p).Ba);
      in = P(p).win;
      P(p).err_E = norm(P(p).E(in) - P(p).Ea(in))/max(norm(P(p).Ea(in)), realmin);
   end
end

print_summary(P, pulse, plasma, has_analytic, opt, K);

% ---------------------------------------------------------------- figures
figs = {}; tags = {};
figs{end+1} = plot_time(P, has_analytic, comp, 'E'); tags{end+1} = 'time_EB';
figs{end+1} = plot_time(P, has_analytic, comp, 'F'); tags{end+1} = 'time_Fpm';
[S, figs{end+1}] = plot_spectra(P, dt, has_analytic, plasma); tags{end+1} = 'spectra';
[RT, f] = fresnel_coefficients(P, S, plasma, opt, K.C0);
if ~isempty(f), figs{end+1} = f; tags{end+1} = 'fresnel_RT'; end
O = struct();
if has_analytic
   [O, f] = packet_observables(P, S, RT, pulse, plasma, K);
   if ~isempty(f), figs{end+1} = f; tags{end+1} = 'carrier'; end
end

if opt.save
   for k = 1:numel(figs)
      print(figs{k}, fullfile(case_dir, sprintf('%s-probes-%s.png', basename, tags{k})), '-dpng', '-r150');
   end
end

res.probes = P; res.spectra = S; res.RT = RT; res.observables = O; res.pulse = pulse; res.plasma = plasma;
res.options = opt;
end

% ===================================================================================================================
% input
% ===================================================================================================================
function opt = parse_options(varargin)
opt = struct('ini', 'input.ini', 'nu', [], 'floor', 1e-3, 'window', 'auto', 'save', false, 'selftest', false, ...
             'planes', true);
for k = 1:2:numel(varargin)
   key = lower(varargin{k});
   if ~isfield(opt, key), error('prism_fresnel_probes: unknown option "%s"', varargin{k}); end
   opt.(key) = varargin{k+1};
end
end

function ini = read_ini(file)
% INI as a map 'section.key' (lower case) -> value string; ';' and '#' start a comment.
fid = fopen(file, 'r');
if fid < 0, error('prism_fresnel_probes: cannot open %s', file); end
ini = containers.Map();
sec = '';
while true
   line = fgetl(fid);
   if ~ischar(line), break; end
   line = strtok_comment(line);
   if isempty(line), continue; end
   if line(1) == '['
      sec = lower(strtrim(line(2:find(line == ']', 1)-1)));
   else
      eq = find(line == '=', 1);
      if isempty(eq), continue; end
      ini([sec '.' lower(strtrim(line(1:eq-1)))]) = strtrim(line(eq+1:end));
   end
end
fclose(fid);
end

function line = strtok_comment(line)
cut = find(line == ';' | line == '#', 1);
if ~isempty(cut), line = line(1:cut-1); end
line = strtrim(strrep(line, char(13), ''));
end

function v = ini_get(ini, sec, key, default)
k = [lower(sec) '.' lower(key)];
if isKey(ini, k)
   v = ini(k);
elseif nargin > 3
   v = default;
else
   error('prism_fresnel_probes: missing [%s] %s', sec, key);
end
end

function v = ini_num(ini, sec, key, default)
if nargin > 3
   v = ini_get(ini, sec, key, num2str(default));
else
   v = ini_get(ini, sec, key);
end
v = str2double(v);
end

function pulse = read_pulse(ini, K)
pulse = struct('ok', false, 'type', lower(ini_get(ini, 'initial_conditions', 'type', 'vacuum')), ...
               'axis', [], 'verse', 1, 'B0', 0, 'r0', 0, 'sigma', 0, 'k0', 0, 'length', 0, 'extent', 0, ...
               'c', K.C0);
if ~any(strcmp(pulse.type, {'gaussian_pulse', 'gaussian_derivative_pulse'}))
   warning('prism_fresnel_probes: IC "%s" has no analytic pulse: numerical results only', pulse.type);
   return
end
pdir = strtrim(ini_get(ini, 'initial_conditions', 'pulse_direction'));
pulse.axis = find('xyz' == lower(pdir(end)));
pulse.verse = 1; if pdir(1) == '-', pulse.verse = -1; end
pulse.B0 = ini_num(ini, 'initial_conditions', 'B0');
pulse.r0 = ini_num(ini, 'initial_conditions', 'r0');
pulse.sigma = ini_num(ini, 'initial_conditions', 'sigma');
if strcmp(pulse.type, 'gaussian_pulse')
   pulse.k0 = ini_num(ini, 'initial_conditions', 'k0');
   pulse.length = 5*pulse.sigma;                          % the code cuts the profile at 5 sigma
else
   pulse.k0 = ini_num(ini, 'initial_conditions', 'k0', 0); % optional carrier, 0 by default (single cycle)
   pulse.length = 8*pulse.sigma;                          % no cut: 8 sigma leaves e^-32 (support of the quadrature)
end
pulse.extent = 5*pulse.sigma;                              % front for the wall bound: |field| < 1e-4 of the peak
pulse.ok = true;
end

function E = pulse_profile(sp, pulse)
% Incident E (the y, z or x component of the pulse polarization) at the abscissa sp = x' - x0' - c t along x'.
% Same formulas as adam_prism_ic_object: the IC coordinate is s = x - r0 = verse*sp.
sv = pulse.verse*sp;
if strcmp(pulse.type, 'gaussian_pulse')
   E = pulse.B0*pulse.c*exp(-0.5*(sv/pulse.sigma).^2).*cos(pulse.k0*sv);
   E(abs(sv/pulse.sigma) > 5) = 0;
else
   % A = A0 exp(-s^2/2 sigma^2) cos(k0 s), E = c dA/dx; A0/sigma = B0 e^(1/2) (k0 = 0), B0/sqrt((k0 sigma)^2 + 1/e)
   u = sv/pulse.sigma; ks = pulse.k0*pulse.sigma;
   if ks == 0, a = exp(0.5); else, a = 1/sqrt(ks^2 + exp(-1)); end
   E = -pulse.B0*pulse.c*a*(u.*cos(ks*u) + ks*sin(ks*u)).*exp(-0.5*u.^2);
end
end

function plasma = read_plasma(ini, K, is_pic, pulse)
plasma = struct('ok', pulse.ok, 'wp', 0, 'wpe', 0, 'ne', 0, 'a', 0, 'b', 0, 'slab', false, 'vth', 0);
if ~is_pic, return; end
if ~strcmpi(strtrim(ini_get(ini, 'PIC', 'problem_type', 'plasma')), 'plasma'), return; end
domain = lower(strtrim(ini_get(ini, 'PIC', 'plasma_domain', '')));
n = ini_num(ini, 'PIC', 'plasma_density');
f = ini_num(ini, 'PIC', 'neutral_fraction', 0);
plasma.ne = (1 - f)*n/2;
plasma.wpe = sqrt(plasma.ne*K.QE^2/(K.EPS0*K.ME));
plasma.wp = sqrt(plasma.ne*K.QE^2/K.EPS0*(1/K.ME + 1/K.MI));
plasma.vth = sqrt(K.KB*ini_num(ini, 'particle_injection', 'electronic_temperature', 0)/K.ME); % as PRISM's VTH
if ~strcmp(domain, 'slab')
   warning('prism_fresnel_probes: plasma_domain "%s" is not a slab: no Fresnel solution', domain);
   plasma.ok = false; return
end
normal = find('xyz' == lower(strtrim(ini_get(ini, 'PIC', 'slab_normal'))));
if pulse.ok && normal ~= pulse.axis
   warning('prism_fresnel_probes: slab normal and pulse direction differ: no normal-incidence solution');
   plasma.ok = false; return
end
c = ini_num(ini, 'PIC', 'slab_center'); w = ini_num(ini, 'PIC', 'slab_thickness');
ab = sort(pulse.verse*[c - w/2, c + w/2]);                 % slab faces along x'
plasma.a = ab(1); plasma.b = ab(2); plasma.slab = true;
end

function Pin = probe_inputs(ini, case_dir, basename)
% Probe points and files in PRISM's global order: point_N (one file each), then the plane_N lattices (one file per
% plane: per save n1*n2 rows, i fastest, over the cell centres of a uniform lattice spanning the domain along the two
% axes of the plane). Former input: probes_number and one [probe_N] section per probe, one file per probe.
Pin = struct('xyz', {}, 'file', {}, 'q', {}, 'nq', {});
if isKey(ini, 'probes.probes_number')
   n = round(ini_num(ini, 'probes', 'probes_number'));
   fmt = sprintf('%%s-probe_%%0%dd.dat', max(2, numel(num2str(n))));
   for p = 1:n
      sec = sprintf('probe_%d', p);
      Pin(p) = struct('xyz', [ini_num(ini, sec, 'x'), ini_num(ini, sec, 'y'), ini_num(ini, sec, 'z')], ...
                      'file', fullfile(case_dir, sprintf(fmt, basename, p)), 'q', 0, 'nq', 0);
   end
   return
end
npt = 0; while isKey(ini, sprintf('probes.point_%d', npt+1)), npt = npt + 1; end
npl = 0; while isKey(ini, sprintf('probes.plane_%d', npl+1)), npl = npl + 1; end
fmt = sprintf('%%s-probe_%%0%dd.dat', max(2, numel(num2str(npt))));
for p = 1:npt
   v = sscanf(ini(sprintf('probes.point_%d', p)), '%f');
   Pin(end+1) = struct('xyz', v(1:3).', 'file', fullfile(case_dir, sprintf(fmt, basename, p)), 'q', 0, 'nq', 0); %#ok<AGROW>
end
lo = [ini_num(ini, 'grid', 'emin_x'), ini_num(ini, 'grid', 'emin_y'), ini_num(ini, 'grid', 'emin_z')];
hi = [ini_num(ini, 'grid', 'emax_x'), ini_num(ini, 'grid', 'emax_y'), ini_num(ini, 'grid', 'emax_z')];
fmt = sprintf('%%s-plane_%%0%dd.dat', max(2, numel(num2str(npl))));
for l = 1:npl
   tok = strsplit(strtrim(ini(sprintf('probes.plane_%d', l))));
   a = find('xyz' == lower(strrep(strrep(tok{1}, '"', ''), '''', '')), 1);
   if isempty(a), error('prism_fresnel_probes: [probes] plane_%d axis "%s" is not x, y or z', l, tok{1}); end
   ax = [a, setdiff(1:3, a)];                             % normal, first, second lattice axis
   x = str2double(tok{2}); n1 = str2double(tok{3}); n2 = str2double(tok{4});
   file = fullfile(case_dir, sprintf(fmt, basename, l));
   q = 0;
   for j = 1:n2
      for i = 1:n1
         q = q + 1;
         xyz = zeros(1, 3);
         xyz(ax(1)) = x;
         xyz(ax(2)) = lo(ax(2)) + (i - 0.5)*(hi(ax(2)) - lo(ax(2)))/n1;
         xyz(ax(3)) = lo(ax(3)) + (j - 0.5)*(hi(ax(3)) - lo(ax(3)))/n2;
         Pin(end+1) = struct('xyz', xyz, 'file', file, 'q', q, 'nq', n1*n2); %#ok<AGROW>
      end
   end
end
% a finished run of the former output (one file per probe, numbered with the digits of the probe count) whose input
% was converted afterwards: read the per-probe files
fmt = sprintf('%%s-probe_%%0%dd.dat', max(2, numel(num2str(numel(Pin)))));
if npl > 0 && ~isfile(fullfile(case_dir, sprintf('%s-plane_%02d.dat', basename, 1))) && ...
   isfile(fullfile(case_dir, sprintf(fmt, basename, 1)))
   for p = 1:numel(Pin)
      Pin(p).file = fullfile(case_dir, sprintf(fmt, basename, p));
      Pin(p).q = 0;
      Pin(p).nq = 0;
   end
end
end

function [names, data] = read_probe(file)
fid = fopen(file, 'r');
if fid < 0, error('prism_fresnel_probes: cannot open %s', file); end
names = {};
while isempty(names)
   line = fgetl(fid);
   if ~ischar(line), error('prism_fresnel_probes: no %%VARIABLES header in %s', file); end
   if strncmp(line, '%VARIABLES', 10)
      names = regexp(line, '"([^"]*)"', 'tokens');
      names = cellfun(@(c) c{1}, names, 'UniformOutput', false);
   end
end
C = textscan(fid, repmat('%f', 1, numel(names)), 'CommentStyle', '%', 'CollectOutput', true);
fclose(fid);
data = C{1};
end

function v = column(names, data, name)
k = find(strcmp(names, name), 1);
if isempty(k), error('prism_fresnel_probes: variable %s not in the probe file', name); end
v = data(:, k);
end

function where = probe_region(xp, plasma)
if ~plasma.slab
   where = 'vacuum';
elseif xp < plasma.a
   where = 'front';
elseif xp > plasma.b
   where = 'rear';
else
   where = 'inside';
end
end

function Q = plane_means(P, s, c)
% Transverse means of the probes that share their abscissa along the pulse axis: a y-z plane of probes, one per
% transverse cell centre, averages to the exact k_perp = 0 mode of the grid. Planes are ordered along x'.
xa = [P.xa];
tol = 1e-9*(max(abs(xa)) + realmin);
[~, ~, g] = unique(round(xa/tol));
Q = P([]);
for k = 1:max(g)
   idx = find(g == k);
   q = P(idx(1));
   q.n = numel(idx);
   q.E = mean([P(idx).E], 2);
   q.B = mean([P(idx).B], 2);
   q.Fp = 0.5*(q.E + s*c*q.B);
   q.Fm = 0.5*(q.E - s*c*q.B);
   if q.n > 1, q.label = sprintf('plane at %.4g m (%s, %d probes)', q.xa, q.where, q.n); end
   Q(end+1) = q; %#ok<AGROW>
end
[~, o] = sort([Q.xp]);
Q = Q(o);
end

function dt = check_uniform_time(P)
dt = mean(diff(P(1).t));
for p = 1:numel(P)
   d = diff(P(p).t);
   if any(abs(d - dt) > 1e-6*dt) || numel(P(p).t) ~= numel(P(1).t)
      error('prism_fresnel_probes: probe %d is not sampled on the uniform time grid of probe 1', p);
   end
end
end

function tw = wall_time(xp, pulse, plasma, xw, c)
% Lower bound of the arrival at xp of a wave reflected by a domain wall: the pulse front starts at x0' + extent and
% nothing travels faster than c (the slab only delays).
if ~pulse.ok, tw = Inf; return; end
front = pulse.verse*pulse.r0 + pulse.extent;
right = (xw(2) - front) + (xw(2) - xp);                    % on to the far wall and back
if plasma.slab
   left = (plasma.a - front) + (plasma.a - xw(1)) + (xp - xw(1)); % reflected by the slab, then by the near wall
else
   left = Inf;
end
tw = min(right, left)/c;
end

function in = time_window(w, t, t_wall)
if ischar(w) && strcmpi(w, 'none')
   in = true(size(t));
elseif ischar(w)                                           % 'auto'
   in = t < t_wall;
else
   in = t >= w(1) & t <= w(2);
end
end

% ===================================================================================================================
% analytic solution
% ===================================================================================================================
function [E, B] = analytic_probe(t, dt, xp, pulse, plasma, nu, c, pad)
% E, B at x' = xp for the times t: incident spectrum at the slab front x' = a, times the exact transfer function.
a = plasma.a; if ~plasma.slab, a = xp; end
x0 = pulse.verse*pulse.r0;
t_lo = min(0, (a - x0 - pulse.length)/c - 4*dt);
m = ceil(-t_lo/dt);                                        % t_lo = -m dt: the grid contains the probe times
span = (t(end) + m*dt);
N = 2^nextpow2(pad*ceil(span/dt));                         % zero padding: no wrap of slowly decaying responses
tg = (-m + (0:N-1)')*dt;
g = pulse_profile(a - x0 - c*tg, pulse);                   % incident field at x' = a
G = N*ifft(g);                                             % spectrum, convention exp(-i w t)
w = 2*pi/(N*dt)*[0:N/2, -N/2+1:-1]';
[HE, HB] = transfer(abs(w), xp, plasma, nu, c, pulse.verse);
neg = w < 0;
HE(neg) = conj(HE(neg)); HB(neg) = conj(HB(neg));          % real signals: H(-w) = conj(H(w))
HE(N/2+1) = real(HE(N/2+1)); HB(N/2+1) = real(HB(N/2+1));
Eg = real(fft(G.*HE))/N;
Bg = real(fft(G.*HB))/N;
idx = round(t/dt) + m + 1;
E = Eg(idx); B = Bg(idx);
end

function [HE, HB, r, tt] = transfer(w, xp, plasma, nu, c, verse)
% Exact normal-incidence field of a unit incident wave exp(i k (x' - a)) on the slab [a, b] (vacuum around it),
% convention exp(-i w t); B is the component paired with E by the IC (B = verse E/c for the forward vacuum wave).
k = w/c;
if ~plasma.slab
   HE = exp(1i*k*0); HB = verse/c*HE; r = 0*k; tt = 1 + 0*k;
   return
end
a = plasma.a; d = plasma.b - plasma.a;
wsafe = max(w, realmin);
epsr = 1 - plasma.wp^2./(wsafe.*(wsafe + 1i*nu));
n = sqrt(epsr);                                             % Im n >= 0: evanescent/absorbed into the slab
r12 = (1 - n)./(1 + n); r23 = -r12; t12 = 2./(1 + n);
e2 = exp(2i*k.*n*d);
den = 1 - r12.^2.*e2;
C = t12./den;
r = (r12 + r23.*e2)./den;
tt = C.*(1 + r23).*exp(1i*k.*n*d);
zero = (w == 0);                                            % w -> 0: eps -> -inf, total reflection
r(zero) = -1; tt(zero) = 0; C(zero) = 0;
if xp < a
   ef = exp(1i*k*(xp - a)); eb = r.*exp(-1i*k*(xp - a));
   HE = ef + eb; HB = verse/c*(ef - eb);
elseif xp > plasma.b
   ef = tt.*exp(1i*k*(xp - plasma.b));
   HE = ef; HB = verse/c*ef;
else
   u = xp - a;
   ef = C.*exp(1i*k.*n*u);
   eb = C.*r23.*exp(1i*k.*n*(2*d - u));                     % D exp(-i k n u), written without overflow
   HE = ef + eb; HB = verse/c*n.*(ef - eb);
   HE(zero) = 0; HB(zero) = 0;
end
end

% ===================================================================================================================
% spectra and Fresnel coefficients
% ===================================================================================================================
function [S, fig] = plot_spectra(P, dt, has_analytic, plasma)
fig = figure('Name', 'Spectra at the probes', 'Color', 'w');
np = numel(P);
N = 2^nextpow2(4*numel(P(1).t));                          % one length for all: common frequency grid
for p = 1:np
   in = P(p).win;
   S(p).w = 2*pi/(N*dt)*(0:N/2)';
   S(p).Fp = spectrum(P(p).Fp(in), N, dt);
   S(p).Fm = spectrum(P(p).Fm(in), N, dt);
   S(p).E = spectrum(P(p).E(in), N, dt);
   if has_analytic
      S(p).Ea = spectrum(P(p).Ea(in), N, dt);
   end
   subplot(np, 1, p);
   [xw, xl] = freq_axis(S(p).w, plasma);
   semilogy(xw, abs(S(p).E), 'b-', 'LineWidth', 1.2); hold on;
   if has_analytic, semilogy(xw, abs(S(p).Ea), 'k--', 'LineWidth', 1.0); end
   ylim_floor(abs(S(p).E));
   cutoff_line(plasma);
   grid on; ylabel('|E(\omega)|');
   title(P(p).label);
   if has_analytic, legend('PRISM', 'linear ref. (Fresnel + IFFT)', 'Location', 'best'); end
   xlim([0 xw_max(xw, plasma)]);
end
xlabel(xl);
end

function X = spectrum(x, N, dt)
X = dt*N*ifft(x, N);                                       % convention exp(-i w t), as the analytic solution
X = X(1:N/2+1);
end

function [RT, fig] = fresnel_coefficients(P, S, plasma, opt, c)
RT = struct(); fig = [];
if ~plasma.slab, return; end
front = find(strcmp({P.where}, 'front')); rear = find(strcmp({P.where}, 'rear'));
if isempty(front)
   warning('prism_fresnel_probes: no probe in front of the slab: R(w) and T(w) are not measured');
   return
end
[~, k] = max([P(front).xp]); f = front(k);                 % the probes closest to the slab
w = S(f).w; kv = w/c;
inc = S(f).Fp;
ok = abs(inc) > opt.floor*max(abs(inc));
RT.w = w;
RT.R = S(f).Fm./inc.*exp(-2i*kv*(plasma.a - P(f).xp));
[~, ~, RT.Ra, RT.Ta] = transfer(w, plasma.a - 1, plasma, opt.nu, c, 1);
RT.T = nan(size(w));
if ~isempty(rear)
   [~, k] = min([P(rear).xp]); b = rear(k);
   RT.T = S(b).Fp./inc.*exp(-1i*kv*(P(b).xp - plasma.b)).*exp(1i*kv*(P(f).xp - plasma.a));
end
RT.R(~ok) = NaN; RT.T(~ok) = NaN;
RT.ok = ok;

fig = figure('Name', 'Fresnel coefficients', 'Color', 'w');
[xw, xl] = freq_axis(w, plasma);
subplot(2, 2, 1);
plot(xw, abs(RT.R), 'b-', xw, abs(RT.Ra), 'k--', 'LineWidth', 1.2); grid on; cutoff_line(plasma);
ylabel('|R(\omega)|'); legend('PRISM', 'Fresnel', 'Location', 'best'); ylim([0 1.1]);
xlim([0 xw_max(xw, plasma)]);
subplot(2, 2, 2);
plot(xw, abs(RT.T), 'r-', xw, abs(RT.Ta), 'k--', 'LineWidth', 1.2); grid on; cutoff_line(plasma);
ylabel('|T(\omega)|'); legend('PRISM', 'Fresnel', 'Location', 'best'); ylim([0 1.1]);
xlim([0 xw_max(xw, plasma)]);
subplot(2, 2, 3);
plot(xw, abs(RT.R).^2 + abs(RT.T).^2, 'm-', xw, abs(RT.Ra).^2 + abs(RT.Ta).^2, 'k--', 'LineWidth', 1.2);
grid on; cutoff_line(plasma); ylabel('|R|^2 + |T|^2'); xlabel(xl); ylim([0 1.2]);
xlim([0 xw_max(xw, plasma)]);
subplot(2, 2, 4);
plot(xw, angle(RT.R), 'b-', xw, angle(RT.Ra), 'k--', 'LineWidth', 1.2); grid on; cutoff_line(plasma);
ylabel('arg R(\omega) [rad]'); xlabel(xl); xlim([0 xw_max(xw, plasma)]);
end

function [O, fig] = packet_observables(P, S, RT, pulse, plasma, K)
% Narrow-band packet (k0 > 0): the physics at the carrier w0 = c k0, the same estimator on PRISM and on the linear
% reference (spectra of the windowed signals, exp(-i w t)), against the ideal cold plasma n = sqrt(1 - wp^2/w0^2).
%   inside the slab, against the depth x' - a: log |E(w0)| (slope -kappa), arg E(w0) (slope Re k, phase velocity
%   w0/Re k) and the group delay d arg E/dw at w0 (slope 1/v_g, v_g = c sqrt(1 - wp^2/w0^2) above the cutoff);
%   faces: R(w0), T(w0) with phases and group delays referred to the faces (R ones omitted above the cutoff, where
%   the cases sit on a Fabry-Perot null of R), energy fractions sum F-^2/sum F+^2 (front), sum F+_rear^2/sum F+^2.
%   Ideal: below the cutoff the half-space reflection phase -2 atan(kappa c/w0) and delay 2/(c kappa); above it one
%   pass through the slab, phase Re k d and delay d/v_g (Orfanidis EWA eq. 1.15.2, sec. 5.4).
O = struct(); fig = [];
if pulse.k0 <= 0 || ~plasma.slab || plasma.wp <= 0, return; end
c = K.C0; w0 = c*pulse.k0; half = 0.5*c/pulse.sigma;
D = c/plasma.wpe; T0 = 1/plasma.wpe; E0 = pulse.B0*c;
front = find(strcmp({P.where}, 'front')); rear = find(strcmp({P.where}, 'rear'));
inside = find(strcmp({P.where}, 'inside'));
if isempty(front) || isempty(rear) || numel(inside) < 2
   warning('prism_fresnel_probes: carrier observables need a front, a rear and two inside probes or planes');
   return
end
[~, k] = max([P(front).xp]); f = front(k);
[~, k] = min([P(rear).xp]); r = rear(k);
[~, o] = sort([P(inside).xp]); inside = inside(o);
depth = ([P(inside).xp]' - plasma.a)/D;
w = S(1).w; above = w0 > plasma.wp;
n0 = sqrt(1 - plasma.wp^2/w0^2 + 0i); kid = n0*w0/c;
O.ideal.kappa_delta = imag(kid)*D; O.ideal.Re_k_delta = real(kid)*D;
O.ideal.v_phase_over_c = NaN; O.ideal.v_group_over_c = NaN;
O.ideal.R_phase = NaN; O.ideal.R_delay_wpe = NaN; O.ideal.T_phase = NaN; O.ideal.T_delay_wpe = NaN;
if above
   O.ideal.v_phase_over_c = w0/(real(kid)*c); O.ideal.v_group_over_c = sqrt(1 - plasma.wp^2/w0^2);
   O.ideal.T_phase = angle(exp(1i*real(kid)*(plasma.b - plasma.a)));
   O.ideal.T_delay_wpe = (plasma.b - plasma.a)/(O.ideal.v_group_over_c*c)/T0;
else
   O.ideal.R_phase = -2*atan(imag(kid)*c/w0); O.ideal.R_delay_wpe = 2/(c*imag(kid))/T0;
end
tags = {'prism', 'ref'}; spec = {'E', 'Ea'}; Rn = {'R', 'Ra'}; Tn = {'T', 'Ta'};
Fp = {'Fp', 'Fpa'}; Fm = {'Fm', 'Fma'};
for j = 1:2
   A = zeros(numel(inside), 1); ph = A; tau = A;
   for m = 1:numel(inside)
      [A(m), ph(m), tau(m)] = carrier_fit(w, S(inside(m)).(spec{j}), w0, half);
   end
   ph = unwrap(ph);
   o = struct('depth_delta', depth, 'amp', A, 'phase', ph, 'group_delay_wpe', tau/T0);
   pa = polyfit(depth, log(A), 1); pp = polyfit(depth, ph, 1); pt = polyfit(depth, tau/T0, 1);
   o.kappa_delta = -pa(1); o.Re_k_delta = pp(1); o.inv_vg_c = pt(1);
   o.v_phase_over_c = NaN; o.v_group_over_c = NaN;
   if above && o.Re_k_delta > 0, o.v_phase_over_c = w0/plasma.wpe/o.Re_k_delta; end
   if above && o.inv_vg_c > 0, o.v_group_over_c = 1/o.inv_vg_c; end
   [o.R_abs, o.R_phase, o.R_delay_wpe] = carrier_fit(w, RT.(Rn{j}), w0, half);
   [o.T_abs, o.T_phase, o.T_delay_wpe] = carrier_fit(w, RT.(Tn{j}), w0, half);
   o.R_delay_wpe = o.R_delay_wpe/T0; o.T_delay_wpe = o.T_delay_wpe/T0;
   if above, o.R_phase = NaN; o.R_delay_wpe = NaN; end
   inc = sum(P(f).(Fp{j})(P(f).win).^2);
   o.energy_R = sum(P(f).(Fm{j})(P(f).win).^2)/inc;
   o.energy_T = sum(P(r).(Fp{j})(P(r).win).^2)/inc;
   O.(tags{j}) = o;
end

fprintf('   carrier observables at w0 = %.4g w_p (PRISM, linear reference, ideal cold plasma):\n', w0/plasma.wp);
rows = {'kappa c/w_pe', 'kappa_delta'; 'Re k c/w_pe', 'Re_k_delta'; 'v_phase/c', 'v_phase_over_c'; ...
        'v_group/c', 'v_group_over_c'; '|R(w0)|', 'R_abs'; 'arg R(w0)', 'R_phase'; 'R delay w_pe', 'R_delay_wpe'; ...
        '|T(w0)|', 'T_abs'; 'arg T(w0)', 'T_phase'; 'T delay w_pe (transit)', 'T_delay_wpe'; ...
        'energy R', 'energy_R'; 'energy T', 'energy_T'};
idm = struct('kappa_delta', O.ideal.kappa_delta, 'Re_k_delta', O.ideal.Re_k_delta, ...
             'v_phase_over_c', O.ideal.v_phase_over_c, 'v_group_over_c', O.ideal.v_group_over_c, ...
             'R_phase', O.ideal.R_phase, 'R_delay_wpe', O.ideal.R_delay_wpe, 'T_phase', O.ideal.T_phase, ...
             'T_delay_wpe', O.ideal.T_delay_wpe);
for m = 1:size(rows, 1)
   key = rows{m, 2}; idv = NaN;
   if isfield(idm, key), idv = idm.(key); end
   fprintf('     %-24s %12.5g %12.5g %12.5g\n', rows{m, 1}, O.prism.(key), O.ref.(key), idv);
end
if above
   fprintf('     (ideal T: one pass at the cold-plasma k and v_g, no Fabry-Perot)\n');
else
   fprintf('     (ideal R: semi-infinite plasma; T delay: tunnelling time through the slab)\n');
end

fig = figure('Name', 'Carrier observables', 'Color', 'w');
subplot(2, 2, 1);
semilogy(depth, O.prism.amp, 'bo', depth, O.ref.amp, 'k--', ...
         depth, O.ref.amp(1)*exp(-O.ideal.kappa_delta*(depth - depth(1))), 'r:', 'LineWidth', 1.2);
grid on; xlabel('depth (x - a) \omega_{pe}/c'); ylabel('|E(\omega_0)|');
legend('PRISM', 'linear ref.', 'ideal e^{-\kappa x}', 'Location', 'best');
subplot(2, 2, 2);
plot(depth, O.prism.phase, 'bo', depth, O.ref.phase, 'k--', ...
     depth, O.ref.phase(1) + O.ideal.Re_k_delta*(depth - depth(1)), 'r:', 'LineWidth', 1.2);
grid on; xlabel('depth (x - a) \omega_{pe}/c'); ylabel('arg E(\omega_0) [rad]');
legend('PRISM', 'linear ref.', 'ideal Re k', 'Location', 'best');
subplot(2, 2, 3);
plot(depth, O.prism.group_delay_wpe, 'bo', depth, O.ref.group_delay_wpe, 'k--', 'LineWidth', 1.2); hold on;
if above
   plot(depth, O.ref.group_delay_wpe(1) + (depth - depth(1))/O.ideal.v_group_over_c, 'r:', 'LineWidth', 1.2);
   legend('PRISM', 'linear ref.', 'ideal 1/v_g', 'Location', 'best');
else
   legend('PRISM', 'linear ref.', 'Location', 'best');
end
grid on; xlabel('depth (x - a) \omega_{pe}/c'); ylabel('group delay \omega_{pe}\tau');
subplot(2, 2, 4);
band = abs(w - w0) <= 3*half; xw = w(band)/plasma.wp;
plot(xw, abs(RT.R(band)).^2, 'b-', xw, abs(RT.Ra(band)).^2, 'b--', ...
     xw, abs(RT.T(band)).^2, 'r-', xw, abs(RT.Ta(band)).^2, 'r--', 'LineWidth', 1.1); hold on;
yl = [-0.05 1.1]; plot(w0/plasma.wp*[1 1], yl, 'g:'); ylim(yl);
grid on; xlabel('\omega/\omega_p'); legend('|R|^2', '|R|^2 ref.', '|T|^2', '|T|^2 ref.', 'Location', 'best');
end

function [A, ph, tau] = carrier_fit(w, X, w0, half)
% |X(w0)|, arg X(w0) and the group delay d arg X/dw at w0: |X|^2-weighted linear fit of the unwrapped phase over
% w0 +- half (NaN samples, below the spectral floor, are skipped).
sel = abs(w - w0) <= half & isfinite(X);
ws = w(sel) - w0; x = X(sel);
ph = unwrap(angle(x)); wt = abs(x).^2;
M = [ones(size(ws)), ws];
cf = (M'*(M.*wt))\(M'*(ph.*wt));
A = interp1(w(isfinite(X)), abs(X(isfinite(X))), w0);
ph = cf(1); tau = cf(2);
end

% ===================================================================================================================
% plots and summary
% ===================================================================================================================
function fig = plot_time(P, has_analytic, comp, what)
np = numel(P);
if what == 'E'
   fig = figure('Name', 'Fields at the probes', 'Color', 'w');
   lab = {sprintf('E_%s [V/m]', comp{1}(2)), sprintf('B_%s [T]', comp{2}(2))};
   num = {'E', 'B'}; ana = {'Ea', 'Ba'};
else
   fig = figure('Name', 'Forward/backward waves at the probes', 'Color', 'w');
   lab = {'F^+ [V/m]', 'F^- [V/m]'};
   num = {'Fp', 'Fm'}; ana = {'Fpa', 'Fma'};
end
for p = 1:np
   for j = 1:2
      subplot(np, 2, 2*(p-1)+j);
      tn = P(p).t*1e9;
      [sc, sl] = axis_scale(P(p).(num{j}));               % the scale goes in the label, not over the title
      plot(tn, P(p).(num{j})/sc, 'b-', 'LineWidth', 1.2); hold on;
      if has_analytic, plot(tn, P(p).(ana{j})/sc, 'k--', 'LineWidth', 1.0); end
      if isfinite(P(p).t_wall) && P(p).t_wall < P(p).t(end)
         yl = ylim; plot(P(p).t_wall*1e9*[1 1], yl, 'r:', 'LineWidth', 1.0); ylim(yl);
      end
      grid on; ylabel(strrep(lab{j}, '[', ['[' sl]));
      if j == 1
         t = P(p).label;
         if has_analytic && what == 'E', t = sprintf('%s, err_{L2}(E) = %.2e', t, P(p).err_E); end
         title(t);
      end
      if p == 1 && j == 2 && has_analytic, legend('PRISM', 'linear ref. (Fresnel + IFFT)', 'Location', 'best'); end
      if p == np, xlabel('t [ns]'); end
   end
end
end

function [sc, sl] = axis_scale(y)
e = floor(log10(max(abs(y)) + realmin)/3)*3;
sc = 10^e; sl = '';
if e ~= 0, sl = sprintf('10^{%d} ', e); end
end

function [xw, xl] = freq_axis(w, plasma)
if plasma.wp > 0
   xw = w/plasma.wp; xl = '\omega/\omega_p';
else
   xw = w/(2*pi)*1e-9; xl = 'f [GHz]';
end
end

function m = xw_max(xw, plasma)
if plasma.wp > 0, m = min(xw(end), 4); else, m = xw(end); end
end

function cutoff_line(plasma)
if plasma.wp > 0
   yl = ylim; hold on; plot([1 1], yl, 'Color', [0.5 0.5 0.5], 'LineStyle', ':'); ylim(yl);
end
end

function ylim_floor(y)
m = max(y);
if m > 0, ylim([m*1e-6, m*2]); end
end

function print_summary(P, pulse, plasma, has_analytic, opt, K)
fprintf('\n== PRISM probes: %d probes\n', sum([P.n]));
if pulse.ok
   sgn = '+'; if pulse.verse < 0, sgn = '-'; end
   axn = 'xyz';
   fprintf('   pulse %s along %s%s: B0 = %g T, r0 = %g m, sigma = %g m', pulse.type, sgn, axn(pulse.axis), ...
           pulse.B0, pulse.r0, pulse.sigma);
   if pulse.k0 ~= 0, fprintf(', k0 = %g 1/m', pulse.k0); end
   fprintf('\n');
end
if plasma.wp > 0
   fprintf('   plasma: n_e = %g m^-3, w_p = %g rad/s (f_p = %g Hz), c/w_p = %g m', plasma.ne, plasma.wp, ...
           plasma.wp/(2*pi), K.C0/plasma.wp);
   if plasma.slab
      fprintf(', slab [%g, %g] m along x'' (%g c/w_p)', plasma.a, plasma.b, (plasma.b - plasma.a)*plasma.wp/K.C0);
   end
   fprintf(', nu = %g rad/s\n', opt.nu);
   if pulse.ok && pulse.k0 > 0
      fprintf('   carrier w0 = c k0 = %g w_p, k0 sigma = %g (relative spectral width ~ 1/(k0 sigma))\n', ...
              K.C0*pulse.k0/plasma.wp, pulse.k0*pulse.sigma);
   elseif pulse.ok
      fprintf('   sigma w_p/c = %g (spectral peak of the derivative pulse at w = c/sigma)\n', pulse.sigma*plasma.wp/K.C0);
   end
end
for p = 1:numel(P)
   if P(p).n > 1
      fprintf('   plane %d: x = %g m (mean of %d probes), %s, %d samples, window t < %g ns', p, P(p).xa, P(p).n, ...
              P(p).where, numel(P(p).t), max(P(p).t(P(p).win))*1e9);
   else
      fprintf('   probe %d: (%g, %g, %g) m, %s, %d samples, window t < %g ns', p, P(p).xyz, P(p).where, ...
              numel(P(p).t), max(P(p).t(P(p).win))*1e9);
   end
   if has_analytic
      fprintf(', L2 error of E vs linear reference = %.3e (reference quadrature, padding 4x vs 8x: %.1e)', ...
              P(p).err_E, P(p).ref_conv);
   end
   fprintf('\n');
end
if ~has_analytic, fprintf('   no linear reference for this case\n'); end
print_validity(P, pulse, plasma, K);
end

function print_validity(P, pulse, plasma, K)
% Hypotheses of the linear cold reference: warm corrections O((v_th/c)^2), slab edges sharp over the run, linear
% response (a0 = e E/(m_e w c) << 1 at the spectral peak of the pulse). Thresholds are indicative.
if plasma.wp <= 0 || ~pulse.ok, return; end
fprintf('   validity of the linear cold reference:\n');
fprintf('     (v_th/c)^2 = %.2e (warm corrections to the transverse response; v_th = %.3g m/s)\n', ...
        (plasma.vth/K.C0)^2, plasma.vth);
if plasma.slab
   spread = plasma.vth*P(1).t(end)/(plasma.b - plasma.a);
   fprintf('     v_th t_max / thickness = %.2e (spread of the slab edges over the run)%s\n', spread, ...
           warn_if(spread > 0.05));
end
w = K.C0/pulse.sigma;                                      % spectral peak of the derivative pulse (k sigma = 1)
if pulse.k0 > 0, w = K.C0*pulse.k0; end                    % the carrier
a0 = K.QE*K.C0*pulse.B0/(K.ME*w*K.C0);
fprintf('     a0 = e E/(m_e w c) = %.2e at w = %.3g rad/s = %.3g w_p (linearity)%s\n', a0, w, w/plasma.wp, ...
        warn_if(a0 > 0.1));
end

function m = warn_if(c)
m = ''; if c, m = '  <-- WARNING: outside the hypotheses'; end
end

% ===================================================================================================================
% self test of the Fresnel solution
% ===================================================================================================================
function out = selftest(K)
% Slab of 3 c/wp: (1) |R|^2 + |T|^2 = 1 above the cutoff, (2) R, T equal to an independent formulation, the
% characteristic matrix of a homogeneous layer (Born & Wolf, Principles of Optics, ch. 1, section not checked; the
% same content as the propagation/matching matrices of Orfanidis Secs. 5.1-5.2), time factor exp(-i w t),
% (3) E and B continuous at both faces. It checks the code, not the physics (both solve the same cold model).
wp = 1e10; nu = 1e-12*wp;
plasma = struct('ok', true, 'wp', wp, 'slab', true, 'a', 0.1, 'b', 0.1 + 3*K.C0/wp);
w = linspace(0.01, 4, 4001)'*wp;
w(abs(w - wp) < 1e-3*wp) = [];
[~, ~, r, t] = transfer(w, 0, plasma, nu, K.C0, 1);
out.energy = max(abs(abs(r(w > wp)).^2 + abs(t(w > wp)).^2 - 1));
n = sqrt(1 - wp^2./(w.*(w + 1i*nu)));
dl = w/K.C0.*n*(plasma.b - plasma.a);
m11 = cos(dl); m12 = -1i*sin(dl)./n; m21 = -1i*n.*sin(dl); m22 = cos(dl);
den = (m11 + m12) + (m21 + m22);
out.born_wolf = max([abs(r - ((m11 + m12) - (m21 + m22))./den); abs(t - 2./den)]);
h = 1e-9*(plasma.b - plasma.a);
out.continuity = 0;
for x = [plasma.a, plasma.b]
   for sgn = [1, -1]
      [E1, B1] = transfer(w, x - h, plasma, nu, K.C0, sgn);
      [E2, B2] = transfer(w, x + h, plasma, nu, K.C0, sgn);
      out.continuity = max([out.continuity; abs(E1 - E2); K.C0*abs(B1 - B2)]);
   end
end
out.below_cutoff_T = max(abs(t(w < 0.5*wp)));
fprintf('selftest: max ||R|^2+|T|^2-1| above cutoff = %.2e\n', out.energy);
fprintf('selftest: max |R,T - characteristic matrix| = %.2e\n', out.born_wolf);
fprintf('selftest: max jump of E, cB at the faces      = %.2e (step 1e-9 of the thickness)\n', out.continuity);
fprintf('selftest: max |T| below 0.5 wp (3 c/wp slab)  = %.2e\n', out.below_cutoff_T);
out.pass = out.energy < 1e-9 && out.born_wolf < 1e-9 && out.continuity < 1e-6;
fprintf('selftest: %s\n', ternary(out.pass, 'PASS', 'FAIL'));
end

function v = ternary(c, a, b)
if c, v = a; else, v = b; end
end
