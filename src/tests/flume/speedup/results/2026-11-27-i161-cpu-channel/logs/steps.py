"""Step cost per pair from i161 A/B outputs: (loop(110) - loop(10)) / 100; medians per label; also the step-end wait."""
import re, sys, statistics as st, collections
for path in sys.argv[1:]:
    runs = collections.defaultdict(dict); wait = collections.defaultdict(dict)
    for line in open(path, errors='replace'):
        m = re.search(r'\d\d:\d\d:\d\d (\S+) r(\d+)-n(\d+) rc 0 wall \S+ loop ([\d.]+), steps ([\d.]+), regrid ([\d.]+)', line)
        if m:
            t, r, n, loop, steps, rg = m.group(1), int(m.group(2)), int(m.group(3)), float(m.group(4)), float(m.group(5)), float(m.group(6))
            runs[t][(r, n)] = loop; wait[t][(r, n)] = rg
    print(f'== {path.split("/")[-1]}')
    meds = {}
    for t in runs:
        pc = [(runs[t][(r, 110)] - runs[t][(r, 10)]) / 100 for r in range(1, 10) if (r, 10) in runs[t] and (r, 110) in runs[t]]
        wt = [(wait[t][(r, 110)] - wait[t][(r, 10)]) / 100 for r in range(1, 10) if (r, 10) in wait[t] and (r, 110) in wait[t]]
        if pc:
            meds[t] = st.median(pc)
            print(f'  {t:8s} pairs ' + ' / '.join(f'{x:.3f}' for x in pc) + f'  median {meds[t]:.3f}  wait/step ' + ' / '.join(f'{x:.3f}' for x in wt))
    if meds:
        ref = next(iter(meds))
        print('  vs ' + ref + ': ' + ', '.join(f'{t} {100*(v/meds[ref]-1):+.1f} %' for t, v in meds.items()))
