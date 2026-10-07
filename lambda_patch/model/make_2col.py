# Make the 2-column twin of a 5-column synthetic set: keep "re im" of every line, drop "Q PiRe PiIm".
# Both writers format re/im with %.7g, so the result is byte-identical to what gen_dumps.py (or
# gen_synthetic_dumps.m with 'Cols', 2) writes for the same seeds: the CIR of the two sets is the same.
#   python3 make_2col.py <root5> <root2>        e.g. ../../synthetic_dumps ../../synthetic_dumps_2col
import os, sys, glob

SUBSETS = ('h0', 'h0_mp20', 'h1')                     # the 5-column-only checks are not copied


def main(root5, root2):
    n = 0
    for sub in SUBSETS:
        for fn in sorted(glob.glob(os.path.join(root5, sub, '*', 'bin', '*_RifCir_AccNum_*.txt'))):
            out = os.path.join(root2, os.path.relpath(fn, root5))
            os.makedirs(os.path.dirname(out), exist_ok=True)
            with open(fn) as f:
                rows = [ln.split()[:2] for ln in f if ln.strip()]
            assert len(rows) == 256, fn
            with open(out, 'w') as f:
                f.write(''.join('%s %s\n' % (a, b) for a, b in rows))
            n += 1
    print('%d dumps written under %s' % (n, root2))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
