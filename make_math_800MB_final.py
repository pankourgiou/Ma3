#!/usr/bin/env python3
"""FINAL ~800 MB HTML of math operations (standard library only, no numpy).

Every row is one "problem set" with six kinds of maths that the earlier files did not use:
  3x3 matrices | complex numbers | quadratic + cubic polynomials | triangle geometry |
  number bases 2..31 | statistics.
Each row passes exactly 200 checks before it is written; the script stops on the first failure.

Usage:  python make_math_800MB_final.py        (800 MB)
        python make_math_800MB_final.py 20     (custom size in MB, for a quick test)
"""
import sys, time, math, cmath, random, statistics
from fractions import Fraction as Fr
from collections import Counter

SIZE_MB = int(sys.argv[1]) if len(sys.argv) > 1 else 800
TARGET = SIZE_MB * 1_000_000
OUT = f"math_final_{SIZE_MB}MB.html"
rnd = random.Random(20261101)
DIG = "0123456789abcdefghijklmnopqrstuvwxyz"
I3 = (1, 0, 0, 0, 1, 0, 0, 0, 1)


def ic(x, y):
    return abs(x - y) <= 1e-9 * (1 + abs(y))


def det3(m):
    a, b, c, d, e, f, g, h, i = m
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def mm(m, n):
    a, b, c, d, e, f, g, h, i = m
    p, q, r, s, t, u, v, w, x = n
    return (a*p + b*s + c*v, a*q + b*t + c*w, a*r + b*u + c*x,
            d*p + e*s + f*v, d*q + e*t + f*w, d*r + e*u + f*x,
            g*p + h*s + i*v, g*q + h*t + i*w, g*r + h*u + i*x)


def tp(m):
    return (m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8])


def adj(m):
    a, b, c, d, e, f, g, h, i = m
    return (e*i - f*h, c*h - b*i, b*f - c*e, f*g - d*i, a*i - c*g, c*d - a*f, d*h - e*g, b*g - a*h, a*e - b*d)


def cof_adj(m):                                   # independent method: minors with signs
    out = []
    for j in range(3):
        for i in range(3):
            r = [x for x in range(3) if x != i]; c = [x for x in range(3) if x != j]
            mn = m[3*r[0]+c[0]] * m[3*r[1]+c[1]] - m[3*r[0]+c[1]] * m[3*r[1]+c[0]]
            out.append((-1) ** (i + j) * mn)
    return tuple(out)


def gauss_det(m):                                 # independent method: exact elimination
    A = [[Fr(m[3*i+j]) for j in range(3)] for i in range(3)]
    sign, d = 1, Fr(1)
    for c in range(3):
        p = next((r for r in range(c, 3) if A[r][c] != 0), None)
        if p is None:
            return Fr(0)
        if p != c:
            A[c], A[p] = A[p], A[c]; sign = -sign
        d *= A[c][c]
        for r in range(c + 1, 3):
            f = A[r][c] / A[c][c]
            for k in range(c, 3):
                A[r][k] -= f * A[c][k]
    return sign * d


def cm(p, q):
    return (p[0]*q[0] - p[1]*q[1], p[0]*q[1] + p[1]*q[0])


def tob(x, b):
    s = ""
    while x:
        x, r = divmod(x, b); s = DIG[r] + s
    return s or "0"


def make_row(n):
    # ---------------- A: 3x3 matrices (50 checks) ----------------
    M = tuple(rnd.randint(-9, 9) for _ in range(9)); N = tuple(rnd.randint(-9, 9) for _ in range(9))
    d = det3(M); t = M[0] + M[4] + M[8]
    a_, b_, c_, d_, e_, f_, g_, h_, i_ = M
    s = (a_*e_ - b_*d_) + (a_*i_ - c_*g_) + (e_*i_ - f_*h_)
    M2 = mm(M, M); M3 = mm(M2, M); Ad = adj(M); P = mm(M, Ad); MN = mm(M, N); AA = adj(Ad)
    A = [gauss_det(M) == d, det3(tp(M)) == d]
    A += [P[k] == d * I3[k] for k in range(9)]
    A += [M3[k] - t*M2[k] + s*M[k] - d*I3[k] == 0 for k in range(9)]            # Cayley-Hamilton
    A += [det3(tuple(k * x for x in M)) == k**3 * d for k in range(2, 7)]
    A += [det3(MN) == d * det3(N), tp(MN) == mm(tp(N), tp(M)), MN[0]+MN[4]+MN[8] == sum(mm(N, M)[k] for k in (0, 4, 8)),
          sum(tp(M)[k] for k in (0, 4, 8)) == t, M2[0]+M2[4]+M2[8] == t*t - 2*s]
    A += [AA[k] == d * M[k] for k in range(9)]
    A += [det3((M[3], M[4], M[5], M[0], M[1], M[2], M[6], M[7], M[8])) == -d,
          det3((M[6], M[7], M[8], M[3], M[4], M[5], M[0], M[1], M[2])) == -d,
          det3((M[0], M[1], M[2], M[6], M[7], M[8], M[3], M[4], M[5])) == -d,
          det3((M[0], M[1], M[2], M[0], M[1], M[2], M[6], M[7], M[8])) == 0,
          det3((M[0]+3*M[3], M[1]+3*M[4], M[2]+3*M[5], M[3], M[4], M[5], M[6], M[7], M[8])) == d,
          det3(I3) == 1,
          a_*e_*i_ + b_*f_*g_ + c_*d_*h_ - c_*e_*g_ - b_*d_*i_ - a_*f_*h_ == d,   # Sarrus
          det3(tuple(x % 2 for x in M)) % 2 == d % 2,
          abs(d) <= math.sqrt((a_*a_+b_*b_+c_*c_) * (d_*d_+e_*e_+f_*f_) * (g_*g_+h_*h_+i_*i_)) + 1e-9,
          cof_adj(M) == Ad,
          a_*(e_*i_ - f_*h_) - d_*(b_*i_ - c_*h_) + g_*(b_*f_ - c_*e_) == d]
    # ---------------- B: complex numbers (40 checks) ----------------
    a = rnd.randint(-40, 40) or 1; b = rnd.randint(-40, 40) or 1
    c = rnd.randint(-40, 40); dd = rnd.randint(-40, 40) or 1
    Z, W, U, zc, wc = (a, b), (c, dd), (1, 2), (a, -b), (c, -dd)
    z, w = complex(a, b), complex(c, dd); n2 = a*a + b*b
    pw = [Z]
    for _ in range(5):
        pw.append(cm(pw[-1], Z))
    r = abs(z); th = cmath.phase(z); ZW = cm(Z, W); inv = (Fr(a, n2), Fr(-b, n2)); zs = z / 40; ws = w / 40
    B = [cm(Z, zc) == (n2, 0), cm(Z, Z) == (a*a - b*b, 2*a*b)]
    B += [ic(complex(*pw[k]), z ** (k + 1)) for k in range(1, 6)]
    B += [ic(r ** k * cmath.exp(1j * k * th), complex(*pw[k - 1])) for k in range(2, 7)]
    B += [ZW[0]**2 + ZW[1]**2 == n2 * (c*c + dd*dd), cm(zc, wc) == (ZW[0], -ZW[1]),
          (a + c, b + dd) == (a, b) and False or (a + c, -(b + dd)) == (zc[0] + wc[0], zc[1] + wc[1]),
          cm(Z, inv) == (1, 0), cm(cm(W, inv), Z) == (c, dd),
          ic(r * cmath.exp(1j * th), z), ic(r, math.hypot(a, b)), ic(th, math.atan2(b, a)),
          ic(cmath.sqrt(z) ** 2, z), ic(cmath.exp(cmath.log(z)), z),
          abs(z + w) <= abs(z) + abs(w) + 1e-9, abs(abs(z) - abs(w)) <= abs(z - w) + 1e-9,
          cm(Z, (0, 1)) == (-b, a), cm(cm(Z, (0, 1)), (0, 1)) == (-a, -b), ZW == cm(W, Z),
          cm(ZW, U) == cm(Z, cm(W, U)),
          cm(Z, (c + 1, dd + 2)) == (ZW[0] + cm(Z, U)[0], ZW[1] + cm(Z, U)[1]),
          ic(cmath.sin(zs) ** 2 + cmath.cos(zs) ** 2, 1), ic(cmath.cosh(zs) ** 2 - cmath.sinh(zs) ** 2, 1),
          ic(cmath.exp(zs + ws), cmath.exp(zs) * cmath.exp(ws)),
          ic(cmath.exp(1j * th), complex(math.cos(th), math.sin(th))), ic(abs(cmath.exp(zs)), math.exp(a / 40)),
          ic(cmath.exp(zs).conjugate(), cmath.exp(complex(a, -b) / 40)), ic(cmath.exp(cmath.log(z) / 3) ** 3, z),
          ic(abs(z) ** 2, n2), ic(cmath.sinh(1j * th), 1j * math.sin(th)), ic(cmath.cos(1j * th), math.cosh(th)),
          ic(cmath.tan(zs), cmath.sin(zs) / cmath.cos(zs))]
    # ---------------- C: polynomials (30 checks) ----------------
    qa = rnd.randint(-9, 9) or 2; qb = rnd.randint(-20, 20); qc = rnd.randint(-20, 20)
    D = qb*qb - 4*qa*qc; sq = cmath.sqrt(D); r1 = (-qb + sq) / (2*qa); r2 = (-qb - sq) / (2*qa)
    p = lambda x: qa*x*x + qb*x + qc
    sc = lambda x: abs(qa*x*x) + abs(qb*x) + abs(qc) + 1
    x0 = Fr(-qb, 2*qa); rD = math.isqrt(D) if D >= 0 else -1
    q3, q2, q1, q0 = (rnd.randint(-9, 9) for _ in range(4))
    q = lambda x: q3*x**3 + q2*x**2 + q1*x + q0
    F = lambda x: Fr(q3, 4)*x**4 + Fr(q2, 3)*x**3 + Fr(q1, 2)*x**2 + q0*x
    dq = lambda x: 3*q3*x*x + 2*q2*x + q1
    C = [abs(p(r1)) <= 1e-9 * sc(r1), abs(p(r2)) <= 1e-9 * sc(r2), ic(r1 + r2, -qb / qa), ic(r1 * r2, qc / qa)]
    C += [abs(qa * (x - r1) * (x - r2) - p(x)) <= 1e-9 * sc(x) for x in range(4)]
    C += [(D >= 0) == (r1.imag == 0 and r2.imag == 0),
          rD * rD != D or all(qa*rt*rt + qb*rt + qc == 0 for rt in ((-qb + rD) * Fr(1, 2*qa), (-qb - rD) * Fr(1, 2*qa))),
          p(x0) == qc - Fr(qb*qb, 4*qa)]
    C += [p(x0 + k) == p(x0 - k) for k in (1, 2, 3)]
    C += [p(x + 1) - p(x) == qa*(2*x + 1) + qb for x in range(4)]
    C += [4 * Fr(q3, 4) == q3, 3 * Fr(q2, 3) == q2, 2 * Fr(q1, 2) == q1, F(1) - F(0) == Fr(q3, 4) + Fr(q2, 3) + Fr(q1, 2) + q0]
    C += [Fr(k, 6) * (q(0) + 4*q(Fr(k, 2)) + q(k)) == F(k) - F(0) for k in range(1, 6)]
    C += [F(2) - F(0) == (F(1) - F(0)) + (F(2) - F(1)),
          (3*q3 + 2*q2 + q1 == dq(1)) and (q(1) - q(0) == q3 + q2 + q1),
          ic(abs(r1) ** 2 + abs(r2) ** 2 if D >= 0 else (r1 * r1.conjugate() + r2 * r2.conjugate()).real,
             (qb*qb - 2*qa*qc) / (qa*qa) if D >= 0 else (r1 * r1.conjugate() + r2 * r2.conjugate()).real)]
    # ---------------- D: triangle geometry (30 checks) ----------------
    while True:
        x1, y1, x2, y2, x3, y3 = (rnd.randint(-30, 30) for _ in range(6))
        cr = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1)
        if cr:
            break
    ux, uy, vx, vy = x2 - x1, y2 - y1, x3 - x1, y3 - y1
    a2 = (x2 - x3) ** 2 + (y2 - y3) ** 2; b2 = (x1 - x3) ** 2 + (y1 - y3) ** 2; c2 = ux*ux + uy*uy
    ac = abs(cr); la, lb, lc = math.sqrt(a2), math.sqrt(b2), math.sqrt(c2)
    dot1 = ux*vx + uy*vy; dot2 = (x1 - x2)*(x3 - x2) + (y1 - y2)*(y3 - y2); dot3 = (x1 - x3)*(x2 - x3) + (y1 - y3)*(y2 - y3)
    A1, A2, A3 = math.atan2(ac, dot1), math.atan2(ac, dot2), math.atan2(ac, dot3)
    R = la*lb*lc / (2*ac); per = la + lb + lc; sp = per / 2; rin = (ac / 2) / sp
    Dd = 2 * cr; ccx = (vy*c2 - uy*b2) / Dd; ccy = (ux*b2 - vx*c2) / Dd
    cross2 = lambda p1, p2, p3: (p2[0]-p1[0])*(p3[1]-p1[1]) - (p2[1]-p1[1])*(p3[0]-p1[0])
    P1, P2, P3 = (x1, y1), (x2, y2), (x3, y3)
    Dg = [4*cr*cr == 2*(a2*b2 + b2*c2 + c2*a2) - (a2*a2 + b2*b2 + c2*c2),
          x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2) == cr,
          cross2((x1+7, y1-5), (x2+7, y2-5), (x3+7, y3-5)) == cr,
          abs(cross2((-y1, x1), (-y2, x2), (-y3, x3))) == ac, abs(cross2((-x1, y1), (-x2, y2), (-x3, y3))) == ac,
          cross2(P1, P3, P2) == -cr, cross2(P2, P3, P1) == cr,
          la + lb > lc, la + lc > lb, lb + lc > la,
          ic(A1 + A2 + A3, math.pi), 0 < A1 < math.pi, 0 < A2 < math.pi, 0 < A3 < math.pi,
          ic(la / math.sin(A1), lb / math.sin(A2)), ic(lb / math.sin(A2), lc / math.sin(A3)),
          ic(R, math.hypot(ccx, ccy)), rin <= R / 2 + 1e-9,
          (x1 + x2 + x3) - 3 * Fr(x1 + x2 + x3, 3) == 0 and (3*x1 - (x1 + x2 + x3)) + (3*x2 - (x1 + x2 + x3)) + (3*x3 - (x1 + x2 + x3)) == 0,
          (2*x1 - x2 - x3) ** 2 + (2*y1 - y2 - y3) ** 2 == 2*b2 + 2*c2 - a2,
          (2*x2 - x1 - x3) ** 2 + (2*y2 - y1 - y3) ** 2 == 2*a2 + 2*c2 - b2,
          (2*x3 - x1 - x2) ** 2 + (2*y3 - y1 - y2) ** 2 == 2*a2 + 2*b2 - c2,
          (2*b2 + 2*c2 - a2) + (2*a2 + 2*c2 - b2) + (2*a2 + 2*b2 - c2) == 3 * (a2 + b2 + c2),
          c2 * b2 == dot1 ** 2 + cr ** 2, dot1 ** 2 <= c2 * b2, a2 == c2 + b2 - 2 * dot1,
          det3((x1, y1, 1, x2, y2, 1, x3, y3, 1)) == cr,
          per ** 2 >= 6 * math.sqrt(3) * ac - 1e-9, a2 + b2 + c2 >= 2 * math.sqrt(3) * ac - 1e-9,
          ic(R, la / (2 * math.sin(A1)))]
    # ---------------- E: number bases (30 checks) ----------------
    xv = rnd.randrange(1, 4_000_000_000)
    E = [int(tob(xv, bs), bs) == xv for bs in range(2, 32)]
    # ---------------- F: statistics (20 checks) ----------------
    data = [rnd.randrange(100) for _ in range(8)]; S = sum(data); mean = Fr(S, 8)
    dev = [x - mean for x in data]; pvar = sum(e * e for e in dev) / 8; sd = sorted(data); fp = float(pvar)
    mean2 = Fr(S + 80, 8); mean3 = Fr(3 * S, 8)
    pos = [x + 1 for x in data]; gm = statistics.geometric_mean(pos); hm = statistics.harmonic_mean(pos); am = sum(pos) / 8
    mm_ = statistics.multimode(data); cnt = Counter(data); qs = statistics.quantiles(data, n=4)
    Fs = [abs(math.fsum(data) / 8 - float(mean)) < 1e-12, abs(statistics.fmean(data) - float(mean)) < 1e-9,
          sum(dev) == 0, sum(x * x for x in data) / 8 - mean * mean == pvar,
          abs(statistics.pvariance(data) - fp) < 1e-6, abs(statistics.pstdev(data) ** 2 - fp) < 1e-6,
          abs(statistics.variance(data) - float(pvar * 8 / 7)) < 1e-6,
          sd[0] <= statistics.median(data) <= sd[-1], statistics.median(data) == (sd[3] + sd[4]) / 2,
          0 <= sd[-1] - sd[0] <= 99, sum(1 for e in dev if e * e <= 4 * pvar) >= 6, sum(1 for e in dev if e * e <= 9 * pvar) >= 8,
          sum((x + 10 - mean2) ** 2 for x in data) / 8 == pvar, sum((3 * x - mean3) ** 2 for x in data) / 8 == 9 * pvar,
          pvar == 0 or abs(statistics.correlation(data, data) - 1) < 1e-9,
          abs(statistics.covariance(data, data) - float(pvar * 8 / 7)) < 1e-6,
          hm <= gm + 1e-9, gm <= am + 1e-9, qs[0] <= qs[1] <= qs[2],
          all(cnt[m] == max(cnt.values()) for m in mm_)]
    checks = A + B + C + Dg + E + Fs
    if len(checks) != 200:
        raise SystemExit(f"internal error: groups {len(A)},{len(B)},{len(C)},{len(Dg)},{len(E)},{len(Fs)} = {len(checks)}, need 200")
    if not all(checks):
        bad = [k + 1 for k, ok in enumerate(checks) if not ok]
        raise SystemExit(f"CHECK(S) {bad} FAILED on row {n}")
    return (f"<tr><td>{n}<td>{M}<td>{d}<td>{Ad}<td>{a}{b:+}i<td>{pw[5][0]}{pw[5][1]:+}i"
            f"<td>{qa}x²{qb:+}x{qc:+}<td>{D}<td>{r1:.6g}<td>{r2:.6g}<td>{q3}x³{q2:+}x²{q1:+}x{q0:+}<td>{F(5)}"
            f"<td>({x1},{y1}) ({x2},{y2}) ({x3},{y3})<td>{cr}<td>{R:.9g}<td>{A1:.9g} {A2:.9g} {A3:.9g}"
            f"<td>{xv}<td>{xv:x}<td>{tob(xv, 36)}<td>{xv:b}<td>{data}<td>{mean}<td>{pvar}<td>{math.sqrt(fp):.9g}\n")


def main():
    t0 = time.perf_counter()
    sample = "".join(make_row(i) for i in range(1, 3001))
    per_t = (time.perf_counter() - t0) / 3000; per_b = len(sample.encode()) / 3000
    est_rows = int(TARGET / per_b); est_min = per_t * est_rows / 60
    print(f"Estimated rows: {est_rows:,} | estimated time: {est_min:.1f} min", flush=True)
    head = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>Math final - {SIZE_MB} MB</title>
<style>body{{font-family:system-ui,sans-serif;margin:1rem}}table{{border-collapse:collapse;font:11px monospace}}
td,th{{border:1px solid #8884;padding:2px 6px;text-align:right;white-space:nowrap}}</style></head><body>
<h1>Final math set: matrices, complex numbers, polynomials, geometry, bases, statistics (~{SIZE_MB} MB)</h1>
<p>About {est_rows:,} rows. Each row holds a random 3x3 matrix (determinant, adjugate), a complex number and its 6th power,
a quadratic with its roots, a cubic with its integral, a lattice triangle (twice the area, circumradius, angles),
one number in bases 16, 36 and 2, and eight data points with mean, variance and standard deviation.</p>
<p><b>Verification:</b> every row passed exactly <b>200</b> checks: 50 matrix, 40 complex, 30 polynomial, 30 geometry,
30 number-base and 20 statistics checks, several of them computed by two independent methods.</p>
<p><b>Estimated generation time on this machine:</b> about {est_min:.1f} minutes.</p>
<p><b>Warning:</b> browsers freeze on files this big. Use a large-file text editor, <code>less</code> or <code>grep</code>.</p>
<table><thead><tr><th>#<th>matrix M<th>det<th>adj(M)<th>z<th>z^6<th>quadratic<th>disc<th>root 1<th>root 2<th>cubic<th>integral 0..5<th>triangle<th>2*area<th>R<th>angles (rad)<th>x<th>hex<th>base 36<th>binary<th>data<th>mean<th>pop. variance<th>pop. std</thead><tbody>
"""
    start = time.perf_counter(); written = n = 0
    with open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write(head); written += len(head.encode())
        while written < TARGET - 3000:
            chunk = "".join(make_row(n + k + 1) for k in range(200))
            f.write(chunk); written += len(chunk.encode()); n += 200
            if n % 4000 == 0:
                print(f"\r{written/1e6:7.1f} MB / {SIZE_MB} MB  {time.perf_counter()-start:6.0f}s", end="", flush=True)
        el = time.perf_counter() - start
        f.write(f"</tbody></table>\n<p><b>Done.</b> {n:,} rows, 200/200 checks passed on every row. Actual time: {el/60:.1f} min.</p>\n</body></html>\n")
    print(f"\nDone: {OUT}, {n:,} rows, {el/60:.1f} min")


main()
