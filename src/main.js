import * as THREE from "three";
import "./style.css";

/* =====================================================================
   PHYSICS (SI units). Verified against brute force Biot Savart sums.
   ===================================================================== */
const MU0 = 1.25663706212e-6, PI = Math.PI;

// Bulirsch cel, as given in Derby and Olbert, Am. J. Phys. 78, 229 (2010)
function cel(kc, p, c, s) {
  if (kc === 0) return NaN;
  const errtol = 1e-9;
  let k = Math.abs(kc), pp = p, cc = c, ss = s, em = 1, f, g, q;
  if (p > 0) { pp = Math.sqrt(p); ss = s / pp; }
  else { f = kc * kc; q = 1 - f; g = 1 - pp; f = f - pp; q = q * (ss - c * pp); pp = Math.sqrt(f / g); cc = (c - ss) / g; ss = -q / (g * g * pp) + cc * pp; }
  f = cc; cc = cc + ss / pp; g = k / pp; ss = 2 * (ss + f * g); pp = g + pp; g = em; em = k + em;
  let kk = k, it = 0;
  while (Math.abs(g - k) > g * errtol && it++ < 40) { k = 2 * Math.sqrt(kk); kk = k * em; f = cc; cc = cc + ss / pp; g = kk / pp; ss = 2 * (ss + f * g); pp = g + pp; g = em; em = k + em; }
  return (PI / 2) * (ss + cc * em) / (em * (em + pp));
}
function ellipKE(m) {
  let a = 1, b = Math.sqrt(1 - m), c2sum = m / 2, pow = 0.5;
  for (let i = 0; i < 30; i++) { const an = (a + b) / 2, bn = Math.sqrt(a * b), cn = (a - b) / 2; pow *= 2; c2sum += pow * cn * cn; a = an; b = bn; if (Math.abs(cn) < 1e-15 * a) break; }
  const K = PI / (2 * a); return [K, K * (1 - c2sum)];
}
// Uniformly magnetized cylinder: [B_rho, B_z] in its own frame
function magnetField(rho, z, a, b, Br) {
  const B0 = Br / PI; if (rho < 1e-12 * a) rho = 1e-12 * a;
  const zp = z + b, zm = z - b, ap = a + rho, am = a - rho;
  const sp = Math.sqrt(zp * zp + ap * ap), sm = Math.sqrt(zm * zm + ap * ap);
  const gam = am / ap;
  let kp = Math.sqrt((zp * zp + am * am) / (zp * zp + ap * ap)), km = Math.sqrt((zm * zm + am * am) / (zm * zm + ap * ap));
  if (kp < 1e-12) kp = 1e-12; if (km < 1e-12) km = 1e-12;
  const Brho = B0 * (a / sp * cel(kp, 1, 1, -1) - a / sm * cel(km, 1, 1, -1));
  const Bz = (B0 * a / ap) * (zp / sp * cel(kp, gam * gam, 1, gam) - zm / sm * cel(km, gam * gam, 1, gam));
  return [Brho, Bz];
}
// Current loop vector potential A_phi and field
function loopA(rho, z, R, I) {
  if (rho < 1e-9 * R) { const d = R * R + z * z; return MU0 * I * R * R * rho / (4 * d * Math.sqrt(d)); }
  const q = (R + rho) * (R + rho) + z * z, m = 4 * R * rho / q; const [K, E] = ellipKE(m);
  return (MU0 * I / (PI * Math.sqrt(m))) * Math.sqrt(R / rho) * ((1 - m / 2) * K - E);
}
function loopB(rho, z, R, I) {
  const q = (R + rho) * (R + rho) + z * z, m = 4 * R * rho / q; const [K, E] = ellipKE(m);
  const d = (R - rho) * (R - rho) + z * z, C = MU0 * I / (2 * PI), sq = Math.sqrt(q);
  const Bz = C / sq * (K + (R * R - rho * rho - z * z) / d * E);
  let Br = 0; if (rho > 1e-9 * R) Br = C * z / (rho * sq) * (-K + (R * R + rho * rho + z * z) / d * E);
  return [Br, Bz];
}
function gaussLegendre(n) {
  const x = [], w = [];
  for (let i = 1; i <= n; i++) {
    let r = Math.cos(PI * (i - 0.25) / (n + 0.5)), dp;
    for (let it = 0; it < 100; it++) { let p0 = 1, p1 = r; for (let k = 2; k <= n; k++) { const p2 = ((2 * k - 1) * r * p1 - (k - 1) * p0) / k; p0 = p1; p1 = p2; } dp = n * (r * p1 - p0) / (r * r - 1); const dr = p1 / dp; r -= dr; if (Math.abs(dr) < 1e-15) break; }
    x.push(r); w.push(2 / ((1 - r * r) * dp * dp));
  }
  return { x, w };
}
function perpBasis(n) {
  const t = Math.abs(n[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
  let u = [n[1] * t[2] - n[2] * t[1], n[2] * t[0] - n[0] * t[2], n[0] * t[1] - n[1] * t[0]];
  const l = Math.hypot(u[0], u[1], u[2]); u = [u[0] / l, u[1] / l, u[2] / l];
  return [u, [n[1] * u[2] - n[2] * u[1], n[2] * u[0] - n[0] * u[2], n[0] * u[1] - n[1] * u[0]]];
}
function makeQuad(nz, nphi) { const gl = gaussLegendre(nz), cs = [], sn = []; for (let j = 0; j < nphi; j++) { cs.push(Math.cos(2 * PI * j / nphi)); sn.push(Math.sin(2 * PI * j / nphi)); } return { gl, cs, sn, nphi }; }
// Flux linkage of the coil from the magnet, by reciprocity (see notes at the end of the file)
function fluxLinkage(mag, coil, quad) {
  const M = mag.Br / MU0, [u, v] = perpBasis(mag.axis), n = coil.axis, { gl, cs, sn, nphi } = quad;
  let total = 0;
  for (let iz = 0; iz < gl.x.length; iz++) {
    const zs = gl.x[iz] * mag.b, wz = gl.w[iz] * mag.b;
    const cx = mag.pos[0] + mag.axis[0] * zs, cy = mag.pos[1] + mag.axis[1] * zs, cz = mag.pos[2] + mag.axis[2] * zs;
    let circ = 0;
    for (let j = 0; j < nphi; j++) {
      const ex = u[0] * cs[j] + v[0] * sn[j], ey = u[1] * cs[j] + v[1] * sn[j], ez = u[2] * cs[j] + v[2] * sn[j];
      const px = cx + mag.a * ex, py = cy + mag.a * ey, pz = cz + mag.a * ez;
      const tx = -u[0] * sn[j] + v[0] * cs[j], ty = -u[1] * sn[j] + v[1] * cs[j], tz = -u[2] * sn[j] + v[2] * cs[j];
      for (const L of coil.loops) {
        const dx = px - coil.pos[0] - n[0] * L.z, dy = py - coil.pos[1] - n[1] * L.z, dz = pz - coil.pos[2] - n[2] * L.z;
        const zl = dx * n[0] + dy * n[1] + dz * n[2];
        const rx = dx - zl * n[0], ry = dy - zl * n[1], rz = dz - zl * n[2];
        const rho = Math.hypot(rx, ry, rz); if (rho < 1e-15) continue;
        const A = loopA(rho, zl, L.r, 1);
        circ += L.w * A * ((n[1] * rz - n[2] * ry) * tx + (n[2] * rx - n[0] * rz) * ty + (n[0] * ry - n[1] * rx) * tz) / rho;
      }
    }
    total += wz * M * circ * (2 * PI * mag.a / nphi);
  }
  return total * coil.N;
}
function coilLoops(rIn, rOut, halfLen, nz, nr) {
  const gz = gaussLegendre(nz), gr = gaussLegendre(nr), loops = [];
  for (let i = 0; i < nz; i++) for (let j = 0; j < nr; j++) loops.push({ z: gz.x[i] * halfLen, r: (rIn + rOut) / 2 + gr.x[j] * (rOut - rIn) / 2, w: gz.w[i] * gr.w[j] / 4 });
  return loops;
}
// Meridional field line tracer, mirrored across z = 0
function traceMeridional(fieldFn, rho0, rhoCut, hFn, rMax) {
  const pts = [[rho0, 0]]; let r = rho0, z = 0, closed = false;
  const dir = (rr, zz) => { const [a, b] = fieldFn(Math.abs(rr), zz); const s = rr < 0 ? -1 : 1; const l = Math.hypot(a, b) || 1; return [s * a / l, b / l]; };
  for (let i = 0; i < 6000; i++) {
    const h = hFn(r, z);
    const k1 = dir(r, z), k2 = dir(r + h / 2 * k1[0], z + h / 2 * k1[1]), k3 = dir(r + h / 2 * k2[0], z + h / 2 * k2[1]), k4 = dir(r + h * k3[0], z + h * k3[1]);
    const nr = r + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]), nz = z + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
    if (i > 2 && z > 0 && nz <= 0 && nr > rhoCut) { const t = z / (z - nz); pts.push([r + t * (nr - r), 0]); closed = true; break; }
    r = nr; z = nz; pts.push([Math.abs(r), z]); if (Math.hypot(r, z) > rMax) break;
  }
  if (closed) { const full = pts.slice(); for (let i = pts.length - 2; i >= 1; i--) full.push([pts[i][0], -pts[i][1]]); return { pts: full, closed: true }; }
  const back = []; for (let i = pts.length - 1; i >= 1; i--) back.push([pts[i][0], -pts[i][1]]);
  return { pts: back.concat(pts), closed: false };
}
function equalFluxRadii(BzFn, rhoEdge, nRings) {
  const n = 400, cum = [0], rs = [0];
  for (let i = 1; i <= n; i++) { const r0 = rhoEdge * (i - 1) / n, r1 = rhoEdge * i / n, rm = (r0 + r1) / 2; cum.push(cum[i - 1] + 2 * PI * rm * BzFn(rm) * (r1 - r0)); rs.push(r1); }
  const total = cum[n], radii = [];
  for (let k = 0; k < nRings; k++) { const target = total * (k + 0.5) / nRings; let i = 1; while (cum[i] < target) i++; const t = (target - cum[i - 1]) / (cum[i] - cum[i - 1]); radii.push(rs[i - 1] + t * (rs[i] - rs[i - 1])); }
  return { radii, total };
}

/* =====================================================================
   APPARATUS
   ===================================================================== */
// Magnet materials: remanence Br and density. Values are typical datasheet figures.
const MAGNETS = {
  n52: { sym: 'N52', nm: 'NdFeB', name: 'Neodymium iron boron, N52', Br: 1.44, rho: 7500 },
  n42: { sym: 'N42', nm: 'NdFeB', name: 'Neodymium iron boron, N42', Br: 1.30, rho: 7500 },
  smco: { sym: 'SmCo', nm: 'Sm₂Co₁₇', name: 'Samarium cobalt', Br: 1.10, rho: 8400 },
  fer: { sym: 'C8', nm: 'Ferrite', name: 'Strontium ferrite, C8', Br: 0.39, rho: 4900 }
};
// Wire materials: resistivity at 20 C
const WIRES = {
  cu: { sym: 'Cu', nm: 'Copper', rho: 1.68e-8, color: 0xc27a45, metal: 1 },
  ag: { sym: 'Ag', nm: 'Silver', rho: 1.59e-8, color: 0xd9dde2, metal: 1 },
  al: { sym: 'Al', nm: 'Aluminium', rho: 2.65e-8, color: 0xb3bac3, metal: 1 },
  nicr: { sym: 'NiCr', nm: 'Nichrome', rho: 1.10e-6, color: 0x8b8c8f, metal: 0.8 }
};
const MAGNET = { a: 0.006, b: 0.0125, Br: 1.30, rho: 7500 };   // cylinder 12 mm across, 25 mm long
const BOBBIN = { bore: 0.015, rIn: 0.016, halfLen: 0.006, flange: 0.001 };
const WIRE_D = 0.5e-3;          // copper wire diameter, m
const R_METER = 1.0;            // current meter in series, ohm
const TRACK = 0.15, SWING = 0.12, S = 100; // metres; render scale 1 unit = 1 cm
const magMoment = () => Math.max(MAGNET.Br, 0.05) / MU0 * PI * MAGNET.a * MAGNET.a * 2 * MAGNET.b;

const state = {
  s: -0.12, v: 0, mode: 'swing', swing: { phase: -PI / 2, c: 0, A: SWING, w: 0, speed: 0.3 },
  sTarget: -0.12, sTargetPrev: -0.12, vTarget: 0, dragging: false,
  psi: 0, psiDot: 0, flip: null,
  theta: 0, thetaDot: 0, thetaTarget: 0, coilWaiting: false,
  mat: 'n42', wire: 'cu', show: { lines: true, map: true }, exp: 'coil', core: false, corePending: false, J: 0, F0: 0,
  N: 200, V: 0, I: 0, emf: 0, lam: 0, Gs: 0, force: 0, forceBatt: 0, t: 0, samples: []
};
let coil = null;
const quad = makeQuad(8, 24);

function buildCoil() {
  const pitch = WIRE_D * 1.07, perLayer = Math.floor(2 * BOBBIN.halfLen / pitch);
  const layers = Math.ceil(state.N / perLayer), depth = layers * pitch;
  const rOut = BOBBIN.rIn + depth, rMean = BOBBIN.rIn + depth / 2;
  const R = WIRES[state.wire].rho * state.N * 2 * PI * rMean / (PI * WIRE_D * WIRE_D / 4);
  const inch = 0.0254, rI = rMean / inch, lI = 2 * BOBBIN.halfLen / inch, dI = depth / inch;
  const L = 0.8 * rI * rI * state.N * state.N / (6 * rI + 9 * lI + 10 * dI) * 1e-6; // Wheeler multilayer
  coil = { perLayer, layers, depth, rOut, rMean, R, L,
    loops: coilLoops(BOBBIN.rIn, rOut, BOBBIN.halfLen, 4, depth > 0.004 ? 3 : 2),
    traceLoops: coilLoops(BOBBIN.rIn, rOut, BOBBIN.halfLen, 2, 2) };
  coreCoil();
}
/* ---------------------------------------------------------------------
   IRON CORE. A soft iron rod in the bore, 27 mm across and 40 mm long.
   Approximation: the rod is magnetized uniformly along the coil axis, with
     J = mu0 M = chi_eff * (average applied axial field), capped at Bs.
   chi_eff = chi / (1 + N chi), with N the rod's demagnetizing factor,
   computed from the exact field of a uniformly magnetized cylinder.
   Forces and EMF come from one co-energy, so they stay consistent:
     W' = I Lambda_m(s) + (V / 2 mu0) chi_eff (Bm(s) + beta I)^2
   Laminated iron is assumed, so the core carries no eddy currents.
   --------------------------------------------------------------------- */
const IRON = { rc: 0.0135, h: 0.02, chi: 5000, Bs: 2.0 };
const core = (() => {
  const ga = gaussLegendre(6), gr = gaussLegendre(4), pts = []; let ws = 0;
  for (let i = 0; i < 6; i++) for (let j = 0; j < 4; j++) { const r = (gr.x[j] + 1) / 2 * IRON.rc, w = ga.w[i] * gr.w[j] * r; pts.push([ga.x[i] * IRON.h, r, w]); ws += w; }
  pts.forEach(p => p[2] /= ws);
  let own = 0; for (const [z, r, w] of pts) own += w * magnetField(r, z, IRON.rc, IRON.h, 1)[1];
  const Nd = 1 - own, chiE = IRON.chi / (1 + Nd * IRON.chi);
  const cs = [], sn = []; for (let k = 0; k < 8; k++) { cs.push(Math.cos(PI * k / 4)); sn.push(Math.sin(PI * k / 4)); }
  return { pts, Nd, chiE, V: PI * IRON.rc * IRON.rc * 2 * IRON.h, cs, sn, beta: 0, K: 0, Ladd: 0 };
})();
function coreCoil() {   // average coil field per amp over the rod, and the rod's flux through the coil per tesla
  let beta = 0; for (const [z, r, w] of core.pts) { let b = 0; for (const L of coil.loops) b += loopB(r, z - L.z, L.r, state.N * L.w)[1]; beta += w * b; }
  core.beta = beta; core.K = core.V * beta / MU0; core.Ladd = core.K * core.chiE * beta;
}
function magB(px, py, pz, ms, m, out) {
  const dx = px - ms, z = dx * m[0] + py * m[1] + pz * m[2], rx = dx - z * m[0], ry = py - z * m[1], rz = pz - z * m[2], rho = Math.hypot(rx, ry, rz);
  const f = magnetField(rho, z, MAGNET.a, MAGNET.b, MAGNET.Br), ir = rho > 1e-12 ? f[0] / rho : 0;
  out[0] = ir * rx + f[1] * m[0]; out[1] = ir * ry + f[1] * m[1]; out[2] = ir * rz + f[1] * m[2];
}
function coreBm(s, psi, th) {   // magnet field along the coil axis, averaged over the rod
  const n = coilAxis(th), [u, v] = perpBasis(n), m = magAxis(psi), B = [0, 0, 0]; let sum = 0;
  for (const [z, r, w] of core.pts) for (let k = 0; k < 8; k++) {
    const c = core.cs[k] * r, q = core.sn[k] * r;
    magB(n[0] * z + u[0] * c + v[0] * q, n[1] * z + u[1] * c + v[1] * q, n[2] * z + u[2] * c + v[2] * q, s, m, B);
    sum += w / 8 * (B[0] * n[0] + B[1] * n[1] + B[2] * n[2]);
  }
  return sum;
}
const coreJ = (Bm, I) => { const j = core.chiE * (Bm + core.beta * I); return Math.max(-IRON.Bs, Math.min(IRON.Bs, j)); };

const magAxis = (psi) => [Math.cos(psi), 0, -Math.sin(psi)];
const coilAxis = (th) => [Math.cos(th), 0, -Math.sin(th)];
const magnetMass = () => MAGNET.rho * PI * MAGNET.a * MAGNET.a * 2 * MAGNET.b;
function lambda(s, psi, th) {
  return fluxLinkage({ pos: [s, 0, 0], axis: magAxis(psi), a: MAGNET.a, b: MAGNET.b, Br: MAGNET.Br },
    { pos: [0, 0, 0], axis: coilAxis(th), N: state.N, loops: coil.loops }, quad);
}
function innerTurnFlux(s, psi, th) {
  return fluxLinkage({ pos: [s, 0, 0], axis: magAxis(psi), a: MAGNET.a, b: MAGNET.b, Br: MAGNET.Br },
    { pos: [0, 0, 0], axis: coilAxis(th), N: 1, loops: [{ r: BOBBIN.rIn, z: 0, w: 1 }] }, quad);
}

// Collision between the magnet and the coil body, by sampling both surfaces
const colSamples = (() => {
  const pts = [];
  for (let i = 0; i <= 10; i++) for (let j = 0; j < 24; j++) { const t = 2 * PI * j / 24; pts.push([-1 + 2 * i / 10, Math.cos(t), Math.sin(t), 1]); }
  for (const e of [-1, 1]) for (const f of [0.35, 0.7]) for (let j = 0; j < 12; j++) { const t = 2 * PI * j / 12; pts.push([e, Math.cos(t), Math.sin(t), f]); }
  return pts;
})();
function collides(s, psi, th, coreOn = state.core) {
  const m = magAxis(psi), n = coilAxis(th), [u, v] = perpBasis(m);
  const zLim = BOBBIN.halfLen + BOBBIN.flange + 0.0005, rLo = BOBBIN.bore - 0.0005, rHi = coil.rOut + 0.002;
  for (const p of colSamples) {
    const zz = p[0] * MAGNET.b, rr = MAGNET.a * p[3];
    const x = s + m[0] * zz + rr * (u[0] * p[1] + v[0] * p[2]), y = m[1] * zz + rr * (u[1] * p[1] + v[1] * p[2]), z = m[2] * zz + rr * (u[2] * p[1] + v[2] * p[2]);
    const zc = x * n[0] + y * n[1] + z * n[2], rc = Math.hypot(x - zc * n[0], y - zc * n[1], z - zc * n[2]);
    if (coreOn && Math.abs(zc) < IRON.h + 0.0005 && rc < IRON.rc + 0.0005) return true;
    if (Math.abs(zc) > zLim) continue;
    if (rc > rLo && rc < rHi) return true;
  }
  const [cu, cv] = perpBasis(n), rings = [[-zLim, BOBBIN.bore], [0, BOBBIN.bore], [zLim, BOBBIN.bore]];
  if (coreOn) rings.push([-IRON.h, IRON.rc], [IRON.h, IRON.rc], [-IRON.h, IRON.rc * 0.5], [IRON.h, IRON.rc * 0.5], [-IRON.h, 0], [IRON.h, 0]);
  for (const [zc, R] of rings) for (let j = 0; j < 48; j++) {
    const t = 2 * PI * j / 48, c = Math.cos(t), sn = Math.sin(t);
    const x = n[0] * zc + R * (cu[0] * c + cv[0] * sn), y = n[1] * zc + R * (cu[1] * c + cv[1] * sn), z = n[2] * zc + R * (cu[2] * c + cv[2] * sn);
    const dx = x - s, zm = dx * m[0] + y * m[1] + z * m[2]; if (Math.abs(zm) > MAGNET.b + 0.0005) continue;
    if (Math.hypot(dx - zm * m[0], y - zm * m[1], z - zm * m[2]) < MAGNET.a + 0.0005) return true;
    if (R === 0) break;
  }
  return false;
}
// How far the magnet can swing from its home side before touching the coil
const freeCache = new Map();
function freeEnd(side, th, psi) {
  const coreOn = state.core || state.corePending, key = side + '|' + (th * 180 / PI).toFixed(1) + '|' + psi.toFixed(3) + '|' + state.N + '|' + coreOn;
  if (freeCache.has(key)) return freeCache.get(key);
  let end = side * SWING;
  for (let k = 0; k <= 120; k++) { const x = side * (SWING - k * 0.002); if (collides(x, psi, th, coreOn)) break; end = x; }
  const r = Math.abs(end + side * SWING) < 1e-9 ? -side * SWING : end - side * 0.004;
  freeCache.set(key, r); return r;
}
function swingRange() {
  const side = state.s <= 0 ? -1 : 1;
  const e1 = freeEnd(side, state.theta, state.psi), e2 = freeEnd(side, state.thetaTarget, state.psi);
  const inner = side < 0 ? Math.min(e1, e2) : Math.max(e1, e2);
  return side < 0 ? [-SWING, Math.max(inner, -SWING + 0.02)] : [Math.min(inner, SWING - 0.02), SWING];
}

/* =====================================================================
   FLUX TABLE along the track: exact values, cubic interpolation for
   dLambda/ds at every substep. Rebuilt when the coil or magnet turns.
   ===================================================================== */
const TAB = { key: null, vals: new Float64Array(301), n: 0, x0: -0.15, h: 0.001, N: 301 };
const tabKey = () => state.mat + '|' + state.N + '|' + state.psi.toFixed(4) + '|' + state.theta.toFixed(4);
const TABM = { vals: new Float64Array(301), n: 0, x0: -0.15, h: 0.001, N: 301 };
function fillTable(budgetMs) {
  const t0 = performance.now(), k = tabKey();
  if (TAB.key !== k) { TAB.key = k; TAB.n = 0; TABM.n = 0; TAB.psi = state.psi; TAB.th = state.theta; }
  if (state.flip || state.exp !== 'coil') return;
  while (TAB.n < TAB.N && performance.now() - t0 < budgetMs) { TAB.vals[TAB.n] = lambda(TAB.x0 + TAB.n * TAB.h, TAB.psi, TAB.th); TAB.n++; }
  if (state.core) while (TABM.n < TABM.N && performance.now() - t0 < budgetMs) { TABM.vals[TABM.n] = coreBm(TABM.x0 + TABM.n * TABM.h, TAB.psi, TAB.th); TABM.n++; }
}
function tabEval(x, T = TAB) {
  const u = (x - T.x0) / T.h, f = T.vals; if (u < 1 || u > T.N - 3) return null;
  const i = Math.floor(u), t = u - i, t2 = t * t, t3 = t2 * t;
  const f0 = f[i], f1 = f[i + 1], m0 = (f[i + 1] - f[i - 1]) / 2, m1 = (f[i + 2] - f[i]) / 2;
  return [(2 * t3 - 3 * t2 + 1) * f0 + (t3 - 2 * t2 + t) * m0 + (-2 * t3 + 3 * t2) * f1 + (t3 - t2) * m1,
    ((6 * t2 - 6 * t) * f0 + (3 * t2 - 4 * t + 1) * m0 + (-6 * t2 + 6 * t) * f1 + (3 * t2 - 2 * t) * m1) / T.h];
}

/* =====================================================================
   SIMULATION STEP
   Circuit: coil, meter and battery in series.
     L dI/dt = V_battery + EMF - R I,   EMF = -dLambda/dt from the magnet
   Magnet: driven when swinging or parked. When you drag it, it rides a
   spring to your cursor and feels the coil force G I, with G = dLambda/ds.
   ===================================================================== */
function simulate(dt) {
  const sT0 = state.sTargetPrev, sT1 = state.sTarget; state.sTargetPrev = sT1;
  if (state.dragging) state.vTarget += ((sT1 - sT0) / dt - state.vTarget) * Math.min(1, dt / 0.05); else state.vTarget = 0;
  // the coil angle follows the slider, but waits while the magnet is in the way
  const w = 8, e = state.thetaTarget - state.theta;
  let thDotC = state.thetaDot + (w * w * e - 2 * w * state.thetaDot) * dt, thNew = state.theta + thDotC * dt;
  state.coilWaiting = false;
  if (Math.abs(thNew - state.theta) > 1e-9 && collides(state.s, state.psi, thNew)) { thNew = state.theta; thDotC = 0; state.coilWaiting = Math.abs(e) > 0.01; }
  state.thetaDot = thDotC;
  const th0 = state.theta, thDot = (thNew - th0) / dt; state.theta = thNew;
  // magnet flip, driven over 0.8 s
  if (state.flip) {
    const f = state.flip; f.t += dt; const x = Math.min(1, f.t / f.dur);
    state.psi = f.from + (f.to - f.from) * x * x * (3 - 2 * x); state.psiDot = (f.to - f.from) * 6 * x * (1 - x) / f.dur;
    if (x >= 1) { state.psi = ((f.to + PI) % (2 * PI) + 2 * PI) % (2 * PI) - PI; state.psiDot = 0; state.flip = null; }
  }
  // swing range shrinks to stay clear of a turned coil
  const sw = state.swing;
  if (state.mode === 'swing') {
    const [lo, hi] = swingRange(), k = Math.min(1, dt / 0.6);
    sw.c += ((lo + hi) / 2 - sw.c) * k; sw.A += ((hi - lo) / 2 - sw.A) * k;
    sw.w = Math.min(2 * PI * 1.5, sw.speed / Math.max(sw.A, 0.005));
  }
  const useT = TAB.n === TAB.N && TAB.key === tabKey();
  const hS = 2e-4, psiM = state.psi - state.psiDot * dt / 2, thM = thNew - thDot * dt / 2, sM = state.s + state.v * dt / 2;
  const GsF = useT ? 0 : (lambda(sM + hS, psiM, thM) - lambda(sM - hS, psiM, thM)) / (2 * hS);
  const GthF = Math.abs(thDot) > 1e-6 ? (lambda(sM, psiM, thM + 1e-3) - lambda(sM, psiM, thM - 1e-3)) / 2e-3 : 0;
  const GpsiF = Math.abs(state.psiDot) > 1e-6 ? (lambda(sM, psiM + 1e-3, thM) - lambda(sM, psiM - 1e-3, thM)) / 2e-3 : 0;
  const Eo = -(GthF * thDot + GpsiF * state.psiDot);
  // iron core: applied field from the magnet, averaged over the rod, and its slope
  if (state.corePending !== state.core && (!state.corePending || !collides(state.s, state.psi, state.theta, true))) setCoreNow(state.corePending);
  const useM = state.core && useT && TABM.n === TABM.N;
  let BmF = 0, BmdF = 0;
  if (state.core && !useM) { BmF = coreBm(sM, psiM, thM); BmdF = (coreBm(sM + 1e-3, psiM, thM) - coreBm(sM - 1e-3, psiM, thM)) / 2e-3; }
  const coreTerms = (xm, Gs, I) => {   // effective slope, inductance and the iron's own pull on the magnet
    if (!state.core) return [Gs, coil.L, 0, 0];
    let Bm = BmF, Bmd = BmdF; if (useM) { const r = tabEval(xm, TABM); if (r) { Bm = r[0]; Bmd = r[1]; } }
    const jl = core.chiE * (Bm + core.beta * I);
    if (Math.abs(jl) < IRON.Bs) return [Gs + core.K * core.chiE * Bmd, coil.L + core.Ladd, core.V / MU0 * core.chiE * Bm * Bmd, Bm];
    return [Gs, coil.L, core.V / MU0 * Math.sign(jl) * IRON.Bs * Bmd, Bm];   // saturated rod
  };
  let Bmlast = 0, F0 = 0;

  const m = magnetMass(), Rt = coil.R + R_METER, L = coil.L, n = 8, h2 = dt / n, w0 = 30, V = state.V;
  let I = state.I, v = state.v, x = state.s, emf = 0, Gs = GsF, Geff = GsF;
  const samples = [];
  for (let k = 0; k < n; k++) {
    if (state.mode === 'drag') {
      if (useT) { const r = tabEval(x + v * h2 / 2); if (r) Gs = r[1]; }
      const [Ge, Le, f0, bm] = coreTerms(x + v * h2 / 2, Gs, I); F0 = f0; Bmlast = bm; Geff = Ge;
      const sT = sT0 + (sT1 - sT0) * k / n;
      const Fsp = m * w0 * w0 * (sT - x) + 2 * m * w0 * (state.vTarget - v) + f0;
      const D = Le / h2 + Rt, alpha = (Le * I / h2 + V + Eo) / D, beta = Ge / D;
      v = (m * v + h2 * Fsp + h2 * Ge * alpha) / (m + h2 * Ge * beta); I = alpha - beta * v; x += v * h2;
    } else {
      let xNew = x;
      if (state.mode === 'swing') {
        const ph = sw.phase + sw.w * h2; xNew = sw.c + sw.A * Math.sin(ph);
        if (collides(xNew, state.psi, state.theta)) { xNew = x; sw.phase = PI - sw.phase; } else sw.phase = ph;
      }
      if (useT) { const r = tabEval((x + xNew) / 2); if (r) Gs = r[1]; }
      const [Ge, Le, f0, bm] = coreTerms((x + xNew) / 2, Gs, I); F0 = f0; Bmlast = bm; Geff = Ge;
      v = (xNew - x) / h2; x = xNew;
      const D = Le / h2 + Rt; I = (Le * I / h2 + V + Eo - Ge * v) / D;
    }
    emf = -Geff * v + Eo;
    samples.push([state.t + (k + 1) * h2, emf, I - V / Rt]);
  }
  if (Math.abs(x) > TRACK) { x = Math.sign(x) * TRACK; v = 0; }
  if (state.mode === 'drag' && collides(x, state.psi, state.theta)) { x = state.s; v = 0; }
  state.s = x; state.v = state.mode === 'swing' || state.mode === 'drag' ? v : 0;
  state.I = I; state.emf = emf; state.Gs = Geff; state.force = Geff * I + F0; state.forceBatt = Geff * V / Rt; state.F0 = F0;
  state.J = state.core ? coreJ(Bmlast, I) : 0;
  const lt = useT ? tabEval(x) : null; state.lam = (lt ? lt[0] : lambda(x, state.psi, state.theta)) + (state.core ? core.K * state.J : 0);
  state.samples = samples; state.t += dt;
}
function flipMagnet() {
  if (state.flip) return;
  for (let i = 1; i <= 12; i++) if (collides(state.s, state.psi + PI * i / 12, state.theta)) { ui.toast('The magnet cannot turn over inside the coil. It will fit once it swings out.'); return; }
  state.flip = { from: state.psi, to: state.psi + PI, t: 0, dur: 0.8 };
}
function setSwing(on) {
  if (on) {
    const [lo, hi] = swingRange(), sw = state.swing; sw.c = (lo + hi) / 2; sw.A = (hi - lo) / 2;
    sw.phase = Math.asin(Math.max(-1, Math.min(1, (state.s - sw.c) / sw.A)));
    state.mode = 'swing';
  } else if (state.mode === 'swing') state.mode = 'park';
  $('swing').setAttribute('aria-pressed', on); $('swing').querySelector('span').textContent = on ? 'Swinging' : 'Swing';
}

/* =====================================================================
   FIELD LINES
   With the battery off, the coil's own field is under 3 percent of the
   magnet's, so the lines are the magnet's exact field, traced once and
   carried with it. With a stronger coil current, the combined field of
   magnet plus coil is traced live in 3D, a few lines per frame.
   ===================================================================== */
const NL = 32, GOLD = PI * (3 - Math.sqrt(5)), LG0 = Math.log10(3e-4), LG1 = Math.log10(1.3);
// Viridis colormap on a log scale of |B|, from 0.3 mT to 1.3 T. The color bar uses the same stops.
const CMAP = [[0, 0x44, 0x01, 0x54], [0.25, 0x3b, 0x52, 0x8b], [0.5, 0x21, 0x91, 0x8c], [0.75, 0x5e, 0xc9, 0x62], [1, 0xfd, 0xe7, 0x25]];
const bT = (B) => Math.min(1, Math.max(0, (Math.log10(B + 1e-12) - LG0) / (LG1 - LG0)));
function ramp(t, out, o) {
  let i = 1; while (i < CMAP.length - 1 && t > CMAP[i][0]) i++;
  const a = CMAP[i - 1], b = CMAP[i], f = Math.min(1, Math.max(0, (t - a[0]) / (b[0] - a[0])));
  out[o] = (a[1] + (b[1] - a[1]) * f) / 255; out[o + 1] = (a[2] + (b[2] - a[2]) * f) / 255; out[o + 2] = (a[3] + (b[3] - a[3]) * f) / 255;
}
const lineColor = (B, out, o) => { ramp(0.1 + 0.9 * bT(B), out, o); }; // lift the floor so weak lines stay visible
const seedRadii = equalFluxRadii(r => magnetField(r, 0, MAGNET.a, MAGNET.b, MAGNET.Br)[1], MAGNET.a * 0.9999, NL).radii;

class TubeSet {
  constructor(maxPts, radius) {
    this.r = radius; this.maxV = maxPts * 4;
    const g = new THREE.BufferGeometry();
    this.pos = new Float32Array(this.maxV * 3); this.col = new Float32Array(this.maxV * 3); this.idx = new Uint32Array(maxPts * 24);
    g.setAttribute('position', new THREE.BufferAttribute(this.pos, 3).setUsage(THREE.DynamicDrawUsage));
    g.setAttribute('color', new THREE.BufferAttribute(this.col, 3).setUsage(THREE.DynamicDrawUsage));
    g.setIndex(new THREE.BufferAttribute(this.idx, 1).setUsage(THREE.DynamicDrawUsage));
    this.geo = g; this.mesh = new THREE.Mesh(g, new THREE.MeshBasicMaterial({ vertexColors: true })); this.mesh.frustumCulled = false;
    this.arrows = new THREE.InstancedMesh(new THREE.ConeGeometry(0.075, 0.27, 8), new THREE.MeshBasicMaterial({ color: 0xf2efe8 }), 900);
    this.arrows.frustumCulled = false; this.arrows.count = 0;
  }
  // lines: [{ p: [x,y,z,...] scene units, b: [|B|...], closed }]
  set(lines) {
    let nv = 0, ni = 0, na = 0; const P = this.pos, C = this.col, X = this.idx, r = this.r;
    const M = new THREE.Matrix4(), q = new THREE.Quaternion(), up = new THREE.Vector3(0, 1, 0), tv = new THREE.Vector3(), pv = new THREE.Vector3(), one = new THREE.Vector3(1, 1, 1);
    for (const Ln of lines) {
      const p = Ln.p, n = p.length / 3; if (n < 3 || nv + n * 4 > this.maxV) continue;
      const base = nv; let nx = 0, ny = 1, nz = 0, arc = 0, nextArrow = 1.3;
      for (let i = 0; i < n; i++) {
        const a = Math.max(0, i - 1), b = Math.min(n - 1, i + 1);
        let tx = p[b * 3] - p[a * 3], ty = p[b * 3 + 1] - p[a * 3 + 1], tz = p[b * 3 + 2] - p[a * 3 + 2];
        const tl = Math.hypot(tx, ty, tz) || 1; tx /= tl; ty /= tl; tz /= tl;
        // carry the normal along the line so the tube does not twist
        const d = nx * tx + ny * ty + nz * tz; nx -= d * tx; ny -= d * ty; nz -= d * tz;
        let nl = Math.hypot(nx, ny, nz);
        if (nl < 1e-3) { if (Math.abs(tx) < 0.9) { nx = 0; ny = -tz; nz = ty; } else { nx = tz; ny = 0; nz = -tx; } nl = Math.hypot(nx, ny, nz); }
        nx /= nl; ny /= nl; nz /= nl;
        const bx = ty * nz - tz * ny, by = tz * nx - tx * nz, bz = tx * ny - ty * nx;
        for (let k = 0; k < 4; k++) {
          const c = k === 0 ? 1 : k === 2 ? -1 : 0, s = k === 1 ? 1 : k === 3 ? -1 : 0, o = (nv + k) * 3;
          P[o] = p[i * 3] + r * (c * nx + s * bx); P[o + 1] = p[i * 3 + 1] + r * (c * ny + s * by); P[o + 2] = p[i * 3 + 2] + r * (c * nz + s * bz);
          lineColor(Ln.b[i], C, o);
        }
        nv += 4;
        if (i > 0) { arc += Math.hypot(p[i * 3] - p[i * 3 - 3], p[i * 3 + 1] - p[i * 3 - 2], p[i * 3 + 2] - p[i * 3 - 1]);
          if (arc > nextArrow && na < 900) { nextArrow += 2.6; pv.set(p[i * 3], p[i * 3 + 1], p[i * 3 + 2]); tv.set(tx, ty, tz); q.setFromUnitVectors(up, tv); M.compose(pv, q, one); this.arrows.setMatrixAt(na++, M); } }
      }
      const segs = Ln.closed ? n : n - 1;
      for (let i = 0; i < segs; i++) {
        const a0 = base + i * 4, b0 = base + ((i + 1) % n) * 4;
        for (let k = 0; k < 4; k++) { const k1 = (k + 1) % 4; X[ni++] = a0 + k; X[ni++] = b0 + k; X[ni++] = b0 + k1; X[ni++] = a0 + k; X[ni++] = b0 + k1; X[ni++] = a0 + k1; }
      }
    }
    this.geo.setDrawRange(0, ni);
    for (const at of ['position', 'color']) { const A = this.geo.attributes[at]; A.needsUpdate = true; A.clearUpdateRanges && A.clearUpdateRanges(); }
    this.geo.index.needsUpdate = true;
    this.arrows.count = na; this.arrows.instanceMatrix.needsUpdate = true;
  }
}

// magnet only lines, in the magnet's own frame
function buildRigidLines() {
  const a = MAGNET.a, b = MAGNET.b, F = (r, z) => magnetField(r, z, a, b, MAGNET.Br);
  const hFn = (r, z) => { const d = Math.hypot(Math.max(Math.abs(Math.abs(r) - a), 0), Math.max(Math.abs(z) - b, 0)); return Math.min(Math.max(0.05 * Math.max(d, 0.002), 1.5e-4), 0.02 * Math.hypot(r, z) + 1e-3, 4e-3); };
  const lines = seedRadii.map((r0, k) => {
    const L = traceMeridional(F, r0, a, hFn, 0.45), ang = k * GOLD, ca = Math.cos(ang), sa = Math.sin(ang);
    const pts = L.closed ? L.pts.slice(0, -1) : L.pts, p = [], bb = [];
    for (const [rho, z] of pts) { p.push(z * S, rho * ca * S, rho * sa * S); const f = F(rho, z); bb.push(Math.hypot(f[0], f[1])); }
    return { p, b: bb, closed: L.closed };
  });
  vis.rigid.set(lines);
}

// combined field in world coordinates
function Bworld(px, py, pz, sn, out) {
  const m = sn.m, dx = px - sn.ms, z = dx * m[0] + py * m[1] + pz * m[2];
  const rx = dx - z * m[0], ry = py - z * m[1], rz = pz - z * m[2], rho = Math.hypot(rx, ry, rz);
  const f = magnetField(rho, z, MAGNET.a, MAGNET.b, MAGNET.Br), ir = rho > 1e-12 ? f[0] / rho : 0;
  let bx = ir * rx + f[1] * m[0], by = ir * ry + f[1] * m[1], bz = ir * rz + f[1] * m[2];
  if (sn.NI !== 0) {
    const n = sn.n, zc = px * n[0] + py * n[1] + pz * n[2], cx = px - zc * n[0], cy = py - zc * n[1], cz = pz - zc * n[2], rc = Math.hypot(cx, cy, cz);
    let sr = 0, sz = 0; for (const L of sn.loops) { const g = loopB(rc, zc - L.z, L.r, sn.NI * L.w); sr += g[0]; sz += g[1]; }
    const ic = rc > 1e-12 ? sr / rc : 0; bx += ic * cx + sz * n[0]; by += ic * cy + sz * n[1]; bz += ic * cz + sz * n[2];
  }
  if (sn.J !== 0) {
    const n = sn.n, zc = px * n[0] + py * n[1] + pz * n[2], cx = px - zc * n[0], cy = py - zc * n[1], cz = pz - zc * n[2], rc = Math.hypot(cx, cy, cz);
    const g = magnetField(rc, zc, IRON.rc, IRON.h, sn.J), ic = rc > 1e-12 ? g[0] / rc : 0;
    bx += ic * cx + g[1] * n[0]; by += ic * cy + g[1] * n[1]; bz += ic * cz + g[1] * n[2];
  }
  out[0] = bx; out[1] = by; out[2] = bz;
}
function stepSize(px, py, pz, sn) {
  const m = sn.m, dx = px - sn.ms, z = dx * m[0] + py * m[1] + pz * m[2];
  const rho = Math.hypot(dx - z * m[0], py - z * m[1], pz - z * m[2]);
  let dm = Math.hypot(Math.max(rho - MAGNET.a, 0), Math.max(Math.abs(z) - MAGNET.b, 0)); if (dm === 0) dm = Math.min(MAGNET.a - rho, MAGNET.b - Math.abs(z));
  let d = Math.max(dm, 0.0015);
  if (sn.NI !== 0) { const n = sn.n, zc = px * n[0] + py * n[1] + pz * n[2], rc = Math.hypot(px - zc * n[0], py - zc * n[1], pz - zc * n[2]);
    d = Math.min(d, Math.max(Math.hypot(rc - coil.rMean, zc) - 0.6 * Math.hypot(coil.depth / 2, BOBBIN.halfLen), 0.002)); }
  if (sn.J !== 0) { const n = sn.n, zc = px * n[0] + py * n[1] + pz * n[2], rc = Math.hypot(px - zc * n[0], py - zc * n[1], pz - zc * n[2]);
    let dc = Math.hypot(Math.max(rc - IRON.rc, 0), Math.max(Math.abs(zc) - IRON.h, 0)); if (dc === 0) dc = Math.min(IRON.rc - rc, IRON.h - Math.abs(zc));
    d = Math.min(d, Math.max(dc, 0.0015)); }
  return Math.min(Math.max(0.05 * d, 3e-4), 0.03 * Math.hypot(dx, py, pz) + 1.5e-3, 5e-3);
}
function traceDir(x0, y0, z0, sign, sn, maxSteps) {
  const p = [x0, y0, z0], bb = [], B = [0, 0, 0], k = [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]];
  let x = x0, y = y0, z = z0, closed = false, throughMag = false, prevSide = 0;
  const dir = (ax, ay, az, o) => { Bworld(ax, ay, az, sn, B); const l = Math.hypot(B[0], B[1], B[2]) || 1e-30; o[0] = sign * B[0] / l; o[1] = sign * B[1] / l; o[2] = sign * B[2] / l; return l; };
  for (let i = 0; i < maxSteps; i++) {
    const h = stepSize(x, y, z, sn);
    const bmag = dir(x, y, z, k[0]); bb.push(bmag);
    dir(x + h / 2 * k[0][0], y + h / 2 * k[0][1], z + h / 2 * k[0][2], k[1]);
    dir(x + h / 2 * k[1][0], y + h / 2 * k[1][1], z + h / 2 * k[1][2], k[2]);
    dir(x + h * k[2][0], y + h * k[2][1], z + h * k[2][2], k[3]);
    x += h / 6 * (k[0][0] + 2 * k[1][0] + 2 * k[2][0] + k[3][0]); y += h / 6 * (k[0][1] + 2 * k[1][1] + 2 * k[2][1] + k[3][1]); z += h / 6 * (k[0][2] + 2 * k[1][2] + 2 * k[2][2] + k[3][2]);
    // did it pass through the magnet's midplane inside the magnet?
    const m = sn.m, dz = (x - sn.ms) * m[0] + y * m[1] + z * m[2], side = dz < 0 ? -1 : 1;
    if (prevSide && side !== prevSide && Math.hypot(x - sn.ms - dz * m[0], y - dz * m[1], z - dz * m[2]) < MAGNET.a) throughMag = true;
    prevSide = side;
    if (i > 25 && Math.hypot(x - x0, y - y0, z - z0) < 1.5 * h) { closed = true; break; }
    if (Math.hypot(x, y, z) > 0.32) { p.push(x, y, z); Bworld(x, y, z, sn, B); bb.push(Math.hypot(B[0], B[1], B[2])); break; }
    p.push(x, y, z);
  }
  if (bb.length < p.length / 3) { Bworld(x, y, z, sn, B); bb.push(Math.hypot(B[0], B[1], B[2])); }
  return { p, bb, closed, throughMag };
}
function traceSeed(seed, sn) {
  const f = traceDir(seed[0], seed[1], seed[2], 1, sn, 900);
  if (seed[3] === 'coil' && f.throughMag) return null;       // already drawn from a magnet seed
  let p = f.p, bb = f.bb;
  if (!f.closed) {
    const bk = traceDir(seed[0], seed[1], seed[2], -1, sn, 900), np = [], nb = [];
    for (let i = bk.p.length / 3 - 1; i >= 1; i--) { np.push(bk.p[i * 3], bk.p[i * 3 + 1], bk.p[i * 3 + 2]); nb.push(bk.bb[i]); }
    p = np.concat(p); bb = nb.concat(bb);
  }
  for (let i = 0; i < p.length; i++) p[i] *= S;
  return { p, b: bb, closed: f.closed };
}
// Strength map: |B| on the horizontal plane through the track, same colormap, alpha rising with strength.
class FieldMap {
  constructor(W, H, wx, wz) {
    this.W = W; this.H = H; this.wx = wx; this.wz = wz; this.data = new Uint8Array(W * H * 4);
    this.tex = new THREE.DataTexture(this.data, W, H, THREE.RGBAFormat); this.tex.magFilter = this.tex.minFilter = THREE.LinearFilter; this.tex.colorSpace = THREE.SRGBColorSpace;
    this.mesh = new THREE.Mesh(new THREE.PlaneGeometry(wx * S, wz * S).rotateX(-PI / 2), new THREE.MeshBasicMaterial({ map: this.tex, transparent: true, depthWrite: false, side: THREE.DoubleSide }));
    this.mesh.renderOrder = 1; this.row = 0; this.job = null;
  }
  // fill rows with fn(x, z) -> |B|; returns true when the whole map is done
  fill(fn, budgetMs) {
    const t0 = performance.now(), c = [0, 0, 0], D = this.data;
    while (this.row < this.H && performance.now() - t0 < budgetMs) {
      const j = this.row, z = (0.5 - (j + 0.5) / this.H) * this.wz; // texture row 0 sits at +z after the plane rotation
      for (let i = 0; i < this.W; i++) {
        const x = ((i + 0.5) / this.W - 0.5) * this.wx, t = bT(fn(x, z)), o = (j * this.W + i) * 4;
        ramp(t, c, 0); D[o] = c[0] * 255; D[o + 1] = c[1] * 255; D[o + 2] = c[2] * 255; D[o + 3] = 200 * Math.pow(t, 1.6);
      }
      this.row++;
    }
    if (this.row >= this.H) { this.tex.needsUpdate = true; return true; }
    return false;
  }
}
function buildRigidMap() {
  const m = vis.rigidMap; m.row = 0;
  m.fill((x, z) => { const f = magnetField(Math.abs(z), x, MAGNET.a, MAGNET.b, MAGNET.Br); return Math.hypot(f[0], f[1]); }, 1e9);
}
const FL = { job: null, doneKey: '', mode: 'rigid', mapKey: '', mapSn: null };
function snapshot() {
  const m = magAxis(state.psi), n = coilAxis(state.theta), NI = state.N * state.I, J = state.core ? state.J : 0;
  const key = [state.s.toFixed(4), state.psi.toFixed(3), state.theta.toFixed(3), (NI).toPrecision(2), state.N, J.toPrecision(2)].join('|');
  return { ms: state.s, m, n, NI, J, loops: coil.traceLoops, key };
}
function coilMomentRatio() { return Math.abs(state.N * state.I) * PI * coil.rMean * coil.rMean / magMoment(); }
function updateFieldLines() {
  if (state.exp === 'tube') { vis.rigid.mesh.visible = vis.rigid.arrows.visible = state.show.lines; vis.live.mesh.visible = vis.live.arrows.visible = false; vis.rigidMap.mesh.visible = state.show.map; vis.liveMap.mesh.visible = false; return; }
  const traced = coilMomentRatio() > 0.03 || Math.abs(state.V) >= 0.5 || (state.core && Math.abs(state.J) > 0.003);
  const L = state.show.lines, Mp = state.show.map;
  if (!traced) {
    if (FL.mode !== 'rigid') { FL.mode = 'rigid'; FL.job = null; FL.doneKey = ''; FL.mapKey = ''; }
    vis.rigid.mesh.visible = vis.rigid.arrows.visible = L; vis.live.mesh.visible = vis.live.arrows.visible = false;
    vis.rigidMap.mesh.visible = Mp; vis.liveMap.mesh.visible = false; return;
  }
  // live strength map, a few rows per frame
  if (Mp) {
    if (!FL.mapSn) { const sn = snapshot(); if (sn.key !== FL.mapKey) { FL.mapSn = sn; vis.liveMap.row = 0; } }
    if (FL.mapSn) { const sn = FL.mapSn, B = [0, 0, 0];
      if (vis.liveMap.fill((x, z) => { Bworld(x, 0, z, sn, B); return Math.hypot(B[0], B[1], B[2]); }, 4)) { FL.mapKey = sn.key; FL.mapSn = null; vis.liveMap.mesh.visible = true; vis.rigidMap.mesh.visible = false; } }
  } else { vis.liveMap.mesh.visible = vis.rigidMap.mesh.visible = false; }
  if (!FL.job) {
    const sn = snapshot(); if (sn.key === FL.doneKey && FL.mode === 'traced') return;
    const [u, v] = perpBasis(sn.m), seeds = [];
    seedRadii.forEach((r, k) => { const c = Math.cos(k * GOLD) * r, s = Math.sin(k * GOLD) * r; seeds.push([sn.ms + u[0] * c + v[0] * s, u[1] * c + v[1] * s, u[2] * c + v[2] * s, 'mag']); });
    const ratio = coilMomentRatio(), nc = Math.abs(state.V) >= 0.5 ? Math.round(6 + 10 * Math.min(1, ratio / 0.3)) : Math.min(16, Math.round(NL * ratio)), [cu, cv] = perpBasis(sn.n);
    for (let k = 0; k < nc; k++) { const r = BOBBIN.rIn * 0.92 * Math.sqrt((k + 0.5) / nc), a = k * GOLD + 0.4, c = Math.cos(a) * r, s = Math.sin(a) * r; seeds.push([cu[0] * c + cv[0] * s, cu[1] * c + cv[1] * s, cu[2] * c + cv[2] * s, 'coil']); }
    FL.job = { sn, seeds, i: 0, lines: [] };
  }
  const J = FL.job, t0 = performance.now();
  while (J.i < J.seeds.length && performance.now() - t0 < 7) { const L = traceSeed(J.seeds[J.i++], J.sn); if (L) J.lines.push(L); }
  if (J.i >= J.seeds.length) {
    vis.live.set(J.lines); FL.doneKey = J.sn.key; FL.job = null;
    if (FL.mode !== 'traced') { FL.mode = 'traced'; }
    vis.rigid.mesh.visible = vis.rigid.arrows.visible = false; vis.live.mesh.visible = vis.live.arrows.visible = true;
  }
  if (!L) { vis.live.mesh.visible = vis.live.arrows.visible = vis.rigid.mesh.visible = vis.rigid.arrows.visible = false; }
}

/* =====================================================================
   FALLING THROUGH A TUBE
   The tube wall is split into 120 rings, 2.5 mm tall. Each ring has its
   resistance, its self inductance, and a mutual inductance to every other
   ring (Maxwell's formula for coaxial loops). The magnet's flux through a
   ring comes from the same exact magnet model. Each substep solves, together,
     (L/h + R) I' = L I / h - G v'      (all ring currents)
     m v' = m v - h m g + h G . I'       (the magnet)
   with G_k = dPhi_k/dy. The field of the eddy currents is about one percent
   of the magnet's field near the wall, so the drawn lines show the magnet only.
   ===================================================================== */
const TUBES = {
  cu: { sym: 'Cu', nm: 'Copper', rho: 1.68e-8, color: 0xc27a45, metal: 1 },
  al: { sym: 'Al', nm: 'Aluminium', rho: 2.65e-8, color: 0xb9c0c8, metal: 1 },
  br: { sym: 'CuZn', nm: 'Brass', rho: 7.0e-8, color: 0xc8a44f, metal: 1 },
  ss: { sym: 'SS', nm: 'Stainless', rho: 7.2e-7, color: 0x9ea3a8, metal: 0.9 },
  pvc: { sym: 'PVC', nm: 'Plastic', rho: Infinity, color: 0xd9d8d0, metal: 0 }
};
const TUBE = { ri: 0.007, len: 0.30, hr: 0.0025, N: 120, top: 0.15, start: 0.19, floor: -0.215 };
const tube = { mat: 'cu', t: 0.0015, y: TUBE.start, v: 0, phase: 'held', I: new Float64Array(TUBE.N), dirty: true, tIn: null, tOut: null, tDrop: 0, ySpring: 0 };
const mutualLoops = (a, b, d) => { const m = 4 * a * b / ((a + b) * (a + b) + d * d), [K, E] = ellipKE(m), k = Math.sqrt(m); return MU0 * Math.sqrt(a * b) * ((2 / k - k) * K - 2 / k * E); };
function buildTube() {
  const N = TUBE.N, rm = TUBE.ri + tube.t / 2, T = TUBES[tube.mat];
  tube.rm = rm; tube.z = new Float64Array(N); for (let k = 0; k < N; k++) tube.z[k] = -TUBE.len / 2 + (k + 0.5) * TUBE.hr;
  tube.cond = isFinite(T.rho);
  tube.Rk = tube.cond ? T.rho * 2 * PI * rm / (tube.t * TUBE.hr) : Infinity;
  const gmd = 0.2235 * (tube.t + TUBE.hr); tube.Ls = MU0 * rm * (Math.log(8 * rm / gmd) - 2);   // thin rectangular ring
  const Mrow = new Float64Array(N); Mrow[0] = tube.Ls; for (let k = 1; k < N; k++) Mrow[k] = mutualLoops(rm, rm, k * TUBE.hr);
  tube.Mrow = Mrow;
  // magnet flux through a ring at axial offset d, as a table
  const gl = gaussLegendre(40), M = MAGNET.Br / MU0, D0 = -0.25, dz = 5e-4, P = new Float64Array(1001);
  for (let i = 0; i <= 1000; i++) { const d = D0 + i * dz; let f = 0; for (let j = 0; j < 40; j++) f += gl.w[j] * MAGNET.b * M * mutualLoops(MAGNET.a, rm, d - gl.x[j] * MAGNET.b); P[i] = f; }
  tube.P = { vals: P, x0: D0, h: dz, N: 1001 };
  // (L/h + R) inverse for the substep length used below
  const h = 1 / 480; tube.h = h;
  if (tube.cond) {
    const A = new Float64Array(N * N), Ai = new Float64Array(N * N);
    for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) { A[i * N + j] = Mrow[Math.abs(i - j)] / h + (i === j ? tube.Rk : 0); Ai[i * N + j] = i === j ? 1 : 0; }
    for (let c = 0; c < N; c++) {   // Gauss Jordan; the matrix is symmetric and diagonally dominant
      const d = A[c * N + c]; for (let j = 0; j < N; j++) { A[c * N + j] /= d; Ai[c * N + j] /= d; }
      for (let r = 0; r < N; r++) if (r !== c) { const f = A[r * N + c]; if (f) for (let j = 0; j < N; j++) { A[r * N + j] -= f * A[c * N + j]; Ai[r * N + j] -= f * Ai[c * N + j]; } }
    }
    tube.Ai = Ai;
  }
  // quasi static terminal speed at mid tube: m g = v * sum G^2 / R
  let sg = 0; if (tube.cond) for (let k = 0; k < N; k++) { const g = tabEval(tube.z[k], tube.P); if (g) sg += g[1] * g[1] / tube.Rk; }
  tube.vTerm = tube.cond ? magnetMass() * 9.81 / sg : Infinity;
  tube.I.fill(0); tube.dirty = false;
}
function dropMagnet() {
  if (tube.dirty) buildTube();
  tube.y = TUBE.start; tube.v = 0; tube.I.fill(0); tube.phase = 'fall'; tube.tIn = tube.tOut = null; tube.tDrop = state.t; scope.buf.length = 0;
}
function simulateTube(dt) {
  if (tube.dirty) buildTube();
  const N = TUBE.N, m = magnetMass(), g = 9.81, n = 8, h = dt / n, sig = tube.up ? 1 : -1;
  const G = tube.G || (tube.G = new Float64Array(N)), rhs = tube.rhs || (tube.rhs = new Float64Array(N)), w = tube.w || (tube.w = new Float64Array(N)), u = tube.u || (tube.u = new Float64Array(N));
  const I = tube.I, Mr = tube.Mrow, Ai = tube.Ai, samples = [];
  const sT0 = state.sTargetPrev, sT1 = state.sTarget; state.sTargetPrev = sT1;
  if (state.dragging) state.vTarget += ((sT1 - sT0) / dt - state.vTarget) * Math.min(1, dt / 0.05); else state.vTarget = 0;
  let F = 0;
  for (let s = 0; s < n; s++) {
    let y = tube.y, v = tube.v;
    if (tube.phase === 'held' || tube.phase === 'landed') { v = 0; }
    // magnet flux slope for every ring at the midpoint of this substep
    const ym = y + v * h / 2;
    for (let k = 0; k < N; k++) { const r = tabEval(tube.z[k] - ym, tube.P); G[k] = r ? -sig * r[1] : 0; }
    let ext = -m * g, Gw = 0, Gu = 0;
    if (tube.phase === 'drag') { const w0 = 30, yT = sT0 + (sT1 - sT0) * s / n; ext = m * w0 * w0 * (yT - y) + 2 * m * w0 * (state.vTarget - v); }   // your hand also carries the weight
    if (tube.cond && tube.phase !== 'held' && tube.phase !== 'landed') {
      for (let i = 0; i < N; i++) { let a = 0; for (let j = 0; j < N; j++) a += Mr[i > j ? i - j : j - i] * I[j]; rhs[i] = a / h; }
      for (let i = 0; i < N; i++) { let a = 0, b = 0; const o = i * N; for (let j = 0; j < N; j++) { a += Ai[o + j] * rhs[j]; b += Ai[o + j] * G[j]; } w[i] = a; u[i] = b; Gw += G[i] * a; Gu += G[i] * b; }
      const vn = (m * v + h * ext + h * Gw) / (m + h * Gu);
      F = 0; for (let i = 0; i < N; i++) { I[i] = w[i] - u[i] * vn; F += G[i] * I[i]; }
      v = vn;
    } else if (tube.phase !== 'held' && tube.phase !== 'landed') { v += h * ext / m; F = 0; I.fill(0); }
    else { if (tube.cond) for (let i = 0; i < N; i++) I[i] *= 0.5; F = 0; }
    y += v * h;
    if (y <= TUBE.floor + MAGNET.b && tube.phase === 'fall') { y = TUBE.floor + MAGNET.b; v = 0; tube.phase = 'landed'; }
    if (tube.phase === 'drag') y = Math.max(TUBE.floor + MAGNET.b, Math.min(0.26, y));
    // transit timing through the tube
    if (tube.phase === 'fall') { if (tube.tIn === null && y - MAGNET.b < TUBE.top) tube.tIn = state.t + (s + 1) * h; if (tube.tOut === null && y + MAGNET.b < -TUBE.top) tube.tOut = state.t + (s + 1) * h; }
    tube.y = y; tube.v = v;
    samples.push([state.t + (s + 1) * h, -v, F]);
  }
  tube.F = F;
  let heat = 0, imax = 0; for (let i = 0; i < N; i++) { heat += I[i] * I[i] * tube.Rk; if (Math.abs(I[i]) > Math.abs(imax)) imax = I[i]; }
  tube.heat = tube.cond ? heat : 0; tube.imax = imax;
  // the ring just below the magnet's lower end, for the meter
  const kb = Math.max(0, Math.min(N - 1, Math.floor((tube.y - MAGNET.b - 0.004 + TUBE.len / 2) / TUBE.hr)));
  tube.iBelow = I[kb];
  state.samples = samples; state.t += dt;
}

/* =====================================================================
   SCENE
   ===================================================================== */
let renderer, scene, camera, magnetGroup, coilGroup, tubeGroup;
const vis = {};
function textSprite(text, color, size) {
  const c = document.createElement('canvas'); c.width = 128; c.height = 128; const g = c.getContext('2d');
  g.font = '600 88px Archivo, Arial, sans-serif'; g.textAlign = 'center'; g.textBaseline = 'middle'; g.fillStyle = color; g.fillText(text, 64, 70);
  const tex = new THREE.CanvasTexture(c); tex.colorSpace = THREE.SRGBColorSpace;
  const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: tex, transparent: true, depthTest: false })); sp.scale.set(size, size, 1); sp.renderOrder = 10; return sp;
}
function labelSprite(text, size, color = '#9a968c') {
  const c = document.createElement('canvas'), g = c.getContext('2d'), font = '500 44px Archivo, Arial, sans-serif'; g.font = font;
  const w = Math.ceil(g.measureText(text).width) + 24; c.width = w; c.height = 64; g.font = font; g.fillStyle = color; g.textBaseline = 'middle'; g.fillText(text, 12, 34);
  const tex = new THREE.CanvasTexture(c); tex.colorSpace = THREE.SRGBColorSpace;
  const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: tex, transparent: true, depthTest: false })); sp.scale.set(size * w / 64, size, 1); sp.renderOrder = 10; return sp;
}
function envMap() {
  const pm = new THREE.PMREMGenerator(renderer), es = new THREE.Scene();
  es.add(new THREE.Mesh(new THREE.SphereGeometry(50, 32, 16), new THREE.ShaderMaterial({ side: THREE.BackSide,
    vertexShader: 'varying vec3 p; void main(){ p = position; gl_Position = projectionMatrix*modelViewMatrix*vec4(position,1.); }',
    fragmentShader: 'varying vec3 p; void main(){ float h = normalize(p).y; gl_FragColor = vec4(mix(vec3(0.06,0.06,0.058), vec3(0.38,0.37,0.35), smoothstep(-0.3,0.8,h)),1.); }' })));
  const panel = (x, y, z, w, hh, i) => { const m = new THREE.Mesh(new THREE.PlaneGeometry(w, hh), new THREE.MeshBasicMaterial({ color: new THREE.Color(i, i, i * 1.03), side: THREE.DoubleSide })); m.position.set(x, y, z); m.lookAt(0, 0, 0); es.add(m); };
  panel(0, 30, 10, 40, 14, 3.2); panel(-30, 12, -12, 14, 20, 1.4); panel(28, 8, 18, 10, 16, 1.0);
  const t = pm.fromScene(es, 0.03).texture; pm.dispose(); return t;
}
function initScene() {
  const view = $('view3d');
  renderer = new THREE.WebGLRenderer({ antialias: true, powerPreference: 'high-performance' });
  renderer.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
  renderer.toneMapping = THREE.ACESFilmicToneMapping;
  view.insertBefore(renderer.domElement, view.firstChild);
  renderer.domElement.setAttribute('aria-label', 'Three dimensional view of the magnet, the coil and the magnetic field');
  scene = new THREE.Scene(); scene.background = new THREE.Color(0x141414); scene.environment = envMap();
  camera = new THREE.PerspectiveCamera(36, 1, 0.05, 600);
  scene.add(new THREE.HemisphereLight(0xdedbd4, 0x1c1b19, 0.55));
  const key = new THREE.DirectionalLight(0xffffff, 1.6); key.position.set(8, 14, 10); scene.add(key);
  const rim = new THREE.DirectionalLight(0xe8e2d6, 0.6); rim.position.set(-10, 4, -12); scene.add(rim);
  const benchY = -4.6, grid = new THREE.GridHelper(60, 60, 0x33322f, 0x222220); grid.position.y = benchY; scene.add(grid);
  const tick = []; for (let x = -20; x <= 20; x++) { const l = x % 5 === 0 ? 0.9 : 0.4; tick.push(x, benchY + 0.01, 2.2, x, benchY + 0.01, 2.2 + l); } tick.push(-20, benchY + 0.01, 2.2, 20, benchY + 0.01, 2.2);
  const tg = new THREE.BufferGeometry(); tg.setAttribute('position', new THREE.Float32BufferAttribute(tick, 3)); scene.add(new THREE.LineSegments(tg, new THREE.LineBasicMaterial({ color: 0x6d6a63 })));
  for (let x = -20; x <= 20; x += 5) { const l = labelSprite(x + ' cm', 0.55); l.position.set(x, benchY + 0.05, 3.7); l.material.depthTest = true; scene.add(l); }

  magnetGroup = new THREE.Group(); scene.add(magnetGroup);
  const half = (color, sign) => {
    const g = new THREE.CylinderGeometry(MAGNET.a * S, MAGNET.a * S, MAGNET.b * S, 48, 1); g.rotateZ(-PI / 2); g.translate(sign * MAGNET.b * S / 2, 0, 0);
    const mesh = new THREE.Mesh(g, new THREE.MeshPhysicalMaterial({ color, metalness: 0.15, roughness: 0.32, clearcoat: 0.6, clearcoatRoughness: 0.25, transparent: true, opacity: 0.62, depthWrite: false }));
    mesh.renderOrder = 2; magnetGroup.add(mesh);
  };
  half(0xc9463d, 1); half(0x3b6fc9, -1);
  const cap = (txt, sign) => {
    const c = document.createElement('canvas'); c.width = c.height = 128; const g = c.getContext('2d');
    g.font = '700 84px Archivo, Arial, sans-serif'; g.textAlign = 'center'; g.textBaseline = 'middle'; g.fillStyle = 'rgba(255,255,255,0.92)'; g.fillText(txt, 64, 70);
    const t = new THREE.CanvasTexture(c); t.colorSpace = THREE.SRGBColorSpace;
    const m = new THREE.Mesh(new THREE.CircleGeometry(MAGNET.a * S * 0.8, 32), new THREE.MeshBasicMaterial({ map: t, transparent: true, depthWrite: false }));
    m.rotation.y = sign * PI / 2; m.position.x = sign * (MAGNET.b * S + 0.01); m.renderOrder = 3; magnetGroup.add(m);
  };
  cap('N', 1); cap('S', -1);
  const nS = textSprite('N', '#ffd9d4', 0.7); nS.position.set(MAGNET.b * S / 2, MAGNET.a * S + 0.6, 0);
  const sS = textSprite('S', '#d6e3ff', 0.7); sS.position.set(-MAGNET.b * S / 2, MAGNET.a * S + 0.6, 0);
  magnetGroup.add(nS, sS); vis.magLabels = [nS, sS];
  vis.pick = new THREE.Mesh(new THREE.CylinderGeometry(MAGNET.a * S * 1.7, MAGNET.a * S * 1.7, MAGNET.b * S * 2.3, 16).rotateZ(-PI / 2), new THREE.MeshBasicMaterial({ visible: false }));
  magnetGroup.add(vis.pick);
  // force on the magnet from the coil
  const fa = new THREE.Group(), fm = new THREE.MeshBasicMaterial({ color: 0xf2f5f8 });
  const shaft = new THREE.Mesh(new THREE.CylinderGeometry(0.06, 0.06, 1, 8).rotateZ(-PI / 2).translate(0.5, 0, 0), fm);
  const head = new THREE.Mesh(new THREE.ConeGeometry(0.2, 0.5, 12).rotateZ(-PI / 2), fm);
  fa.add(shaft, head); scene.add(fa); vis.force = { g: fa, shaft, head };
  vis.forceLabel = { push: labelSprite('push', 0.5, '#f2f5f8'), pull: labelSprite('pull', 0.5, '#f2f5f8') }; scene.add(vis.forceLabel.push, vis.forceLabel.pull);

  // field lines: magnet only set rides with the magnet, live set lives in the world
  vis.rigid = new TubeSet(40000, 0.03); magnetGroup.add(vis.rigid.mesh, vis.rigid.arrows);
  vis.live = new TubeSet(70000, 0.03); scene.add(vis.live.mesh, vis.live.arrows); vis.live.mesh.visible = vis.live.arrows.visible = false;
  buildRigidLines();
  vis.rigidMap = new FieldMap(256, 150, 0.44, 0.26); magnetGroup.add(vis.rigidMap.mesh); buildRigidMap();
  vis.liveMap = new FieldMap(150, 88, 0.46, 0.27); scene.add(vis.liveMap.mesh); vis.liveMap.mesh.visible = false;

  coilGroup = new THREE.Group(); scene.add(coilGroup);
  rebuildCoilVisuals();
  vis.bench = [grid]; scene.children.forEach(o => { if (o !== grid && (o.isLineSegments || (o.isSprite && o.position.y < -4))) vis.bench.push(o); });
  tubeGroup = new THREE.Group(); tubeGroup.visible = false; scene.add(tubeGroup);
  rebuildTubeVisuals();
}
function rebuildCoilVisuals() {
  while (coilGroup.children.length) { const o = coilGroup.children.pop(); o.traverse(c => { if (c.geometry) c.geometry.dispose(); }); }
  const rIn = BOBBIN.rIn * S, rOut = coil.rOut * S, hl = BOBBIN.halfLen * S, fl = BOBBIN.flange * S, bore = BOBBIN.bore * S;
  const wind = new THREE.LatheGeometry([new THREE.Vector2(rIn, -hl), new THREE.Vector2(rOut, -hl), new THREE.Vector2(rOut, hl), new THREE.Vector2(rIn, hl)], 96).rotateZ(-PI / 2);
  const c = document.createElement('canvas'); c.width = 8; c.height = 256; const g = c.getContext('2d'), turns = Math.max(4, coil.perLayer);
  for (let i = 0; i < 256; i++) { const vv = 150 + 80 * Math.sin((i / 256 * turns * 3) % 1 * PI); g.fillStyle = `rgb(${vv},${vv},${vv})`; g.fillRect(0, i, 8, 1); }
  const bump = new THREE.CanvasTexture(c); bump.wrapS = bump.wrapT = THREE.RepeatWrapping;
  coilGroup.add(new THREE.Mesh(wind, new THREE.MeshPhysicalMaterial({ color: WIRES[state.wire].color, metalness: WIRES[state.wire].metal, roughness: 0.33, bumpMap: bump, bumpScale: 0.6, clearcoat: 0.4, clearcoatRoughness: 0.3 })));
  const acrylic = new THREE.MeshPhysicalMaterial({ color: 0xdfe8f2, roughness: 0.12, transparent: true, opacity: 0.2, depthWrite: false, side: THREE.DoubleSide });
  const bp = [[bore, -hl - fl], [rOut + 0.15, -hl - fl], [rOut + 0.15, -hl], [rIn, -hl], [rIn, hl], [rOut + 0.15, hl], [rOut + 0.15, hl + fl], [bore, hl + fl], [bore, -hl - fl]].map(q => new THREE.Vector2(q[0], q[1]));
  const bob = new THREE.Mesh(new THREE.LatheGeometry(bp, 96).rotateZ(-PI / 2), acrylic); bob.renderOrder = 3; coilGroup.add(bob);
  // glow: how much of the magnet's field passes through the coil
  vis.glowMat = new THREE.MeshBasicMaterial({ color: 0x35b779, transparent: true, opacity: 0, side: THREE.DoubleSide, depthWrite: false, blending: THREE.AdditiveBlending });
  const glow = new THREE.Mesh(new THREE.CircleGeometry(bore * 0.98, 64).rotateY(PI / 2), vis.glowMat); vis.glow = glow; glow.renderOrder = 4; coilGroup.add(glow);
  // current direction chevrons around the outside of the winding
  vis.chev = new THREE.InstancedMesh(new THREE.ConeGeometry(0.13, 0.4, 10), new THREE.MeshBasicMaterial({ color: 0xffd9b0 }), 24);
  vis.chev.frustumCulled = false; coilGroup.add(vis.chev); vis.chevPhase = vis.chevPhase || 0; vis.chevR = rOut + 0.32;
  vis.coilPoleP = textSprite('N', '#ffd9d4', 0.85); vis.coilPoleM = textSprite('S', '#d6e3ff', 0.85); coilGroup.add(vis.coilPoleP, vis.coilPoleM);
  vis.coilPoleX = Math.max(hl + fl, IRON.h * S * (state.core ? 1 : 0)) + 0.75;
  // laminated soft iron rod
  const lc = document.createElement('canvas'); lc.width = 256; lc.height = 8; const lg = lc.getContext('2d');
  for (let i = 0; i < 256; i++) { const vv = (i % 8 < 1) ? 90 : 170 + 20 * Math.sin(i * 0.9); lg.fillStyle = `rgb(${vv},${vv},${vv})`; lg.fillRect(i, 0, 1, 8); }
  const lt = new THREE.CanvasTexture(lc); lt.wrapS = lt.wrapT = THREE.RepeatWrapping;
  vis.coreMesh = new THREE.Mesh(new THREE.CylinderGeometry(IRON.rc * S, IRON.rc * S, 2 * IRON.h * S, 48, 1).rotateZ(-PI / 2),
    new THREE.MeshPhysicalMaterial({ color: 0x6b6f74, metalness: 0.85, roughness: 0.42, bumpMap: lt, bumpScale: 0.4 }));
  vis.coreMesh.visible = state.core; coilGroup.add(vis.coreMesh);
  vis.glowRef = Math.abs(innerTurnFlux(0, 0, 0));
}

function rebuildTubeVisuals() {
  while (tubeGroup.children.length) { const o = tubeGroup.children.pop(); o.traverse(c => { if (c.geometry) c.geometry.dispose(); }); }
  const T = TUBES[tube.mat], ri = TUBE.ri * S, ro = (TUBE.ri + tube.t) * S, hl = TUBE.len / 2 * S;
  // cut away: only the back half of the wall is drawn, so you can see inside
  const prof = [[ri, -hl], [ro, -hl], [ro, hl], [ri, hl], [ri, -hl]].map(q => new THREE.Vector2(q[0], q[1]));
  const mat = new THREE.MeshPhysicalMaterial({ color: T.color, metalness: T.metal, roughness: T.metal ? 0.35 : 0.6, side: THREE.DoubleSide, transparent: !T.metal, opacity: T.metal ? 1 : 0.55 });
  tubeGroup.add(new THREE.Mesh(new THREE.LatheGeometry(prof, 64, PI, PI), mat));
  // eddy current bands on the inner wall, one per ring
  const pts = []; for (let i = 0; i <= 32; i++) { const t = PI + PI * i / 32; pts.push(new THREE.Vector3((ri - 0.03) * Math.sin(t), 0, (ri - 0.03) * Math.cos(t))); }
  const band = new THREE.TubeGeometry(new THREE.CatmullRomCurve3(pts), 32, TUBE.hr * S * 0.36, 4, false);
  vis.bands = new THREE.InstancedMesh(band, new THREE.MeshBasicMaterial({ color: 0xffffff }), TUBE.N); vis.bands.frustumCulled = false;
  const M = new THREE.Matrix4(), c0 = new THREE.Color(0x2a2926);
  for (let k = 0; k < TUBE.N; k++) { M.makeTranslation(0, (-TUBE.len / 2 + (k + 0.5) * TUBE.hr) * S, 0); vis.bands.setMatrixAt(k, M); vis.bands.setColorAt(k, c0); }
  tubeGroup.add(vis.bands);
  vis.ringArrows = new THREE.InstancedMesh(new THREE.ConeGeometry(0.13, 0.36, 10), new THREE.MeshBasicMaterial({ color: 0xf2efe8 }), 60); vis.ringArrows.count = 0; vis.ringArrows.frustumCulled = false; tubeGroup.add(vis.ringArrows);
  // vertical ruler beside the tube, 1 cm ticks
  const tk = [], zr = -(ro + 1.1);
  for (let y = -21; y <= 20; y++) { const l = y % 5 === 0 ? 0.8 : 0.35; tk.push(0, y, zr, 0, y, zr - l); } tk.push(0, -21, zr, 0, 20, zr);
  const tg = new THREE.BufferGeometry(); tg.setAttribute('position', new THREE.Float32BufferAttribute(tk, 3)); tubeGroup.add(new THREE.LineSegments(tg, new THREE.LineBasicMaterial({ color: 0x6d6a63 })));
  for (let y = -20; y <= 20; y += 5) { const l = labelSprite(y + ' cm', 0.55); l.position.set(0, y, zr - 1.9); l.material.depthTest = true; tubeGroup.add(l); }
  const pad = new THREE.Mesh(new THREE.CylinderGeometry(2.2, 2.2, 0.4, 40), new THREE.MeshStandardMaterial({ color: 0x2b2a27, roughness: 0.9 })); pad.position.y = TUBE.floor * S - 0.2; tubeGroup.add(pad);
  // free body diagram on the magnet: weight down, magnetic drag up, drawn to the same scale
  const arrow = (color, label) => { const g = new THREE.Group(), m = new THREE.MeshBasicMaterial({ color });
    const sh = new THREE.Mesh(new THREE.CylinderGeometry(0.06, 0.06, 1, 8).translate(0, 0.5, 0), m), hd = new THREE.Mesh(new THREE.ConeGeometry(0.2, 0.5, 12), m);
    const lb = labelSprite(label, 0.5, '#e9e6df'); g.add(sh, hd, lb); tubeGroup.add(g); return { g, sh, hd, lb }; };
  vis.fbd = { w: arrow(0xe9e6df, 'weight'), d: arrow(0x35b779, 'drag') };
}
function updateTubeVisuals(dt) {
  coilGroup.visible = false; tubeGroup.visible = true; vis.bench.forEach(o => o.visible = false); vis.force.g.visible = vis.forceLabel.push.visible = vis.forceLabel.pull.visible = false;
  magnetGroup.position.set(0, tube.y * S, 0); magnetGroup.rotation.set(0, 0, tube.up ? PI / 2 : -PI / 2);
  // the camera follows the magnet down the tube
  const ty = Math.max(TUBE.floor * S + 4, Math.min(TUBE.start * S, tube.y * S)); cam.target.y += (ty - cam.target.y) * Math.min(1, dt * 4);
  updateFieldLines();
  // band colors: yellow for counterclockwise seen from above, violet for clockwise, brightness with current
  const imax = Math.max(1e-6, ...Array.from(tube.I, Math.abs)), c = new THREE.Color(), cp = new THREE.Color(0xfde725), cn = new THREE.Color(0xb388ff), c0 = new THREE.Color(0x2a2926);
  const M = new THREE.Matrix4(), q = new THREE.Quaternion(), p = new THREE.Vector3(), one = new THREE.Vector3(1, 1, 1), upv = new THREE.Vector3(0, 1, 0), dir = new THREE.Vector3();
  let na = 0; const ref = Math.max(imax, 0.2);
  for (let k = 0; k < TUBE.N; k++) {
    const i = tube.I[k], f = Math.pow(Math.min(1, Math.abs(i) / ref), 0.6);
    c.copy(c0).lerp(i >= 0 ? cp : cn, f); vis.bands.setColorAt(k, c);
    if (k % 3 === 1 && Math.abs(i) > 0.3 * imax && imax > 0.05 && na < 60) {
      p.set(-(TUBE.ri * S - 0.25), (-TUBE.len / 2 + (k + 0.5) * TUBE.hr) * S, 0); dir.set(0, 0, Math.sign(i)); q.setFromUnitVectors(upv, dir); M.compose(p, q, one); vis.ringArrows.setMatrixAt(na++, M);
    }
  }
  vis.bands.instanceColor.needsUpdate = true; vis.ringArrows.count = na; vis.ringArrows.instanceMatrix.needsUpdate = true;
  const mg = magnetMass() * 9.81, len = (f) => 2.2 * f / mg, yb = tube.y * S, x0 = 3.4;
  const put = (a, f, down) => { const L = Math.max(0.01, len(Math.abs(f))); a.g.visible = Math.abs(f) > mg * 0.02 || down;
    a.sh.scale.set(1, L, 1); a.hd.position.y = L + 0.2; a.g.position.set(0, down ? yb - MAGNET.b * S - 0.2 : yb + MAGNET.b * S + 0.2, x0); a.g.rotation.set(down ? PI : 0, 0, 0); a.lb.position.set(0, L + 0.9, 0); };
  put(vis.fbd.w, tube.phase === 'drag' ? 0 : mg, true); put(vis.fbd.d, Math.max(0, tube.F), false);
  vis.fbd.w.g.visible = tube.phase !== 'drag';
}
function updateVisuals(dt) {
  if (state.exp === 'tube') { updateTubeVisuals(dt); return; }
  coilGroup.visible = true; tubeGroup.visible = false; vis.bench.forEach(o => o.visible = true); magnetGroup.rotation.set(0, 0, 0);
  vis.coreMesh.visible = state.core;
  magnetGroup.position.set(state.s * S, 0, 0); magnetGroup.rotation.y = state.psi; coilGroup.rotation.y = state.theta;
  updateFieldLines();
  const phi = innerTurnFlux(state.s, state.psi, state.theta);
  vis.glowMat.opacity = 0.03 + 0.42 * Math.pow(Math.min(1, Math.abs(phi) / vis.glowRef), 0.7);
  const vd = new THREE.Vector3(); camera.getWorldDirection(vd);
  const fade = (d) => 1 - Math.min(1, Math.max(0, (d - 0.75) / 0.17));
  const alongM = Math.abs(vd.x * Math.cos(state.psi) - vd.z * Math.sin(state.psi)), alongC = Math.abs(vd.x * Math.cos(state.theta) - vd.z * Math.sin(state.theta));
  vis.magLabels.forEach(l => l.material.opacity = fade(alongM));
  // coil poles: positive current makes north on the +x face of the coil
  const Ia = Math.abs(state.I), on = Ia > 1e-3, nC = coilAxis(state.theta);
  vis.coilPoleP.visible = vis.coilPoleM.visible = Math.abs(state.V) >= 0.5 || coilMomentRatio() > 0.02;
  vis.coilPoleP.position.set(state.I >= 0 ? vis.coilPoleX : -vis.coilPoleX, 0.9, coil.rOut * S + 0.2);
  vis.coilPoleM.position.set(state.I >= 0 ? -vis.coilPoleX : vis.coilPoleX, 0.9, coil.rOut * S + 0.2);
  for (const l of [vis.coilPoleP, vis.coilPoleM]) { const sx = Math.sign(l.position.x); l.material.opacity = sx * (nC[0] * -vd.x + nC[2] * -vd.z) > 0 ? 1 : fade(alongC); }
  vis.chev.visible = on;
  if (on) {
    const sg = Math.sign(state.I); vis.chevPhase += sg * Math.min(3.5, 0.5 + 2 * Math.log10(1 + Ia * 1000) / 3) * dt;
    const M = new THREE.Matrix4(), q = new THREE.Quaternion(), p = new THREE.Vector3(), one = new THREE.Vector3(1, 1, 1), tn = new THREE.Vector3(), up = new THREE.Vector3(0, 1, 0);
    for (let i = 0; i < 24; i++) { const t = vis.chevPhase + 2 * PI * i / 24; p.set(0, vis.chevR * Math.cos(t), vis.chevR * Math.sin(t)); tn.set(0, -Math.sin(t) * sg, Math.cos(t) * sg); q.setFromUnitVectors(up, tn); M.compose(p, q, one); vis.chev.setMatrixAt(i, M); }
    vis.chev.instanceMatrix.needsUpdate = true;
  }
  // force arrow: length grows with the log of the force
  const F = state.force, fa = vis.force, show = Math.abs(F) > 0.005;
  fa.g.visible = show; vis.forceLabel.push.visible = vis.forceLabel.pull.visible = false;
  if (show) {
    const len = 0.8 + 1.1 * Math.log10(Math.abs(F) / 0.005), dir = Math.sign(F), x0 = state.s * S + dir * (MAGNET.b * S + 0.5);
    fa.shaft.scale.set(len, 1, 1); fa.head.position.set(len + 0.2, 0, 0);
    fa.g.position.set(x0, 0, 0); fa.g.rotation.set(0, dir > 0 ? 0 : PI, 0);
    const lbl = dir * -Math.sign(state.s || -1) > 0 ? vis.forceLabel.pull : vis.forceLabel.push;
    lbl.visible = Math.abs(F) > 0.02; lbl.position.set(x0 + dir * (len * 0.5 + 0.3), -1.1, 0);
  }
}

/* =====================================================================
   CAMERA and POINTER
   ===================================================================== */
const cam = { az: 0.62, el: 0.36, r: 32, target: new THREE.Vector3(0, 0, 0), anim: null };
function placeCamera() {
  const t = cam.target, r = cam.r * Math.max(1, Math.pow(1.3 / camera.aspect, 1.0));
  camera.position.set(t.x + r * Math.cos(cam.el) * Math.sin(cam.az), t.y + r * Math.sin(cam.el), t.z + r * Math.cos(cam.el) * Math.cos(cam.az)); camera.lookAt(t);
}
function flyTo(az, el, r) {
  let d = ((az - cam.az + PI) % (2 * PI) + 2 * PI) % (2 * PI) - PI;
  cam.anim = { from: { az: cam.az, el: cam.el, r: cam.r }, to: { az: cam.az + d, el, r }, t: 0, dur: window.matchMedia('(prefers-reduced-motion: reduce)').matches ? 0.01 : 0.8 };
}
function setView(name) {
  document.querySelectorAll('.views button').forEach(b => b.setAttribute('aria-pressed', b.dataset.view === name));
  if (state.exp === 'tube') {
    cam.target.set(0, -1, 0);
    if (name === 'oblique') flyTo(PI / 2 + 0.75, 0.2, 30);
    if (name === 'side') flyTo(PI / 2, 0.02, 30);
    if (name === 'full') flyTo(PI / 2 + 0.5, 0.12, 62);
    if (name === 'top' || name === 'axis') flyTo(PI / 2, 1.5, 50);
    return;
  }
  cam.target.set(0, 0, 0);
  if (name === 'oblique') flyTo(0.62, 0.36, 32);
  if (name === 'side') flyTo(0, 0.05, 27);
  if (name === 'top') flyTo(0, 1.52, 30);
  if (name === 'axis') { const n = coilAxis(state.theta), side = state.s * n[0] > 0 ? 1 : -1; flyTo(Math.atan2(side * n[0], side * n[2]), 0.08, 20); }
}
function initPointer() {
  const el = renderer.domElement, ray = new THREE.Raycaster(), ndc = new THREE.Vector2(), pointers = new Map();
  let mode = null, last = null, pinch = 0;
  const rayAt = (e) => { const r = el.getBoundingClientRect(); ndc.set((e.clientX - r.left) / r.width * 2 - 1, -(e.clientY - r.top) / r.height * 2 + 1); ray.setFromCamera(ndc, camera); return ray.ray; };
  const trackPoint = (rr) => { const o = rr.origin, d = rr.direction, c = state.exp === 'tube' ? 'y' : 'x', a = 1 - d[c] * d[c]; return a < 1e-6 ? null : (o[c] - d[c] * o.dot(d)) / a / S; };
  el.addEventListener('pointerdown', (e) => {
    el.setPointerCapture(e.pointerId); pointers.set(e.pointerId, { x: e.clientX, y: e.clientY });
    if (pointers.size === 2) { const [p, q] = [...pointers.values()]; pinch = Math.hypot(p.x - q.x, p.y - q.y); mode = 'pinch'; return; }
    rayAt(e);
    if (ray.intersectObject(vis.pick, false).length && !state.flip) {
      mode = 'magnet'; el.classList.add('dragging'); state.dragging = true; state.vTarget = 0;
      if (state.exp === 'tube') { tube.phase = 'drag'; state.sTarget = state.sTargetPrev = tube.y; const s = trackPoint(rayAt(e)); vis.dragOffset = s === null ? 0 : tube.y - s; }
      else { setSwing(false); state.mode = 'drag'; state.sTarget = state.sTargetPrev = state.s; const s = trackPoint(rayAt(e)); vis.dragOffset = s === null ? 0 : state.s - s; }
    } else { mode = 'orbit'; cam.anim = null; }
    last = { x: e.clientX, y: e.clientY };
  });
  el.addEventListener('pointermove', (e) => {
    if (!pointers.has(e.pointerId)) return; pointers.set(e.pointerId, { x: e.clientX, y: e.clientY });
    if (mode === 'pinch' && pointers.size === 2) { const [p, q] = [...pointers.values()], d = Math.hypot(p.x - q.x, p.y - q.y); cam.r = Math.min(80, Math.max(6, cam.r * pinch / d)); pinch = d; return; }
    if (mode === 'magnet') { const s = trackPoint(rayAt(e)); if (s !== null) state.sTarget = state.exp === 'tube' ? Math.max(-0.2, Math.min(0.26, s + vis.dragOffset)) : Math.max(-TRACK, Math.min(TRACK, s + vis.dragOffset)); }
    else if (mode === 'orbit' && last) { cam.az -= (e.clientX - last.x) * 0.006; cam.el = Math.max(-1.45, Math.min(1.53, cam.el + (e.clientY - last.y) * 0.005)); document.querySelectorAll('.views button').forEach(b => b.setAttribute('aria-pressed', 'false')); }
    last = { x: e.clientX, y: e.clientY };
  });
  const up = (e) => { pointers.delete(e.pointerId); if (mode === 'magnet') { state.dragging = false; el.classList.remove('dragging'); if (state.exp === 'tube') { tube.phase = 'fall'; tube.tIn = tube.tOut = null; } else { state.mode = 'park'; state.v = 0; } } if (!pointers.size) mode = null; };
  el.addEventListener('pointerup', up); el.addEventListener('pointercancel', up);
  el.addEventListener('wheel', (e) => { e.preventDefault(); cam.r = Math.min(80, Math.max(6, cam.r * Math.exp(e.deltaY * 0.0012))); }, { passive: false });
}

/* =====================================================================
   SCOPE, METER, CAPTION
   ===================================================================== */
const scope = { buf: [], win: 5, eR: 0.01, iR: 0.001, hold: 0, u1: 'V', u2: 'A' };
function nice(x) { const e = Math.pow(10, Math.floor(Math.log10(x))), f = x / e; return (f <= 1 ? 1 : f <= 2 ? 2 : f <= 5 ? 5 : 10) * e; }
function fmt(x, unit, dig = 3) {
  if (!isFinite(x)) return '–'; if (Math.abs(x) < 1e-11) return '0 ' + unit;
  const a = Math.abs(x), pre = [[1e-9, 'n'], [1e-6, 'µ'], [1e-3, 'm'], [1, ''], [1e3, 'k']];
  let p = pre[0]; for (const q of pre) if (a >= q[0] * 0.9995) p = q;
  const v = x / p[0], d = Math.max(0, dig - 1 - Math.floor(Math.log10(Math.abs(v) || 1)));
  return v.toFixed(Math.min(d, 3)) + ' ' + p[1] + unit;
}
function cssv(n) { return getComputedStyle(document.documentElement).getPropertyValue(n).trim(); }
function drawScope() {
  const c = $('scope'), r = c.getBoundingClientRect(), dpr = Math.min(2, window.devicePixelRatio || 1);
  if (c.width !== Math.round(r.width * dpr) || c.height !== Math.round(r.height * dpr)) { c.width = Math.round(r.width * dpr); c.height = Math.round(r.height * dpr); }
  const g = c.getContext('2d'); g.setTransform(dpr, 0, 0, dpr, 0, 0);
  const W = r.width, H = r.height, soft = cssv('--ink-soft'), field = cssv('--ch1'), copper = cssv('--ch2');
  g.clearRect(0, 0, W, H);
  const L = 54, R = 56, T = 10, B = 22, pw = W - L - R, ph = H - T - B, now = state.t, data = scope.buf;
  let me = 0, mi = 0; for (const p of data) { me = Math.max(me, Math.abs(p.e)); mi = Math.max(mi, Math.abs(p.i)); }
  const wE = nice(Math.max(me * 1.15, 1e-3)), wI = nice(Math.max(mi * 1.15, 1e-4));
  if (wE > scope.eR || (now - scope.hold > 2 && wE < scope.eR)) { scope.eR = wE; scope.hold = now; }
  if (wI > scope.iR || (now - scope.hold > 2 && wI < scope.iR)) scope.iR = wI;
  // phosphor screen with divisions
  g.fillStyle = cssv('--screen'); g.beginPath(); g.moveTo(L, T); g.lineTo(L + pw - 8, T); g.lineTo(L + pw, T + 8); g.lineTo(L + pw, T + ph); g.lineTo(L + 8, T + ph); g.lineTo(L, T + ph - 8); g.closePath(); g.fill();
  g.strokeStyle = 'rgba(200,196,186,0.09)'; g.lineWidth = 1;
  for (let k = 1; k < 4; k++) { const y = T + ph * k / 4; g.beginPath(); g.moveTo(L, y); g.lineTo(L + pw, y); g.stroke(); }
  for (let k = 1; k < 10; k++) { const x = L + pw * k / 10; g.beginPath(); g.moveTo(x, T); g.lineTo(x, T + ph); g.stroke(); }
  g.strokeStyle = 'rgba(200,196,186,0.24)'; g.beginPath(); g.moveTo(L, T + ph / 2); g.lineTo(L + pw, T + ph / 2); g.stroke();
  for (let k = 0; k <= 50; k++) { const x = L + pw * k / 50; g.beginPath(); g.moveTo(x, T + ph / 2 - 2.5); g.lineTo(x, T + ph / 2 + 2.5); g.stroke(); }
  g.font = '11px Archivo, Arial, sans-serif'; g.textAlign = 'right'; g.fillStyle = field;
  g.fillText(fmt(scope.eR, scope.u1, 1), L - 6, T + 9); g.fillText(fmt(-scope.eR, scope.u1, 1), L - 6, T + ph);
  g.textAlign = 'left'; g.fillStyle = copper; g.fillText(fmt(scope.iR, scope.u2, 1), L + pw + 6, T + 9); g.fillText(fmt(-scope.iR, scope.u2, 1), L + pw + 6, T + ph);
  g.save(); g.beginPath(); g.rect(L, T, pw, ph); g.clip();
  const plot = (key, range, color, w) => { g.strokeStyle = color; g.lineWidth = w; g.shadowColor = color; g.shadowBlur = 6; g.beginPath(); let st = false;
    for (const p of data) { const x = L + pw * (p.t - (now - scope.win)) / scope.win, y = T + ph / 2 - Math.max(-1.05, Math.min(1.05, p[key] / range)) * ph / 2; if (!st) { g.moveTo(x, y); st = true; } else g.lineTo(x, y); } g.stroke(); };
  plot('i', scope.iR, copper, 1.4); plot('e', scope.eR, field, 1.8); g.restore();
  g.fillStyle = soft; g.textAlign = 'center'; g.font = '11px Archivo, Arial, sans-serif'; g.fillText('0.5 s / div', L + pw / 2, H - 6);
}
function drawRings(g, W, H) {
  const soft = cssv('--ink-soft'), faint = cssv('--tick-soft'), ccw = cssv('--ccw'), cw = cssv('--cw'), ink = cssv('--ink');
  const L = 50, R = 12, T = 12, B = 34, pw = W - L - R, ph = H - T - B, X = (y) => L + (y + 0.2) / 0.45 * pw;
  g.font = '10.5px Archivo, Arial, sans-serif'; g.fillStyle = soft; g.strokeStyle = faint; g.lineWidth = 1;
  for (let y = -0.2; y <= 0.2501; y += 0.05) { g.beginPath(); g.moveTo(X(y), T); g.lineTo(X(y), T + ph); g.stroke(); g.textAlign = 'center'; g.fillText(Math.round(y * 100), X(y), T + ph + 13); }
  g.fillText('height, cm (tube runs from −15 to +15)', L + pw / 2, T + ph + 27);
  g.fillStyle = 'rgba(128,128,120,0.12)'; g.fillRect(X(-0.15), T, X(0.15) - X(-0.15), ph);
  const imax = Math.max(...Array.from(tube.I, Math.abs)), yr = nice(Math.max(imax * 1.15, 0.05)), Y = (i) => T + ph / 2 - i / yr * ph / 2;
  g.strokeStyle = soft; g.globalAlpha = 0.5; g.beginPath(); g.moveTo(L, Y(0)); g.lineTo(L + pw, Y(0)); g.stroke(); g.globalAlpha = 1;
  g.fillStyle = soft; g.textAlign = 'right'; g.fillText(fmt(yr, 'A', 1), L - 5, T + 8); g.fillText(fmt(-yr, 'A', 1), L - 5, T + ph);
  const bw = Math.max(1, (X(TUBE.hr) - X(0)) * 0.8);
  for (let k = 0; k < TUBE.N; k++) { const i = tube.I[k], x = X(tube.z[k]); g.fillStyle = i >= 0 ? ccw : cw; g.fillRect(x - bw / 2, Math.min(Y(0), Y(i)), bw, Math.abs(Y(i) - Y(0))); }
  const xm = X(tube.y); g.fillStyle = 'rgba(217,139,76,0.18)'; g.fillRect(X(tube.y - MAGNET.b), T, X(tube.y + MAGNET.b) - X(tube.y - MAGNET.b), ph);
  g.strokeStyle = ink; g.setLineDash([2, 3]); g.beginPath(); g.moveTo(xm, T); g.lineTo(xm, T + ph); g.stroke(); g.setLineDash([]);
  g.textAlign = 'left'; g.fillStyle = ccw; g.fillText('counterclockwise seen from above', L + 8, T + 11); g.fillStyle = cw; g.fillText('clockwise', L + 8, T + 24);
}
function drawLam() {
  const c = $('lam'), r = c.getBoundingClientRect(), dpr = Math.min(2, window.devicePixelRatio || 1);
  if (!r.width) return;
  if (c.width !== Math.round(r.width * dpr) || c.height !== Math.round(r.height * dpr)) { c.width = Math.round(r.width * dpr); c.height = Math.round(r.height * dpr); }
  const g = c.getContext('2d'); g.setTransform(dpr, 0, 0, dpr, 0, 0);
  const W = r.width, H = r.height, soft = cssv('--ink-soft'), faint = cssv('--tick-soft'), c1 = cssv('--ch1'), ink = cssv('--ink');
  g.clearRect(0, 0, W, H);
  if (state.exp === 'tube') { drawRings(g, W, H); return; }
  const L = 54, R = 12, T = 12, B = 34, pw = W - L - R, ph = H - T - B;
  g.font = '10.5px Archivo, Arial, sans-serif'; g.fillStyle = soft; g.strokeStyle = faint; g.lineWidth = 1;
  const X = (x) => L + (x + 0.15) / 0.3 * pw;
  for (let x = -0.15; x <= 0.1501; x += 0.05) { g.beginPath(); g.moveTo(X(x), T); g.lineTo(X(x), T + ph); g.stroke(); g.textAlign = 'center'; g.fillText(Math.round(x * 100), X(x), T + ph + 13); }
  g.fillText('magnet position x, cm', L + pw / 2, T + ph + 27);
  if (TAB.n < TAB.N || TAB.key !== tabKey()) { g.textAlign = 'center'; g.fillText('computing Λ(x)', L + pw / 2, T + ph / 2); return; }
  const withCore = state.core && TABM.n === TABM.N;   // total flux linkage with the iron magnetized by the magnet alone
  const mk = TAB.key + '|' + state.core; if (vis.maskKey !== mk) { vis.mask = new Uint8Array(TAB.N); for (let i = 0; i < TAB.N; i++) vis.mask[i] = collides(TAB.x0 + i * TAB.h, TAB.psi, TAB.th) ? 0 : 1; vis.maskKey = mk; }
  const lamAt = (i) => TAB.vals[i] + (withCore ? core.K * coreJ(TABM.vals[i], 0) : 0);
  let mx = 0; for (let i = 0; i < TAB.N; i++) if (vis.mask[i]) mx = Math.max(mx, Math.abs(lamAt(i)));
  const yr = nice(Math.max(mx * 1.1, 1e-6)), Y = (v) => T + ph / 2 - v / yr * ph / 2;
  g.beginPath(); g.moveTo(L, Y(0)); g.lineTo(L + pw, Y(0)); g.strokeStyle = soft; g.globalAlpha = 0.5; g.stroke(); g.globalAlpha = 1;
  g.textAlign = 'right'; g.fillText(fmt(yr, 'Wb', 1), L - 5, T + 8); g.fillText(fmt(-yr, 'Wb', 1), L - 5, T + ph);
  g.beginPath(); let pen = false; for (let i = 0; i < TAB.N; i++) { if (!vis.mask[i]) { pen = false; continue; } const x = X(TAB.x0 + i * TAB.h), y = Y(lamAt(i)); pen ? g.lineTo(x, y) : g.moveTo(x, y); pen = true; }
  for (let i = 0; i < TAB.N; i++) if (!vis.mask[i]) { g.fillStyle = 'rgba(128,128,120,0.14)'; g.fillRect(X(TAB.x0 + i * TAB.h) - 1, T, pw / 300 + 2, ph); }
  g.strokeStyle = cssv('--ch3') || ink; g.lineWidth = 1.6; g.stroke();
  // tangent at the magnet: its slope times the speed gives the EMF from motion
  const e = tabEval(state.s); if (!e) return;
  const em = withCore ? tabEval(state.s, TABM) : null, Gs = state.Gs, lam0 = e[0] + (em ? core.K * coreJ(em[0], 0) : 0), x0 = X(state.s), y0 = Y(lam0), dx = 0.03; e[0] = lam0;
  g.beginPath(); g.moveTo(X(state.s - dx), Y(e[0] - Gs * dx)); g.lineTo(X(state.s + dx), Y(e[0] + Gs * dx)); g.strokeStyle = c1; g.lineWidth = 2; g.stroke();
  g.beginPath(); g.moveTo(x0, T); g.lineTo(x0, T + ph); g.strokeStyle = soft; g.setLineDash([2, 3]); g.lineWidth = 1; g.stroke(); g.setLineDash([]);
  g.fillStyle = ink; g.beginPath(); g.arc(x0, y0, 3.5, 0, 2 * PI); g.fill();
  const emfM = -Gs * state.v;
  g.textAlign = 'left'; g.fillStyle = c1; const tx = L + 8;
  g.fillText(`slope ${fmt(Gs, 'Wb/m', 2)}`, tx, T + 11); g.fillStyle = soft;
  g.fillText(`EMF = −slope × ${state.v.toFixed(2)} m/s = ${fmt(emfM, 'V', 2)}`, tx, T + 24);
}
const meter = { range: 0.01 };
function initMeter() {
  const cx = 105, cy = 140, R = 88, a0 = -60 * PI / 180, a1 = 60 * PI / 180, pt = (a, r) => [cx + r * Math.sin(a), cy - r * Math.cos(a)];
  const [x0, y0] = pt(a0, R), [x1, y1] = pt(a1, R);
  $('arc').setAttribute('d', `M ${x0} ${y0} A ${R} ${R} 0 0 1 ${x1} ${y1}`);
  let t = ''; for (let i = -5; i <= 5; i++) { const a = i / 5 * a1, [p, q] = pt(a, R), [u, v] = pt(a, R - (i % 5 === 0 ? 12 : 6)); t += `<line x1="${p}" y1="${q}" x2="${u}" y2="${v}"/>`; }
  for (let i = -50; i <= 50; i += 1) { if (i % 10 === 0) continue; const a = i / 50 * a1, [p, q] = pt(a, R), [u, v] = pt(a, R - 3); t += `<line x1="${p}" y1="${q}" x2="${u}" y2="${v}" opacity="0.5"/>`; }
  $('ticks').innerHTML = t;
}
function drawMeter() {
  const want = nice(Math.max(Math.abs(state.exp === 'tube' ? (tube.iBelow || 0) : state.I) * 1.2, 1e-3));
  if (want > meter.range) meter.range = want; else if (want < meter.range / 4) meter.range = want;
  const ink = cssv('--ink'), soft = cssv('--ink-soft'), cop = cssv('--copper');
  const val = state.exp === 'tube' ? (tube.iBelow || 0) : state.I, f = Math.max(-1.08, Math.min(1.08, val / meter.range)), a = f * 60 * PI / 180;
  const n = $('needle'); n.setAttribute('x2', 105 + 84 * Math.sin(a)); n.setAttribute('y2', 140 - 84 * Math.cos(a)); n.setAttribute('stroke', cop);
  $('arc').setAttribute('stroke', soft); $('ticks').setAttribute('stroke', soft);
  document.querySelectorAll('.meter text').forEach(t => t.setAttribute('fill', soft)); document.querySelector('.meter circle').setAttribute('fill', ink);
  $('tmin').textContent = fmt(-meter.range, 'A', 1); $('tmax').textContent = fmt(meter.range, 'A', 1);
}
function valUnit(x, unit, dig) { const t = fmt(x, unit, dig), k = t.lastIndexOf(' '); return `${t.slice(0, k)}<small>${t.slice(k + 1)}</small>`; }
function updateTubeReadouts() {
  const mg = magnetMass() * 9.81, F = tube.F || 0, ff = Math.sqrt(2 * 0.30 / 9.81);
  $('kEmf').textContent = 'Speed'; $('vEmf').innerHTML = valUnit(Math.abs(tube.v), 'm/s', 3);
  $('xEmf').textContent = tube.cond ? `terminal about ${fmt(tube.vTerm, 'm/s', 2)}` : 'no drag: free fall';
  $('kCur').textContent = 'Magnetic drag'; $('vCur').innerHTML = valUnit(F, 'N', 3); $('xCur').textContent = `${Math.round(F / mg * 100)}% of the weight, ${fmt(mg, 'N', 3)}`;
  $('hot').classList.remove('on');
  $('kB').textContent = 'Largest ring current'; $('vB').innerHTML = valUnit(Math.abs(tube.imax || 0), 'A', 3); $('xB').textContent = tube.cond ? `ring resistance ${fmt(tube.Rk, 'Ω', 3)}` : 'plastic does not conduct';
  const tr = tube.tIn !== null && tube.tOut !== null ? tube.tOut - tube.tIn : null;
  $('kF').textContent = 'Heat in the tube'; $('vF').innerHTML = valUnit(tube.heat || 0, 'W', 3); $('xF').textContent = tr !== null ? `${tr.toFixed(2)} s through the tube` : 'transit time pending';
  const T = TUBES[tube.mat], spec = [['Tube', `${T.nm}, ID 14 mm, wall ${(tube.t * 1000).toFixed(1)} mm`], ['Resistivity', isFinite(T.rho) ? `${(T.rho * 1e9).toPrecision(3)} nΩ·m` : 'insulator'],
    ['Length', '30 cm, 120 rings of 2.5 mm'], ['Ring L/R', tube.cond ? fmt(tube.Ls / tube.Rk, 's', 2) : 'none'], ['Magnet', `Ø 12 × 25 mm, ${(magnetMass() * 1000).toFixed(1)} g`], ['Free fall, 30 cm', `${ff.toFixed(2)} s from rest`]];
  const sh = spec.map(r => `<dt>${r[0]}</dt><dd>${r[1]}</dd>`).join(''); if (sh !== vis.specHtml) { $('spec').innerHTML = sh; vis.specHtml = sh; }
}
function updateReadouts() {
  if (state.exp === 'tube') { updateTubeReadouts(); return; }
  $('kEmf').textContent = 'Induced EMF'; $('kCur').textContent = 'Current'; $('kB').textContent = state.core ? 'Coil and core field, center' : 'Coil field, center'; $('kF').textContent = 'Force on magnet';
  const Rt = coil.R + R_METER, Ib = state.V / Rt, Ii = state.I - Ib, P = state.I * state.I * coil.R;
  const F = state.force, toward = F * -Math.sign(state.s || -1) > 0;
  $('vEmf').innerHTML = valUnit(state.emf, 'V', 3); $('xEmf').textContent = 'flux ' + fmt(state.lam, 'Wb', 3);
  $('vCur').innerHTML = valUnit(state.I, 'A', 3); $('xCur').textContent = state.V !== 0 ? `battery ${fmt(Ib, 'A', 2)}, induced ${fmt(Ii, 'A', 1)}` : `heat ${fmt(P, 'W', 2)}`;
  $('hot').classList.toggle('on', P > 5); $('hot').textContent = fmt(P, 'W', 2) + ' hot';
  $('vB').innerHTML = valUnit(vis.bPerAmp * state.I + (state.core ? magnetField(1e-9, 0, IRON.rc, IRON.h, state.J)[1] : 0), 'T', 3); $('xB').textContent = `${state.N} turns, ${fmt(coil.R, 'Ω', 3)}`;
  const wireLen = state.N * 2 * PI * coil.rMean, spec = [
    ['Magnet', `Ø 12 × 25 mm, ${(magnetMass() * 1000).toFixed(1)} g`], ['Remanence', `${MAGNET.Br.toFixed(2)} T`],
    ['Coil', `bore 30 mm, ${state.N} turns, ${coil.layers} layers`], ['Wire', `${WIRES[state.wire].sym}, Ø 0.5 mm, ${wireLen.toFixed(1)} m`],
    ['Resistance', `${fmt(coil.R, 'Ω', 3)} + 1 Ω meter`], ['Inductance', state.core ? `${fmt(coil.L + core.Ladd, 'H', 3)} with core` : `${fmt(coil.L, 'H', 3)} (Wheeler)`], ['Time constant', fmt((coil.L + (state.core ? core.Ladd : 0)) / (coil.R + R_METER), 's', 2)]];
  if (state.core) spec.push(['Iron core', `Ø 27 × 40 mm, N ${core.Nd.toFixed(3)}, χ eff ${core.chiE.toFixed(2)}`], ['Core magnetization', `J ${fmt(state.J, 'T', 2)} of 2 T`]);
  const sh = spec.map(r => `<dt>${r[0]}</dt><dd>${r[1]}</dd>`).join(''); if (sh !== vis.specHtml) { $('spec').innerHTML = sh; vis.specHtml = sh; }
  $('vF').innerHTML = Math.abs(F) < 1e-4 ? '0<small>N</small>' : valUnit(Math.abs(F), 'N', 3); $('xF').textContent = Math.abs(F) < 1e-4 ? 'no force' : toward ? 'pull toward the coil' : 'push away from the coil';
}
// One plain sentence for what is happening right now
function updateTubeCaption() {
  const mg = magnetMass() * 9.81, F = tube.F || 0, y = tube.y; let msg, tone = '';
  if (tube.phase === 'drag') { msg = 'Move it fast and the tube pushes back. Hold it still and nothing happens.'; tone = Math.abs(F) > 1e-3 ? 'live' : ''; }
  else if (tube.phase === 'landed') { const tr = tube.tIn !== null && tube.tOut !== null ? tube.tOut - tube.tIn : null; msg = tr !== null ? `Landed. It took ${tr.toFixed(2)} s to pass the 30 cm tube.` : 'Landed. Press Drop to run it again.'; }
  else if (y - MAGNET.b > TUBE.top) msg = 'Falling freely toward the tube.';
  else if (y + MAGNET.b < -TUBE.top) msg = 'Out of the tube and falling freely again.';
  else if (!tube.cond) msg = 'Plastic does not conduct, so no eddy currents slow the fall.';
  else if (F > 0.9 * mg) { msg = `Terminal speed, ${fmt(Math.abs(tube.v), 'm/s', 2)}: eddy current drag balances the weight.`; tone = 'live'; }
  else { msg = 'Eddy currents circle the tube above and below the magnet and push back on it.'; tone = 'live'; }
  $('status').textContent = msg; $('dot').className = 'led ' + tone;
}
function updateCaption() {
  if (state.exp === 'tube') { updateTubeCaption(); return; }
  const Rt = coil.R + R_METER, Ib = state.V / Rt, e = state.emf, rising = -Math.sign(state.lam) * e > 0;
  let msg, tone = '';
  if (state.corePending !== state.core) msg = 'The iron core goes in once the magnet is clear of the coil.';
  else if (state.coilWaiting) msg = 'The coil turns once the magnet is clear of it.';
  else if (state.core && state.V === 0 && Math.abs(state.emf) <= 1e-3 && Math.abs(state.F0) > 2e-3) { msg = 'The magnet magnetizes the iron core, and the iron pulls the magnet in.'; tone = 'warm'; }
  else if (state.V !== 0 && Math.abs(state.forceBatt) > 2e-3) {
    const nFacing = (Ib > 0) === (state.s * coilAxis(state.theta)[0] > 0), toward = state.forceBatt * -Math.sign(state.s || -1) > 0;
    msg = `Electromagnet on: its ${nFacing ? 'north' : 'south'} pole faces the magnet, so it ${toward ? 'pulls it in' : 'pushes it away'}.`; tone = 'warm';
  } else if (Math.abs(e) > 1e-3) { msg = `Flux ${rising ? 'rising' : 'falling'}, so ${fmt(Math.abs(e), 'V', 2)} is induced.`; tone = 'live'; }
  else if (Math.abs(Math.cos(state.theta)) < 0.25) msg = 'Coil side on: field crosses it both ways and cancels.';
  else if (state.V !== 0) { msg = 'Electromagnet on, but too far or too turned to pull the magnet.'; tone = 'warm'; }
  else if (Math.abs(state.v) > 0.02) msg = 'Magnet far away, so the flux barely changes.';
  else msg = Math.abs(state.lam) > 1e-4 ? 'Field through the coil, but steady, so nothing is induced.' : 'Nothing changes, so nothing is induced.';
  $('status').textContent = msg; $('dot').className = 'led ' + tone;
}

/* =====================================================================
   UI and MAIN LOOP
   ===================================================================== */
const $ = (id) => document.getElementById(id);
const ui = { toast(msg) { const t = $('toast'); t.textContent = msg; t.classList.add('show'); clearTimeout(ui._t); ui._t = setTimeout(() => t.classList.remove('show'), 2600); } };
function setTurns(N) {
  state.N = N; $('turnsOut').textContent = N; buildCoil();
  if (collides(state.s, state.psi, state.theta)) { state.theta = 0; }
  freeCache.clear(); rebuildCoilVisuals(); computeBPerAmp(); FL.doneKey = '';
}
function computeBPerAmp() { let b = 0; for (const L of coil.loops) b += loopB(1e-9, -L.z, L.r, state.N * L.w)[1]; vis.bPerAmp = b; }
function setCoreNow(on) {
  state.core = on; TABM.n = 0; FL.doneKey = ''; FL.mapKey = ''; freeCache.clear();
  vis.coilPoleX = Math.max(BOBBIN.halfLen * S + BOBBIN.flange * S, on ? IRON.h * S : 0) + 0.75;
}
function setCore(on) {
  state.corePending = on; freeCache.clear();
  $('coreAir').setAttribute('aria-pressed', !on); $('coreIron').setAttribute('aria-pressed', on);
  if (!on) setCoreNow(false);
}
function setExp(e) {
  state.exp = e; scope.buf.length = 0; scope.eR = 0.01; scope.iR = 0.001;
  document.querySelectorAll('[data-exp]').forEach(el => el.hidden = el.dataset.exp !== e);
  $('expCoil').setAttribute('aria-pressed', e === 'coil'); $('expTube').setAttribute('aria-pressed', e === 'tube');
  if (e === 'tube') {
    scope.u1 = 'm/s'; scope.u2 = 'N'; $('capScope').textContent = 'CH1 falling speed, CH2 magnetic drag'; $('capLam').innerHTML = 'Fig. 3&nbsp; Eddy current in each ring'; $('capLamR').textContent = 'opposite signs above and below';
    $('vcapTxt').textContent = 'Magnet falling through a conducting tube, cut away to show the eddy currents'; $('mlabel').textContent = 'Eddy current, ring below magnet';
    if (tube.dirty) buildTube(); rebuildTubeVisuals(); dropMagnet(); setView('oblique');
    state.show.map = false; $('showMap').setAttribute('aria-pressed', false);   // the map plane would hide the eddy bands
  } else {
    scope.u1 = 'V'; scope.u2 = 'A'; $('capScope').textContent = 'CH1 induced EMF, CH2 induced current'; $('capLam').innerHTML = 'Fig. 3&nbsp; Flux linkage Λ<span class="long"> against magnet position</span>'; $('capLamR').textContent = 'EMF = −slope × speed';
    $('vcapTxt').textContent = 'Field of the magnet and coil, horizontal plane through the axis'; $('mlabel').textContent = 'Galvanometer, coil current';
    state.dragging = false; setView('oblique');
  }
  vis.specHtml = '';
}
function setTubeMat(k) {
  tube.mat = k; tube.dirty = true; document.querySelectorAll('#tubes .cell').forEach(t => t.setAttribute('aria-pressed', t.dataset.k === k));
  $('tubeName').textContent = TUBES[k].nm + ' tube'; buildTube(); rebuildTubeVisuals(); $('tubeRho').textContent = isFinite(TUBES[k].rho) ? fmt(tube.Rk, 'Ω', 2) + ' per ring' : 'insulator';
}
function setMaterial(k) {
  state.mat = k; MAGNET.Br = MAGNETS[k].Br; MAGNET.rho = MAGNETS[k].rho;
  document.querySelectorAll('#mats .cell').forEach(t => t.setAttribute('aria-pressed', t.dataset.k === k));
  $('matName').textContent = MAGNETS[k].name; $('matBr').textContent = MAGNETS[k].Br.toFixed(2) + ' T';
  buildRigidLines(); buildRigidMap(); vis.glowRef = Math.abs(innerTurnFlux(0, 0, 0)); FL.doneKey = ''; FL.mapKey = ''; tube.dirty = true; if (state.exp === 'tube') buildTube();
}
function setWire(k) {
  state.wire = k; buildCoil(); rebuildCoilVisuals(); computeBPerAmp(); FL.doneKey = ''; FL.mapKey = '';
  document.querySelectorAll('#wires .cell').forEach(t => t.setAttribute('aria-pressed', t.dataset.k === k));
  $('wireName').textContent = WIRES[k].nm + ' wire, 0.5 mm'; $('wireR').textContent = fmt(coil.R, 'Ω', 3);
}
function initUI() {
  $('mats').innerHTML = Object.entries(MAGNETS).map(([k, m]) => `<button class="cell" data-k="${k}" title="${m.name}, remanence ${m.Br} T"><span class="z">${m.Br.toFixed(2)} T</span><span class="sym">${m.sym}</span><span class="nm">${m.nm}</span></button>`).join('');
  $('wires').innerHTML = Object.entries(WIRES).map(([k, w]) => `<button class="cell" data-k="${k}" title="${w.nm}, resistivity ${(w.rho * 1e9).toPrecision(3)} nΩ·m"><span class="z">${w.rho * 1e9 >= 100 ? Math.round(w.rho * 1e9) : (w.rho * 1e9).toFixed(1)} nΩm</span><span class="sym">${w.sym}</span><span class="nm">${w.nm}</span></button>`).join('');
  document.querySelectorAll('#mats .cell').forEach(t => t.addEventListener('click', () => setMaterial(t.dataset.k)));
  document.querySelectorAll('#wires .cell').forEach(t => t.addEventListener('click', () => setWire(t.dataset.k)));
  // ruler sliders: fill from the minimum, or from zero for the battery
  const fillRule = (el) => { const lo = +el.min, hi = +el.max, f = (v) => ((v - lo) / (hi - lo) * 100).toFixed(2) + '%', v = +el.value;
    if (el.classList.contains('bipolar')) { el.style.setProperty('--a', f(Math.min(0, v))); el.style.setProperty('--b', f(Math.max(0, v))); }
    else { el.style.setProperty('--a', '0%'); el.style.setProperty('--b', f(v)); } };
  document.querySelectorAll('input.rule').forEach(el => { el.addEventListener('input', () => fillRule(el)); fillRule(el); });
  $('swing').addEventListener('click', () => setSwing(state.mode !== 'swing'));
  $('flip').addEventListener('click', flipMagnet);
  $('speed').addEventListener('input', e => { state.swing.speed = +e.target.value; $('speedOut').textContent = (+e.target.value).toFixed(2) + ' m/s'; });
  $('volt').addEventListener('input', e => { state.V = +e.target.value; $('voltOut').textContent = state.V === 0 ? 'off' : (state.V > 0 ? '+' : '') + state.V.toFixed(1) + ' V'; });
  $('angle').addEventListener('input', e => { state.thetaTarget = +e.target.value * PI / 180; $('angleOut').textContent = e.target.value + '°'; });
  $('turns').addEventListener('input', e => { setTurns(+e.target.value); $('wireR').textContent = fmt(coil.R, 'Ω', 3); });
  for (const [id, k] of [['showLines', 'lines'], ['showMap', 'map']]) $(id).addEventListener('click', e => { state.show[k] = !state.show[k]; e.currentTarget.setAttribute('aria-pressed', state.show[k]); FL.mapKey = ''; });
  document.querySelectorAll('.views button').forEach(b => b.addEventListener('click', () => setView(b.dataset.view)));
  $('speedOut').textContent = '0.30 m/s'; $('voltOut').textContent = 'off'; $('angleOut').textContent = '0°'; $('turnsOut').textContent = state.N;
  // color bar ticks at decades of |B|
  $('bar').innerHTML = [[1, '1 T'], [0.1, '100 mT'], [0.01, '10 mT'], [0.001, '1 mT']].map(([b, l]) => `<span class="tk" style="left:${(bT(b) * 100).toFixed(1)}%">${l}</span>`).join('');
  $('tubes').innerHTML = Object.entries(TUBES).map(([k, w]) => `<button class="cell" data-k="${k}" title="${w.nm}"><span class="z">${isFinite(w.rho) ? (w.rho * 1e9 >= 100 ? Math.round(w.rho * 1e9) : (w.rho * 1e9).toFixed(1)) + ' nΩm' : 'insulator'}</span><span class="sym">${w.sym}</span><span class="nm">${w.nm}</span></button>`).join('');
  document.querySelectorAll('#tubes .cell').forEach(t => t.addEventListener('click', () => { setTubeMat(t.dataset.k); dropMagnet(); }));
  $('wall').addEventListener('input', e => { tube.t = +e.target.value / 1000; $('wallOut').textContent = (+e.target.value).toFixed(1) + ' mm'; tube.dirty = true; });
  $('wall').addEventListener('change', () => { buildTube(); rebuildTubeVisuals(); setTubeMat(tube.mat); dropMagnet(); });
  $('wallOut').textContent = '1.5 mm';
  $('drop').addEventListener('click', dropMagnet);
  $('tflip').addEventListener('click', () => { tube.up = !tube.up; tube.dirty = true; dropMagnet(); });
  $('expCoil').addEventListener('click', () => setExp('coil')); $('expTube').addEventListener('click', () => setExp('tube'));
  $('coreAir').addEventListener('click', () => setCore(false)); $('coreIron').addEventListener('click', () => setCore(true));
  setMaterial('n42'); setWire('cu'); setTubeMat('cu');
}
function resize() {
  const v = $('view3d'), w = v.clientWidth, h = v.clientHeight; if (!w || !h) return;
  renderer.setSize(w, h, false); camera.aspect = w / h; camera.updateProjectionMatrix();
}
let lastT = performance.now(), roT = 1;
function frame(now) {
  if (window.__paused) return;
  const dt = Math.min(1 / 30, Math.max(1e-4, (now - lastT) / 1000)); lastT = now;
  if (state.exp === 'tube') simulateTube(dt); else { fillTable(4); simulate(dt); }
  if (cam.anim) { const A = cam.anim; A.t += dt; const x = Math.min(1, A.t / A.dur), e = x * x * (3 - 2 * x); cam.az = A.from.az + (A.to.az - A.from.az) * e; cam.el = A.from.el + (A.to.el - A.from.el) * e; cam.r = A.from.r + (A.to.r - A.from.r) * e; if (x >= 1) cam.anim = null; }
  placeCamera(); updateVisuals(dt);
  for (const q of state.samples) scope.buf.push({ t: q[0], e: q[1], i: q[2] });
  while (scope.buf.length && scope.buf[0].t < state.t - scope.win) scope.buf.shift();
  drawScope(); drawMeter(); drawLam();
  $('flip').setAttribute('aria-pressed', !!state.flip); $('clock').textContent = `t ${state.t.toFixed(1)} s   step 2.1 ms`;
  if (state.exp === 'tube') $('drop').setAttribute('aria-pressed', tube.phase === 'fall');
  roT += dt; if (roT > 0.2) { roT = 0; updateReadouts(); updateCaption(); }
  renderer.render(scene, camera);
  requestAnimationFrame(frame);
}
function main() {
  buildCoil(); initScene(); computeBPerAmp(); initPointer(); initUI(); initMeter();
  fillTable(2000); setSwing(true);
  new ResizeObserver(resize).observe($('view3d')); resize();
  $('loading').remove(); updateReadouts(); updateCaption();
  window.__bench = { state, tube, TUBE, core, IRON, simulateTube, dropMagnet, buildTube, setExp, setCore, setTubeMat, setMaterial, TAB, TABM, simulate, fillTable, updateVisuals, updateFieldLines, updateCaption, updateReadouts, setTurns, setSwing, flipMagnet, collides, lambda, FL, vis, cam, scope, coil: () => coil,
    renderOnce: () => { placeCamera(); updateVisuals(1 / 60); drawScope(); drawMeter(); drawLam(); updateReadouts(); updateCaption(); renderer.render(scene, camera); }, pause: () => { window.__paused = true; } };
  requestAnimationFrame(frame);
}
requestAnimationFrame(main);

/* Notes on method
   Magnet field: exact solution for a uniformly magnetized cylinder (Derby and Olbert 2010).
   Flux linkage: by reciprocity, a line integral of the coil's vector potential around the magnet's cross
   sections, valid at any tilt or offset. Coil: N turns over the winding section, Gauss Legendre sampling.
   Resistance from copper resistivity; inductance from Wheeler's multilayer formula.
   Combined field lines: magnet field plus the coil's field from its present current, traced with RK4. */
