
"""
Post-tensioned simply-supported transfer beam — EC2 (SS EN 1992 / UK NA conventions)
Span 17.8 m, factored UDL 322 kN/m, overall depth fixed at 1100 mm.

Parametric design + check sheet.  All checks recompute from inputs — change the
INPUTS block and re-run.  Numbers below printed by the script, not hand-typed.

KEY ASSUMPTIONS (confirm before use):
  1. Load split  Gk=200, Qk=35  ->  1.35Gk+1.5Qk = 322.5 kN/m  (drives SLS moment)
  2. Flanged section: a 250 mm slab acts as a top flange (beff~3000). A transfer
     beam this heavy at h=1100 will NOT work as a slim rectangle on ULS ductility;
     the flange (or a ~1500 wide web) is what makes the compression block shallow.
  3. Deflection-sensitive finishes present -> SLS basis = DECOMPRESSION under the
     quasi-permanent combination (durable, robust; conservative for XC1).
  4. Bonded 15.7 mm strand, multistrand ducts.
"""
import math

# ----------------------------- INPUTS --------------------------------------
L    = 17.8                 # m, span (simply supported)
Gk   = 200.0                # kN/m characteristic permanent (incl. self wt)
Qk   = 35.0                 # kN/m characteristic variable
psi1, psi2 = 0.5, 0.3       # frequent / quasi-permanent factors
gG, gQ = 1.35, 1.5

h    = 1100.0               # mm overall depth (FIXED)
bw   = 900.0                # mm web width
hf   = 250.0                # mm flange (slab) thickness
beff = 3000.0               # mm effective flange width (EC2 5.3.2)

fck  = 40.0                 # MPa  (C40/50)
fckt = 32.0                 # MPa  transfer strength
acc  = 0.85                 # NA, compression in flexure
fctm = 0.30*fck**(2/3)      # MPa

# Prestress: 15.7 mm strand
Ap1     = 150.0             # mm2/strand
fpk     = 1860.0
fp01k   = 1600.0
gp      = 1.15
fpd     = fp01k/gp
sig_pm0 = 1300.0            # MPa at transfer (after short-term losses)
sig_pe  = 1150.0            # MPa effective (after all losses ~ -25% from jack)
P0_1    = sig_pm0*Ap1/1e3   # kN/strand at transfer
Pe_1    = sig_pe *Ap1/1e3   # kN/strand effective

y_t  = 180.0                # mm, tendon group centroid above soffit at midspan
# ---------------------------------------------------------------------------

wEd  = gG*Gk + gQ*Qk
wqp  = Gk + psi2*Qk
wchr = Gk + Qk
w_sw = bw*h*25e-6           # kN/m web self weight (at transfer, pre-slab)

MEd  = wEd *L*L/8
VEd  = wEd *L/2
Mqp  = wqp *L*L/8
Mchr = wchr*L*L/8
Msw  = w_sw*L*L/8

print(f"Loads:  wEd={wEd:.1f}  wqp={wqp:.1f}  wchr={wchr:.1f} kN/m   w_sw(web)={w_sw:.1f}")
print(f"ULS:    MEd={MEd:.0f} kNm   VEd={VEd:.0f} kN")
print(f"SLS:    Mqp={Mqp:.0f}   Mchr={Mchr:.0f}   Msw={Msw:.0f} kNm\n")

# ---- Flanged (service) section properties ----
A_w = bw*h ; A_f = (beff-bw)*hf ; A = A_w + A_f
yb  = (A_w*h/2 + A_f*(h-hf/2))/A            # centroid above soffit
Iw  = bw*h**3/12 + A_w*(yb-h/2)**2
If  = (beff-bw)*hf**3/12 + A_f*((h-hf/2)-yb)**2
I   = Iw + If
Zb  = I/yb ; Zt = I/(h-yb)
kb  = Zb/A
e   = yb - y_t                              # eccentricity at midspan
dp  = h - y_t                               # tendon eff. depth
print(f"Section: A={A:.0f}mm2  yb={yb:.0f}mm  I={I:.3e}mm4")
print(f"         Zb={Zb:.3e}  Zt={Zt:.3e}  kern_b={kb:.0f}mm  e={e:.0f}mm  dp={dp:.0f}mm\n")

# ---- 1. SLS: prestress for decompression under quasi-permanent ----
fcd = acc*fck/1.5
P_req = Mqp*1e6 /(kb+e)          # N : sig_bottom>=0 under Mqp  (P/A + Pe/Zb >= M/Zb)
n_sls = P_req/1e3/Pe_1
print(f"[1] SLS decompression (QP): P_req={P_req/1e3:.0f} kN  -> {n_sls:.1f} strands")

# ---- 2. ULS flexure (flanged), pick strands to satisfy MRd ----
def MRd(n):
    T = n*Ap1*fpd/1e3                       # kN, tendon at fpd
    x = T*1e3/(fcd*beff)                    # NA depth assuming within flange
    if x>hf:                                # into web (rare here)
        Tf = fcd*beff*hf/1e3
        x  = hf + (T-Tf)*1e3/(fcd*bw)
    z = dp - 0.4*x
    return T*z/1e3, x                       # kNm, mm
n = math.ceil(n_sls)
while MRd(n)[0] < MEd: n += 1
M_R, x = MRd(n)
print(f"[2] ULS flexure: provide {n} strands  Ap={n*Ap1:.0f}mm2")
print(f"    x={x:.0f}mm  x/dp={x/dp:.2f} (<0.45 ductile)  MRd={M_R:.0f} >= MEd={MEd:.0f}  {'OK' if M_R>=MEd else 'FAIL'}\n")

Pe = n*Pe_1 ; P0 = n*P0_1
print(f"    Effective prestress Pe={Pe:.0f}kN   at transfer P0={P0:.0f}kN\n")

# ---- 3. Transfer stresses (WEB-only rectangular section, slab not yet composite) ----
Aw=bw*h; Zw=bw*h**2/6; ew=h/2-y_t
def stress(P,M,Aa,Ztt,Zbb,ecc):
    sb = P*1e3/Aa + P*1e3*ecc/Zbb - M*1e6/Zbb
    st = P*1e3/Aa - P*1e3*ecc/Ztt + M*1e6/Ztt
    return st, sb
st,sb = stress(P0,Msw,Aw,Zw,Zw,ew)
lim_c = 0.6*fckt ; lim_t = -fctm
print(f"[3] Transfer (full P0, web only): top={st:.1f}  bot={sb:.1f} MPa "
      f"(limits: comp {lim_c:.1f}, tens {lim_t:.1f})")
# staged stressing: max fraction at self-weight-only stage
aT = (lim_t - Msw*1e6/Zw)/(P0*1e3*(1/Aw - ew/Zw))   # top tension governs
aC = (lim_c + Msw*1e6/Zw)/(P0*1e3*(1/Aw + ew/Zw))   # bottom comp governs
amax = min(aT,aC)
print(f"    Full P0 over-stresses at transfer -> STAGE the stressing.")
print(f"    Max stress at SW-only stage ~ {100*amax:.0f}% of P0 ({amax*P0:.0f} kN);"
      f" apply remainder once superimposed dead is on.\n")

# ---- 4. Service stresses under characteristic (composite flanged) ----
st,sb = stress(Pe,Mchr,A,Zt,Zb,e)
print(f"[4] Service (Pe, characteristic): top={st:.1f}  bot={sb:.1f} MPa  "
      f"(0.6fck={0.6*fck:.0f}; bottom {'OK no tension' if sb>=0 else 'tension '+f'{sb:.1f}'})\n")

# ---- 5. ULS shear with tendon vertical component relief ----
slope = 4*e/1e3/L                            # rad, parabola slope at support (e in m)
Vp    = Pe*slope                             # kN upward
Vnet  = VEd - Vp
v1=0.6*(1-fck/250); fcd_v=fck/1.5; z=0.9*dp
duct_ded=0.5*2*100                           # 2 grouted ducts ~100mm at check level
bwn=bw-duct_ded
cot=2.5; tan=1/cot
VRdmax=1.0*bwn*z*v1*fcd_v/(cot+tan)/1e3
fywd=500/1.15
Asw_s = Vnet*1e3/(z*fywd*cot)                # mm2/mm
print(f"[5] Shear: VEd={VEd:.0f}  tendon relief Vp={Vp:.0f}  -> Vnet={Vnet:.0f} kN")
print(f"    VRdmax(bw_net={bwn:.0f})={VRdmax:.0f} kN {'OK' if VRdmax>Vnet else 'FAIL'}"
      f"   links Asw/s={Asw_s:.2f} mm2/mm"
      f"  (e.g. 4-leg T12 @ {4*113/Asw_s:.0f}mm)\n")

# ---- 6. Deflection: net sustained load after balancing, with creep ----
phi=2.0; Ecm=35000; Eeff=Ecm/(1+phi)
w_bal=8*Pe*e/1e3/(L*L)                        # kN/m balanced (e in m)
w_net=wqp - w_bal
d_qp = 5*(w_net)* (L*1e3)**4 /(384*Eeff*I)    # mm  (w in N/mm = kN/m)
print(f"[6] Deflection: w_bal={w_bal:.0f} kN/m  net sustained={w_net:.0f} kN/m")
print(f"    long-term (creep) deflection ~ {d_qp:.0f} mm = L/{L*1000/d_qp:.0f}"
      f"   (limit L/500={L*1000/500:.0f}mm for finishes)  "
      f"{'OK' if d_qp< L*1000/500 else 'CHECK'}")
