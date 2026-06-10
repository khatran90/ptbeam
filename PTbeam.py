import streamlit as st
import math

# 1. CẤU HÌNH GIAO DIỆN WEB
st.set_page_config(page_title="EC2 Post-Tensioned Beam", layout="wide")
st.title("🌉 Post-Tensioned Simply-Supported Transfer Beam — EC2 BY TWA")
st.caption("Conventions: SS EN 1992 / UK NA")
st.markdown("---")

# 2. TẠO THANH SIDEBAR NHẬP SỐ LIỆU ĐỘNG (INPUTS)
st.sidebar.header("📊 DESIGN INPUTS")

with st.sidebar.expander("📐 Geometry & Span", expanded=True):
    L = st.number_input("Span length L (m)", value=17.8, step=0.1)
    h = st.number_input("Overall depth h (mm)", value=1100.0, step=50.0)
    bw = st.number_input("Web width bw (mm)", value=900.0, step=50.0)
    hf = st.number_input("Flange thickness hf (mm)", value=250.0, step=10.0)
    beff = st.number_input("Effective flange width beff (mm)", value=3000.0, step=100.0)
    y_t = st.number_input("Tendon centroid above soffit y_t (mm)", value=180.0, step=10.0)

with st.sidebar.expander("💥 Loading Parameters", expanded=True):
    Gk = st.number_input("Permanent Load Gk (kN/m)", value=200.0, step=5.0)
    Qk = st.number_input("Variable Load Qk (kN/m)", value=35.0, step=5.0)
    psi1 = 0.5
    psi2 = 0.3
    gG, gQ = 1.35, 1.5

with st.sidebar.expander("🧪 Materials & Strands", expanded=True):
    fck = st.number_input("Concrete Strength fck (MPa)", value=40.0, step=5.0)
    fckt = st.number_input("Transfer Strength fckt (MPa)", value=32.0, step=5.0)
    acc = 0.85
    fctm = round(0.30 * fck**(2/3), 2)
    
    st.markdown("**PT Strand (15.7mm):**")
    Ap1 = 150.0 # mm2
    fpk = 1860.0
    fp01k = 1600.0
    gp = 1.15
    fpd = fp01k / gp
    sig_pm0 = 1300.0
    sig_pe = 1150.0

# 3. THUẬT TOÁN LOGIC TÍNH TOÁN KẾT CẤU
P0_1 = sig_pm0 * Ap1 / 1e3   
Pe_1 = sig_pe * Ap1 / 1e3    

wEd = gG*Gk + gQ*Qk
wqp = Gk + psi2*Qk
wchr = Gk + Qk
w_sw = bw*h*25e-6           

MEd = wEd *L*L/8
VEd = wEd *L/2
Mqp = wqp *L*L/8
Mchr = wchr*L*L/8
Msw = w_sw*L*L/8

# Section properties
A_w = bw*h ; A_f = (beff-bw)*hf ; A = A_w + A_f
yb = (A_w*h/2 + A_f*(h-hf/2))/A            
Iw = bw*h**3/12 + A_w*(yb-h/2)**2
If = (beff-bw)*hf**3/12 + A_f*((h-hf/2)-yb)**2
I = Iw + If
Zb = I/yb ; Zt = I/(h-yb)
kb = Zb/A
e = yb - y_t                              
dp = h - y_t                               

# [1] SLS Prestress Req
fcd = acc*fck/1.5
P_req = Mqp*1e6 /(kb+e)          
n_sls = P_req/1e3/Pe_1

# [2] ULS Flexure
def MRd(n_strands):
    T = n_strands*Ap1*fpd/1e3                       
    x_depth = T*1e3/(fcd*beff)                    
    if x_depth > hf:                                
        Tf = fcd*beff*hf/1e3
        x_depth = hf + (T-Tf)*1e3/(fcd*bw)
    z_arm = dp - 0.4*x_depth
    return T*z_arm/1e3, x_depth

n = math.ceil(n_sls)
while MRd(n)[0] < MEd: n += 1
M_R, x = MRd(n)

Pe = n*Pe_1 ; P0 = n*P0_1

# [3] Transfer stresses
Aw=bw*h; Zw=bw*h**2/6; ew=h/2-y_t
def stress(P, M, Aa, Ztt, Zbb, ecc):
    sb = P*1e3/Aa + P*1e3*ecc/Zbb - M*1e6/Zbb
    st = P*1e3/Aa - P*1e3*ecc/Ztt + M*1e6/Ztt
    return st, sb
st_trans, sb_trans = stress(P0, Msw, Aw, Zw, Zw, ew)
lim_c = 0.6*fckt ; lim_t = -fctm

aT = (lim_t - Msw*1e6/Zw)/(P0*1e3*(1/Aw - ew/Zw))   
aC = (lim_c + Msw*1e6/Zw)/(P0*1e3*(1/Aw + ew/Zw))   
amax = min(aT, aC)

# [4] Service stresses
st_serv, sb_serv = stress(Pe, Mchr, A, Zt, Zb, e)

# [5] ULS Shear
slope = 4*e/1e3/L                            
Vp = Pe*slope                             
Vnet = VEd - Vp
v1 = 0.6*(1-fck/250); fcd_v = fck/1.5; z_shear = 0.9*dp
bwn = bw - (0.5*2*100)                           
cot = 2.5; tan = 1/cot
VRdmax = 1.0*bwn*z_shear*v1*fcd_v/(cot+tan)/1e3
fywd = 500/1.15
Asw_s = Vnet*1e3/(z_shear*fywd*cot)                

# [6] Deflection
phi=2.0; Ecm=35000; Eeff=Ecm/(1+phi)
w_bal = 8*Pe*e/1e3/(L*L)                        
w_net = wqp - w_bal
d_qp = 5*(w_net)* (L*1e3)**4 /(384*Eeff*I)    

# 4. TRÌNH BÀY KẾT QUẢ RA GIAO DIỆN WEB (UI/UX)
col1, col2, col3 = st.columns(3)
col1.metric("Design Load wEd", f"{round(wEd, 1)} kN/m")
col2.metric("ULS Moment MEd", f"{int(MEd)} kNm")
col3.metric("Required PT Strands", f"{n} Strands", f"SLS basis: {round(n_sls,1)}")

st.markdown("---")

c_left, c_right = st.columns(2)

with c_left:
    st.subheader("💡 1. Section & Flexural Capacity")
    st.write(f"**Section Area:** {int(A)} $mm^2$ | **Centroid $y_b$:** {int(yb)} mm")
    st.write(f"**Moment Capacity $M_{{Rd}}$:** {int(M_R)} kNm")
    if M_R >= MEd:
        st.success(f"✅ ULS Flexure OK (x/dp = {round(x/dp, 2)} ≤ 0.45 ductile)")
    else:
        st.error("❌ ULS Flexure FAIL")
        
    st.subheader("🚧 2. Transfer Stage Check (Web Only)")
    st.write(f"**Top Stress:** {round(st_trans, 2)} MPa (Limit: {round(lim_t, 2)})")
    st.write(f"**Bottom Stress:** {round(sb_trans, 2)} MPa (Limit: {round(lim_c, 2)})")
    if st_trans < lim_t or sb_trans > lim_c:
        st.warning(f"⚠️ Staging Required! Max safe initial prestress: **{round(amax*100, 0)}%** of $P_0$")
    else:
        st.success("✅ Transfer stresses OK without staging")

with c_right:
    st.subheader("✂️ 3. ULS Shear Verification")
    st.write(f"**Tendon Relief $V_p$:** {int(Vp)} kN ➡️ **Net Shear $V_{{net}}$:** {int(Vnet)} kN")
    st.write(f"**Max Concrete Strut Capacity $V_{{Rd,max}}$:** {int(VRdmax)} kN")
    if VRdmax > Vnet:
        st.success(f"✅ Concrete Strut OK. Required links $A_{{sw}}/s$: **{round(Asw_s, 2)}** $mm^2/mm$")
    else:
        st.error("❌ Concrete Strut Crushing Failure! Increase web width or concrete grade.")

    st.subheader("📉 4. Long-term Deflection (Creep)")
    st.write(f"**Balanced Load:** {round(w_bal, 1)} kN/m | **Net Sustained:** {round(w_net, 1)} kN/m")
    limit_defl = L * 1000 / 500
    if d_qp < limit_defl:
        st.success(f"✅ Deflection OK: {round(d_qp, 1)} mm (Limit L/500 = {round(limit_defl, 1)} mm)")
    else:
        st.error(f"❌ Deflection CRITICAL: {round(d_qp, 1)} mm (Exceeds {round(limit_defl, 1)} mm)")
