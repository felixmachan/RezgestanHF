from sympy import symbols, expand, simplify, diff, integrate, latex
import math
from math import radians, cos, sqrt, atan, sin, degrees, log



def kerek(szam):
    return round(szam, 6)

# =================== ADATOK ================== #

d = 0.02
ro = 7800
m = 0.08
e = 0.4
l_1 = 0.18
l_2 = 0.12
l_3 = 0.12
l_4 = 0.4

k_1 = 52
k_2 = 110
k_t = 10.8
c_1 = 1.5
F_0 = 0
M_0 = 0
r_0 = 0.05
omega = 28

pi = math.pi
g = 9.81

# =================== 1. FELADAT ================ #
print("1. FELADAT\n")
m_1 = ro * ((d**2 * pi )/4 * (l_1 + l_2))
m_2 = ro * ((d**2 * pi )/4 * (l_3 + l_4))

G_1 = m_1 * g
G_2 = m_2 *g

theta_A = m_1 * (0.25* (d/2)**2 + (1/12)*(l_1 + l_2)**2) + m_1 * ((l_1+l_2) / 2)**2 + m_2 * (0.25 * (d/2)**2 + (1/12) *(l_3 + l_4)**2) + m_2 * (((l_3+l_4)/2)**2)
print(latex(theta_A))
print(f"m_1 = {round(m_1,4)} [kg]")
print(f"m_2 = {round(m_2,4)} [kg]")
print(f"G_1 = {round(G_1,4)} [N]")
print(f"G_2 = {round(G_2,4)} [N")
print(f"theta_A = {kerek(theta_A)} [kgm^2]")

print("\n")

# ----------------- 2. FELADAT ----------------- #


m_12 = (l_1 + l_2) * (((d**2)*pi)/4) * ro
m_34 = (l_3 + l_4) * (((d**2)*pi)/4) * ro

print(f"m_12 = {kerek(m_12)} [kg]")
print(f"m_34 = {kerek(m_34)} [kg]")

omega_n = sqrt((k_2*l_1**2 + k_1*l_3**2 +k_t + m_12 * g * (l_1 + l_2) *0.5) / theta_A)
print(f"omega_n = {kerek(omega_n)} [rad/s]")



zeta = (c_1 * (l_3 + l_4)**2 / (2* theta_A * omega_n))

print(f"zeta = {kerek(zeta)} [-]")
print("\n")

# ----------------- 3. FELADAT ----------------- #

print("3. FELADAT\n")
n = 1
A_1 = 1
A_n1 = 0.43


logdek = (1/n)*log(A_1 / A_n1)
zeta_becs =(logdek / (sqrt(4*pi**2 + logdek**2)))
T = 0.57



omega_dbecs = (2*pi) / T
omega_becs = omega_dbecs / (sqrt(1-zeta_becs**2))

H_zeta =   100* math.fabs((zeta - zeta_becs) / zeta)
H_omega =  100* math.fabs((omega_n - omega_becs)) / omega_n

print(f"logdek = {kerek(logdek)} [-]")
print(f"zeta_becs = {kerek(zeta_becs)} [-]")
print(f"omega_dbecs = {kerek(omega_dbecs)} [rad/s]")
print(f"omega_becs = {kerek(omega_becs)} [rad/s]")

print(f"H_zeta = {kerek(H_zeta)} [%]")
print(f"H_omega = {kerek(H_omega)} [%]")

print("\n")

# --------------------- 4. FELADAT ------------------#
print("4. FELADAT\n")
c_t1 = 0
radian = 30 * pi / 180
AT = l_1 + l_2
FT = l_3 * cos(radian)

theta_F = m*l_3**2
m_red1 = theta_A / (AT**2)
m_red2 = theta_F / (FT**2)

omega_2 = math.sqrt((2*m*g*l_3*cos(radian)) / theta_F)
ct_2 = omega_2 * l_3
ct_2n = omega_2 * FT

c_S = (c_t1 * m_red1 + ct_2n * m_red2) / (m_red1 + m_red2)
v_S = c_S
v_t1 = (1+e) * v_S

v_sn = v_S + v_S*e
omega_t = v_sn / (l_1 + l_2)

fipont_0 = -v_t1 / AT
omega_t1 = fipont_0



print(f"AT = {round(AT,4)} [m]")
print(f"FT = {round(FT,4)} [m]")
print(f"theta_F = {round(theta_F,4)} [kgm^2]")
print(f"m_red1 = {round(m_red1,4)} [kg]")
print(f"m_red2 = {round(m_red2,4)} [kg]")
print(f"omega_2 = {round(omega_2,4)} [rad/s]")
print(f"ct_2 = {round(ct_2,4)} [m/s]")
print(f"ct_2n = {round(ct_2n,4)} [m/s]")
print(f"c_S = v_S = {round(c_S,4)} [m]")
print(f"v_t1 = {round(v_t1,4)} [m/s]")
print(f"fipont_0 = omega_t1 = {kerek(fipont_0)} [rad/s]")

print("\n")

# --------------- 5. FELADAT --------------- #
print("5. FELADAT\n")
omega_d = omega_n * math.sqrt(1-zeta**2)

lambda_1 = omega / omega_n
N = sqrt((1-lambda_1**2)**2 + 4*zeta**2 * lambda_1**2)**-1
f_0 = (k_1*l_3 * r_0) / (theta_A * omega_n**2)
fi = N*f_0
kisthetavesszo = atan((2*zeta*lambda_1) / (1- lambda_1**2))
kistheta = kisthetavesszo + pi

C_1 = -fi* sin(kistheta)
C_2 = -(omega_t + zeta *omega_n * C_1 - fi * omega * cos(-kistheta)) / omega_d

print(f"omega_d = {kerek(omega_d)} [rad/s]\n")
print(f"lambda = {kerek(lambda_1)} [-]")
print(f"N = {kerek(N)} [-]")
print(f"f_0 = {kerek(f_0)} [rad]")
print(f"fi = {kerek(fi)} [rad]")
print(f"theta' = {kerek(kisthetavesszo)} [rad]")
print(f"theta = {kerek(kistheta)} [rad]")
print(f"C_1 = {kerek(C_1)} [rad]")
print(f"C_2 = {kerek(C_2)} [rad]")





