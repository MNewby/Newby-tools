import astro_coordinates as co
import math as ma

deg = 180.0/ma.pi

# milkyway thetas and phis
mw_t = [-1.50, -1.53, -1.31, -1.70, -1.91, -1.91, -2.33, -1.92, -1.97, -2.40, -2.36, -2.62, 0.60, -2.20]
mw_p = [-0.06, -0.05, -0.01, 0.15, -0.22, -0.14, -0.27, -0.22, -0.32, -0.10, -0.19, -0.30, -1.20, -2.10]

# BlueGene thetas and phis
bg_t = [1.1, 1.2, 1.3, 1.5, 1.8, 1.9, 2.3, 2.0, 2.1, 2.4, 2.3, 2.6, 2.4, 2.7]
bg_ts = [0.5, 1.3, 0.5, 0.2, 0.8, 0.7, 0.7, 0.8, 0.7, 1.6, 0.3, 0.9, 1.0, 1.6]
bg_p = [3.0, -2.9, 3.1, -3.0, 2.9, 3.0, 2.9, 2.9, 3.0, 3.0, 3.0, 3.0, 1.9, 1.4]
bg_ps = [2.0, 0.5, 0.3, 1.0, 0.3, 0.6, 0.5, 0.3, 1.0, 2.0, 0.7, 2.0, 2.5, 1.6]

for i in range(len(mw_t)):
    print bg_t[i], bg_ts[i], bg_p[i], bg_ps[i], mw_t[i], mw_p[i]
    x, y = co.angle_bounds3(mw_t[i]*deg, mw_p[i]*deg, phi_max=180.0)
    theta, phi = x/deg, y/deg
    dt, dp = abs(theta-bg_t[i])/bg_ts[i], abs(phi-bg_p[i])/bg_ps[i]
    print theta, dt, phi, dp, "\n"
