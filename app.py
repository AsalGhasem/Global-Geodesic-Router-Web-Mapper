from flask import Flask, render_template, request
import numpy as np
from scipy.integrate import solve_ivp
from math import *
import folium

app = Flask(__name__)

phis = [-35.26666667, 48.2, -15.78333333, 45.41666667, 39.91666667, 30.05, 60.16666667, 35.7, 55.75, -25.7]
lams = [149.133333, 16.366667, -47.916667, -75.7, 116.383333, 31.25, 24.933333, 51.416667, 37.6, 28.216667]
caps = ["Canberra", "Vienna", "Brasilia", "Ottawa", "Beijing", "Cairo", "Helsinki", "Tehran", "Moscow", "Pretoria"]
country_map = {
    "Canberra": "Australia", "Vienna": "Austria", "Brasilia": "Brazil", 
    "Ottawa": "Canada", "Beijing": "China", "Cairo": "Egypt", 
    "Helsinki": "Finland", "Tehran": "Iran", "Moscow": "Russia", "Pretoria": "SouthAfrica"
}

phis_rad = [np.deg2rad(phi) for phi in phis]
lams_rad = [np.deg2rad(lam) for lam in lams]
tehran_index = caps.index("Tehran")
tehran_phi, tehran_lam = phis_rad[tehran_index], lams_rad[tehran_index]
e, a, R_G = 0.0818191908426215, 6378137, 6371000 # in meters

def alpha0_calculator(phi1, lam1, phi2, lam2):
    return np.arctan2(np.sin(lam2-lam1), ((np.cos(phi1)*np.tan(phi2)) - (np.sin(phi1)*np.cos(lam2-lam1))))

def s0_calculator(phi1, lam1, phi2, lam2):
    return R_G*np.arccos((np.sin(phi1)*np.sin(phi2)) + (np.cos(phi1)*np.cos(phi2)*np.cos(lam2-lam1)))

def M(phi): 
    return ((a * (1 - e**2)) / ((1 - (e**2 * (sin(phi)**2)))) ** (3 / 2))

def N(phi): 
    return ((a) / ((1 - (e**2 * (sin(phi)**2)))) ** (1 / 2))

def geodesic(s, y):
    phi, lam, alpha = y
    return [np.cos(alpha)/M(phi), np.sin(alpha)/(N(phi)*np.cos(phi)), np.sin(alpha)*np.tan(phi)/N(phi)]

def grid_search(alpha0, da_deg, alpha_sq_num, s0, ds, s_sq_num, phi_city, lam_city):
    da_rad      = np.deg2rad(da_deg)
    alpha_range = np.linspace(alpha0 - da_rad, alpha0 + da_rad, alpha_sq_num)
    s_range     = np.linspace(s0-ds, s0+ds, s_sq_num)
    best_error, best_phi, best_lam, best_alpha, best_s = np.inf, None, None, None, None

    for alpha in alpha_range:
        for s in s_range:
            sol = solve_ivp(geodesic, (0, s), [tehran_phi, tehran_lam, alpha], rtol=1e-10, atol=1e-10, method="RK45")
            p_e, l_e = sol.y[0, -1], sol.y[1, -1]
            err = np.sqrt(((p_e - phi_city)**2) + ((l_e - lam_city)**2))
            if err < best_error:
                best_error, best_phi, best_lam, best_alpha, best_s = err, p_e, l_e, alpha, s
    return np.rad2deg(best_phi), np.rad2deg(best_lam), best_alpha, best_s, best_error

# calculate data for the app
RESULTS = {}
for r in range(len(phis)):
    if caps[r] == "Tehran": continue
    # round 1
    a0 = alpha0_calculator(tehran_phi, tehran_lam, phis_rad[r], lams_rad[r])
    s0 = s0_calculator    (tehran_phi, tehran_lam, phis_rad[r], lams_rad[r])
    b_p, b_l, b_a, b_s, b_e = grid_search(a0, 2, 40, s0, 5000, 10, phis_rad[r], lams_rad[r])
    # round 2 (final answer)
    res = grid_search(b_a, 2, 40, b_s, 5000, 10, phis_rad[r], lams_rad[r])
    RESULTS[caps[r]] = {
        "phi": res[0], "lam": res[1], "azimuth": np.rad2deg(res[2]), 
        "distance": res[3], "err": res[4], "country": country_map[caps[r]]
    }

@app.route('/', methods=['GET', 'POST'])
def index():
    selected_cap = request.form.get('capital')
    data = RESULTS.get(selected_cap)
    
    # create Map
    m = folium.Map(location=[20, 40], zoom_start=2, tiles="cartodb dark_matter")
    if data:
        start = [phis[tehran_index], lams[tehran_index]]
        end = [phis[caps.index(selected_cap)], lams[caps.index(selected_cap)]]

        best = RESULTS[selected_cap]
        best_alpha = np.deg2rad(best["azimuth"])
        best_s = best["distance"]

        sol = solve_ivp(
            geodesic,
            (0, best_s),
            [tehran_phi, tehran_lam, best_alpha],
            max_step=best_s/200
        )

        path = list(zip(np.rad2deg(sol.y[0]), np.rad2deg(sol.y[1])))
        folium.PolyLine(path, color="#2d5a27", weight=3, opacity=0.8).add_to(m)

        folium.Marker(start, popup="Tehran", icon=folium.Icon(color='green')).add_to(m)
        folium.Marker(end,   popup=selected_cap, icon=folium.Icon(color='lightgray')).add_to(m)

    map_html = m._repr_html_()
    return render_template('index.html', caps=RESULTS.keys(), data=data, 
                           sel=selected_cap, map_html=map_html)

if __name__ == '__main__':
    app.run(debug=True)