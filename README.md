# Global-Geodesic-Router-Web-Mapper

I built this web application to tackle coordinate transformations and geodetic analysis by computing the exact, mathematically shortest paths (geodesics) between Tehran and other major global capitals. Rather than relying on simple, flat-map approximations, this project utilizes a custom grid search algorithm to evaluate spatial discrepancies and navigate the Earth's true ellipsoidal shape. It’s a practical, interactive way to visualize complex spatial navigation and distance calculations on the web.

## What's Inside?

* **Geodetic Grid Search:** A custom algorithm using Runge-Kutta numerical integration (`scipy.integrate.solve_ivp`) to pinpoint exact azimuths and spatial trajectories.


* **WGS84 Ellipsoid Math:** Calculates the exact radii of curvature using Earth's semi-major axis ($a = 6378137$ meters) and eccentricity ($e \approx 0.0818$) for maximum geodetic precision.


* **Interactive Web Mapper:** A Flask-driven backend paired with a Folium front-end, plotting the calculated shortest paths dynamically on an interactive "Dark Matter" map.



## Tech Stack

* `flask` (For the web server and application routing)


* `numpy` & `scipy` (For the heavy lifting: mathematical arrays and solving differential equations)


* `folium` (For rendering the sleek, interactive map visualizations)



## How to Run It

1. Clone this repository to your local machine.
2. Install the necessary dependencies via your terminal:
`pip install flask numpy scipy folium`
3. Ensure you have an `index.html` file set up in your `templates/` folder to handle the front-end rendering.


4. Run the application: `python main.py` (or the specific name of your script).


5. Open your web browser and navigate to the localhost address provided (typically `[http://127.0.0.1:5000/](http://127.0.0.1:5000/)`). Select a target capital city from the interface to instantly generate and map the geodesic path from Tehran!
