# LION 🌕
**L**unar **I**njection **O**ptimization and **N**avigation

A comprehensive MATLAB-based simulation and optimization tool for Trans-Lunar Injection (TLI) maneuvers and Earth-Moon trajectory analysis. This project numerically integrates spacecraft trajectories under the gravitational influence of both the Earth and the Moon, optimizing departure conditions to minimize fuel consumption while providing detailed 3D visualizations of the mission.

## 🚀 Features
* **TLI Fuel Optimization:** Automatically calculates the optimal fraction of escape velocity required for Trans-Lunar Injection to minimize fuel usage, utilizing the Tsiolkovsky rocket equation and MATLAB's `fminbnd` bounded optimization.
* **N-Body Numerical Integration:** Propagates the spacecraft's state vector (position and velocity) using the high-order `ode113` solver, accounting for Earth's central gravity and the Moon's gravitational pull (including indirect terms).
* **Harmonic Lunar Ephemeris:** Incorporates a custom harmonic series-based ephemeris model to accurately predict the Moon's position and velocity vectors at any given Julian Date.
* **Advanced Reference Frames:** Computes and visualizes trajectories not only in the Earth-Centered Inertial (ECI) frame but also applies Direction Cosine Matrices (DCM) to transform and plot the data in a Moon-fixed rotating reference frame.
* **Rich 3D Visualizations:** Includes fully animated 3D plots of the spacecraft and Moon trajectories, alongside static visualizations tracking orbital inclination variations, perilune approach, and fuel consumption profiles over time.
* **Astronomical State Tracking:** Continuously monitors and outputs critical mission parameters such as Right Ascension (RA), Declination (Dec), Flight Path Angle, and target point errors upon lunar arrival.

## 🛠️ Technology Stack
This project is developed entirely in **MATLAB** and relies on its core mathematical and differential equation solving capabilities (specifically `ode113` and `fminbnd`). No external toolboxes are strictly required for the core physics engine.

## 📂 Repository Structure
* `LION.m` - The main execution script. Handles parameter initialization, optimization, integration, and triggers all plots.
* `\Functions` - Contains all auxiliary functions.
  * `rates.m` - Contains the differential equations of motion defining the N-body physics for the ODE solver.
  * `simpsonsLunarEphemeris.m` - Calculates the Moon's state vector based on harmonic coefficients and Julian Date.
  * **Astrodynamics & Math Functions:**
    * `launchConditions2RV.m` - Converts launch angles and altitude into initial ECI position and velocity vectors.
    * `RAandDecFromR.m` - Extracts Right Ascension and Declination from Cartesian state vectors.
    * `julianDay.m` - Converts standard calendar dates and Universal Time (UT) into Julian Days.
  * **Propulsion & Optimization:**
    * `computeTotalDeltavandfuel.m` / `computeFuelUsed.m` - Analyzes ΔV requirements and mass properties.
    * `fuelUsedFac.m` - The objective function used by the optimizer to minimize fuel consumption.
  * **Visualization Tools:**
    * `plotAnimatedXYZ.m` - Renders the animated 3D motion of the Earth, Moon, and spacecraft in the ECI frame.
    * `plotXYZ.m` - Plots the trajectory relative to a Moon-fixed rotating coordinate system.
    * `plotFuelvsTime.m` - Graphs the cumulative fuel consumption profile during the main engine burn.

## ⚙️ Installation & Usage
Since this project is developed entirely in MATLAB, no external compilation or complex dependency management is required.

1. **Clone the repository:**
   ```bash
   git clone https://github.com/BFG508/LION
2. **Open MATLAB** and navigate to the cloned `LION` directory.
3. **Add to Path**: Ensure that `/Functions` subdirectory is added to your MATLAB path to allow the main scripts to access necessary functions and data files. You can do this by right-clicking the `/Functions` folder in the Current Folder browser and selecting Add to Path > Selected Folders and Subfolders.
4. **Run the Simulation**: Open `LION.m` and run the script. The console will output the optimization results, departure conditions, and lunar arrival data, followed by the generation of the 3D trajectory animations and analytical plots.