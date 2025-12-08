# 3D Missile vs Aircraft Engagement Simulator

This project is an interactive 3D simulation of a heat seeking missile chasing an aircraft over a 50 km by 50 km airspace. It runs in Python using NumPy and Matplotlib and animates the full engagement in real time.

The original open source code provided basic missile tracking and 3D plotting. I extended it into a full combat-style scenario with radar, different aircraft types, bombs, flares, altitude changes, crash and escape conditions, and a persistent scoreboard.

---

## Screenshots

> Place screenshots inside an `images/` folder.  
> Recommended screenshots are described below under **Suggested Images to Capture**.

### Main Engagement View
![3D engagement overview](images/sim_overview.png)

### Bomb Impact Example
![Bomb hit near enemy base](images/bomb_hit.png)

### Successful Flare Spoof
![Flares spoofing missile](images/flares_spoof.png)

---

## Technologies

- Python  
- NumPy for vector and trajectory computation  
- Matplotlib 3D + `FuncAnimation` for real time visualization  

---

## Scenario Overview

Each engagement simulates:

1. An aircraft spawns on a random edge of the map and flies toward a randomly located enemy base.  
2. Radar attempts to detect the aircraft each simulated second using a per aircraft probability.  
3. Upon detection, a heat seeking missile launches and begins guiding toward the aircraft.  
4. When near the base, the aircraft drops a bomb that falls vertically and may hit the target.  
5. After bomb release, the aircraft switches to a more aggressive escape maneuver with turns and altitude modulation.  
6. The engagement ends if:
   - The missile hits the aircraft.  
   - The aircraft escapes the XY airspace.  
   - The aircraft crashes by dropping below 500 m altitude.  
7. The scoreboard updates and a new engagement begins automatically.

---

## My Contributions

The original repository only contained basic missile tracking and 3D plotting.  
I implemented all of the following systems.

---

### 1. Multiple Aircraft Types

Created a structured `AIRCRAFT_TYPES` system defining:

- Color and aircraft identity  
- Speed ranges (`speed_min`, `speed_max`)  
- Maneuverability (`curve_min`, `curve_max`, `max_turn_rate`)  
- Radar detection probability (`radar_detection_prob`)  
- Flare behavior (`flare_trigger_distance`, `flare_success_prob`)  

Aircraft currently implemented:

- F22  
- A10  
- B2  
- X15  

Each type flies differently, is detected differently, and has unique flare performance.

---

### 2. Radar Detection and Missile Launch Logic

In the original script the missile simply chased the aircraft immediately.

My implementation adds:

- A radar sweep every simulated second  
- Probability based detection depending on the aircraft type  
- Missile launch only when detection occurs  
- A radar indicator that changes from green to red when radar lock is achieved  

This models realistic detectability—for example, the B2 has the lowest detection probability.

---

### 3. Bomb and Enemy Base Mechanics

I added a fully new subsystem:

- A random enemy base appears between 10 km and 40 km in both X and Y  
- The aircraft flies toward the base using constrained heading changes  
- When within 50 m of the base, the aircraft drops a bomb  
- The bomb falls straight down and explodes upon hitting the ground  
- A hit is recorded if the explosion is within 200 m of the base  

If a hit occurs, `payload_recorded` becomes true and later increments the global score.

---

### 4. Post Bomb Escape Maneuver and Altitude Changes

After dropping the bomb, the aircraft transitions into a complex evasive maneuver:

- The escape path is composed of:  
  **straight → turn → straight → turn → straight**
- Turn durations and straight intervals come from aircraft agility ranges  
- Heading uses a dynamic angle `psi` updated each timestep  
- Altitude changes follow cosine shaped vertical bumps synchronized with turns  
- Altitude is clamped to the airspace constraints  

This gives each aircraft a unique, realistic escape pattern.

---

### 5. Flares and Missile Spoofing Logic

Added full countermeasure behavior:

- When the missile is within `flare_trigger_distance`, flares deploy  
- A probability check determines success  
- If flares **fail**:
  - Missile continues tracking the aircraft  
- If flares **succeed**:
  - Missile locks onto the flare location  
  - Missile is frozen in its spoofed position  
  - A large flare marker appears in the 3D view  

Outcome logic incorporates flare success in determining whether the aircraft can later escape or crash.

---

### 6. Crash and Escape Detection

Added two ways for engagements to end besides a missile hit:

- **Escape**: aircraft leaves `[0, 50000]` in X or Y  
- **Crash**: aircraft altitude drops below 500 m  

Crash and escape indices are calculated after trajectory generation and marked visually.

---

### 7. Outcome Logic + Scoreboard

Instead of only detecting missile hits, I added a full decision system:

The possible outcomes:

- **missile** – missile hits aircraft  
- **aircraft** – aircraft escapes airspace  
- **crash** – aircraft crashes  
- **payload delivered** – bomb successfully hits base  

Global counters track:

- `missile_neutralized_target`  
- `aircraft_escaped`  
- `payload_delivered`  

Displayed live on the animation window.

---

### 8. Visualization Enhancements

I added:

- Aircraft start marker  
- Missile start marker  
- Base marker  
- Flares marker  
- Bomb explosion marker  
- Crash marker  
- Intercept marker  
- Timeline, speed, distance, and aircraft-type text  
- Full scoreboard overlay  

I also added adaptive frame stepping so the simulation stays smooth even with dense time resolution.

---

## Code Structure

### `rotate_towards(current, desired, max_angle)`
Clamps heading changes so aircraft cannot instantly rotate.

### `generate_new_engagement()`
Responsible for:

- Resetting all state  
- Randomizing aircraft, base, bomb, missile, speeds  
- Computing full aircraft, missile, and bomb trajectories  
- Running radar, flare, bomb, and post bomb logic  
- Evaluating hit, escape, crash, and payload success  
- Updating markers and legend  

### `init()`
Initializes all Matplotlib artists.

### `update(frame)`
Runs the animation:

- Moves aircraft, missile, bomb markers  
- Draws trails  
- Shows flares, explosions, crashes, intercepts  
- Updates text and scoreboard  
- Stops movement after outcome time  
- Automatically generates a new engagement after delay  

---

## How To Run

```bash
pip install numpy matplotlib
python missile_simulation.py
