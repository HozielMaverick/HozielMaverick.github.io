# 3D Missile vs Aircraft Engagement Simulator

This project is an interactive 3D simulation of a heat seeking missile chasing an aircraft over a 50 km by 50 km airspace. It runs in Python using NumPy and Matplotlib and animates the full engagement in real time.

The original open source code provided basic missile tracking and 3D plotting. I extended it into a full combat style scenario with radar, different aircraft types, bombs, flares, altitude changes, crash and escape conditions, probabilistic bomb accuracy, and a persistent scoreboard.

---

## Screenshots
### Main Engagement View
![3D engagement overview](/images/sim_overview.jpg)
### Missile Intercept Succesful
![Missile Intercept](/images/sim_intercept.jpg)
### Missile Intercept Before Plane Drops Bomb
![Missile Intercept Before Bomb](/images/sim_intercept_no_bomb.jpg)
### Bomb Impact Example
![Bomb hit near enemy base](/images/bomb_hit.jpg)
### Successful Flare Spoof
![Flares spoofing missile](/images/flares_spoof.jpg)
### Plane not Detected by Missile - Stealth
![Stealth](/images/stealth.jpg)

---

## Technologies

- Python  
- NumPy for vector and trajectory computation  
- Matplotlib 3D and `FuncAnimation` for real time visualization  

---

## Scenario Overview

Each engagement simulates:

1. An aircraft spawns on a random edge of the map and flies toward a randomly located enemy base.  
2. Radar attempts to detect the aircraft each simulated second using a per aircraft probability.  
3. Upon detection, a heat seeking missile launches and begins guiding toward the aircraft.  
4. When near the base, the aircraft drops a bomb that falls vertically. The bomb may or may not count as a hit depending on a per aircraft bomb hit probability.  
5. After bomb release, the aircraft switches to a more aggressive escape maneuver with turns and altitude modulation.  
6. The engagement ends if:  
   - The missile hits the aircraft.  
   - The aircraft escapes the XY airspace.  
   - The aircraft crashes by dropping below 500 m altitude.  
7. The scoreboard updates and a new engagement begins automatically. Bomb hits only contribute to the payload counter when both the geometry and the probability check indicate a successful strike.

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
- Bomb accuracy (`bomb_hit_prob`)  

Aircraft currently implemented:

- F22  
- A10  
- B2  
- X15  

Each type flies differently, is detected differently, and has unique flare and bombing performance. For example, the B2 has very low radar detection probability and a high bomb hit probability, while the X15 is much harder to bomb accurately with due to its extreme speed.

---

### 2. Radar Detection and Missile Launch Logic

In the original script the missile simply chased the aircraft immediately.

My implementation adds:

- A radar sweep every simulated second  
- Probability based detection depending on the aircraft type  
- Missile launch only when detection occurs  
- A radar indicator that changes from green to red when radar lock is achieved  

This models realistic detectability. For instance, the B2 has the lowest detection probability.

---

### 3. Bomb, Enemy Base, and Hit Probability Mechanics

I added a fully new subsystem:

- A random enemy base appears between 10 km and 40 km in both X and Y  
- The aircraft flies toward the base using constrained heading changes  
- When within 50 m of the base, the aircraft drops a bomb  
- The bomb falls straight down and reaches the ground at high speed  
- The geometric hit region is defined as within 200 m of the base in the XY plane  

On impact, the bomb logic applies a probabilistic hit model:

- Each aircraft has a `bomb_hit_prob` between 0 and 1.  
- If the bomb lands within 200 m of the base, the code draws a random number.  
- Only if this random draw is below `bomb_hit_prob` does the sim:  
  - Mark `payload_recorded = True`.  
  - Set `explosion_index`.  
  - Render the explosion marker.  
- If the probability test fails, the impact is treated as a miss, no explosion marker is shown, and the payload counter is not incremented.

This lets different aircraft feel distinct as bombers, even when they reach the same release point.

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
- If flares fail:  
  - Missile continues tracking the aircraft  
- If flares succeed:  
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

### 7. Outcome Logic and Scoreboard

Instead of only detecting missile hits, I added a full decision system.

The possible outcomes:

- **missile** – missile hits aircraft  
- **aircraft** – aircraft escapes airspace  
- **crash** – aircraft crashes  
- **payload delivered** – bomb successfully hits base, both geometrically and in the probabilistic model  

Global counters track:

- `missile_neutralized_target`  
- `aircraft_escaped`  
- `payload_delivered`  

`payload_delivered` increments only when `payload_recorded` is true and the explosion time is reached. These values are displayed live on the animation window.

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
- Timeline, speed, distance, and aircraft type text  
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
- Applying probabilistic bomb hit checks based on `bomb_hit_prob`  
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
- Increments the payload counter only once per successful bomb hit  
- Automatically generates a new engagement after a short delay  

---

## Example Aircraft Configuration

```python
AIRCRAFT_TYPES = {
    "A10": {
        "color": "orange",
        "flare_trigger_distance": 100.0,
        "flare_success_prob": 0.70,
        "bomb_hit_prob": 0.85,
        "speed_min": 2000,
        "speed_max": 2100,
        "curve_min": 5,
        "curve_max": 10,
        "max_turn_rate": 0.7,
        "radar_detection_prob": 0.50,
    },
    ...
}
```
---

## Access to full code
[Full code here](Missile_Game_code.md)
