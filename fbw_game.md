# Fly-By-Wire Bombing Mission

This project is an interactive Python / Pygame flight simulation demonstrating how a simple **Fly-By-Wire (FBW) attitude control system** stabilizes an unstable aircraft during a bombing mission inside hostile airspace.

The full simulation includes:

- An unstable aircraft model  
- A stabilizing FBW system  
- Bomb trajectory physics  
- Heat-seeking missiles  
- Flares and spoofing  
- HUD, artificial horizon, top-view map, rear-view camera  
- FBW on/off comparison  

---

# FBW Block Diagram

**FBW Block Diagram Used for Code**

![FBW Block Diagram](/images/fbw_diagram.jpg)

Below, each block in the diagram is mapped directly to the Python code that implements it.

---

# 1. Pilot Input

The pilot controls:

- Pitch & roll commands (Arrow keys)  
- Bomb release (**B**)  
- Flare deployment (**F**)  
- FBW toggle (**C**)  
- Quit (**ESC**)  

### Keyboard Event Handling

```python
for event in pygame.event.get():
    if event.type == pygame.QUIT:
        self.running = False
    if event.type == pygame.KEYDOWN and not self.game_over:
        if event.key == pygame.K_b and self.bombs_remaining > 0:
            sensors = self.aircraft.get_sensors()
            new_bomb = Bomb(sensors["x"], sensors["y"], sensors["h"], sensors["psi"])
            self.bombs.append(new_bomb)
            self.bombs_remaining -= 1
            self.bomb_away_time = self.sim_time

        if event.key == pygame.K_c:
            self.fbw_on = not self.fbw_on
            self.fbw_mode_text = "FBW: ON (attitude hold)" if self.fbw_on else "FBW: OFF (direct surfaces)"

        if event.key == pygame.K_f:
            sensors = self.aircraft.get_sensors()
            self.deploy_flares(sensors)

keys = pygame.key.get_pressed()
if keys[pygame.K_ESCAPE]:
    self.running = False
```

### Mapping Keys → Attitude Commands (FBW ON)

```python
def update_commands_from_keys(self, keys, dt):
    cmd_rate_deg = 20.0

    if keys[pygame.K_UP]:
        self.theta_cmd += math.radians(cmd_rate_deg) * dt
    if keys[pygame.K_DOWN]:
        self.theta_cmd -= math.radians(cmd_rate_deg) * dt

    if keys[pygame.K_RIGHT]:
        self.phi_cmd += math.radians(cmd_rate_deg) * dt
    if keys[pygame.K_LEFT]:
        self.phi_cmd -= math.radians(cmd_rate_deg) * dt

    max_cmd = math.radians(30.0)
    self.theta_cmd = max(-max_cmd, min(max_cmd, self.theta_cmd))
    self.phi_cmd   = max(-max_cmd, min(max_cmd, self.phi_cmd))
```

---

# 2. Attitude Control Law

A PD FBW controller that stabilizes pitch and roll.

### Control Gains

```python
KP_THETA = 2.5
KD_THETA = 1.2
KP_PHI   = 2.5
KD_PHI   = 1.2
MAX_DEFLECTION = 1.0
```

### PD Control Logic

```python
def compute_controls(self, sensors):
    theta = sensors["theta"]
    q     = sensors["q"]
    phi   = sensors["phi"]
    p     = sensors["p"]

    e_theta = self.theta_cmd - theta
    delta_e = KP_THETA * e_theta - KD_THETA * q

    e_phi   = self.phi_cmd - phi
    delta_a = KP_PHI * e_phi - KD_PHI * p

    delta_e = max(-MAX_DEFLECTION, min(MAX_DEFLECTION, delta_e))
    delta_a = max(-MAX_DEFLECTION, min(MAX_DEFLECTION, delta_a))

    return delta_e, delta_a
```

---

# 3. Actuator Command Block

Chooses whether FBW or direct pilot input drives the control surfaces.

```python
if self.fbw_on:
    self.controller.update_commands_from_keys(keys, DT)
    sensors = self.aircraft.get_sensors()
    delta_e, delta_a = self.controller.compute_controls(sensors)
else:
    delta_e = 0.0
    delta_a = 0.0
    step = 0.5
    if keys[pygame.K_UP]:
        delta_e += step
    if keys[pygame.K_DOWN]:
        delta_e -= step
    if keys[pygame.K_RIGHT]:
        delta_a += step
    if keys[pygame.K_LEFT]:
        delta_a -= step

    delta_e = max(-MAX_DEFLECTION, min(MAX_DEFLECTION, delta_e))
    delta_a = max(-MAX_DEFLECTION, min(MAX_DEFLECTION, delta_a))
```

FBW OFF makes the aircraft nearly unflyable because of disturbances.

---

# 4. Aircraft Dynamics

The dynamic model simulates:

- Pitch & roll motion  
- Disturbances  
- Heading from roll  
- Altitude from pitch  
- Constant forward speed  

### Dynamics Implementation

```python
def step(self, delta_e, delta_a, dt):
    # pitch dynamics
    theta_dot = self.q
    q_dot = -PITCH_DAMP * self.q + PITCH_GAIN * delta_e
    q_dot += PITCH_DISTURB * random.gauss(0.0, 1.0)

    # roll dynamics
    phi_dot = self.p
    p_dot = -ROLL_DAMP * self.p + ROLL_GAIN * delta_a
    p_dot += ROLL_DISTURB * random.gauss(0.0, 1.0)

    # integrate
    self.theta += theta_dot * dt
    self.q     += q_dot    * dt
    self.phi   += phi_dot  * dt
    self.p     += p_dot    * dt

    # altitude from pitch
    h_dot = FORWARD_SPEED * math.sin(self.theta)
    self.h += h_dot * dt
    self.h = max(self.h, 0)

    # heading from roll
    psi_dot = -TURN_GAIN * self.phi
    self.psi += psi_dot * dt

    # ground track
    self.x += FORWARD_SPEED * math.sin(self.psi) * dt
    self.y += FORWARD_SPEED * math.cos(self.psi) * dt
```

---

# 5. Inertial Motion Sensors

The IMU provides clean measurements to the controller, missiles, flares, and HUD.

```python
def get_sensors(self):
    return {
        "theta": self.theta,
        "q":     self.q,
        "phi":   self.phi,
        "p":     self.p,
        "h":     self.h,
        "x":     self.x,
        "y":     self.y,
        "psi":   self.psi,
    }
```

---

# 6. Aircraft Response (Displays & HUD)

This includes:

- Artificial horizon  
- Rear-view missile camera  
- Top-view ground map  
- Bomb impact visualization  
- HUD text  

### Example: Horizon

```python
def draw_horizon(self, theta, phi):
    ...
    pygame.draw.line(self.screen, (255, 255, 255), pts_screen[0], pts_screen[1], 3)
```

### Example: Top Map

```python
pygame.draw.rect(self.screen, (200, 50, 50), (tl[0], tl[1], rect_w, rect_h), 2)
pygame.draw.circle(self.screen, (0, 255, 0), (ax, ay), 4)
```

### Example: HUD

```python
lines = [
    f"Pitch: {theta_deg:.1f} deg   Cmd: {math.degrees(self.controller.theta_cmd):.1f} deg",
    f"Roll : {phi_deg:.1f} deg   Cmd: {math.degrees(self.controller.phi_cmd):.1f} deg",
    f"Alt  : {h:.1f} m",
    ...
]
```

---

# Gameplay Overview

**Gameplay**

![Gameplay](/images/fbw_game.jpg)

### 1. Ingress
Fly stably with FBW ON while disturbances try to destabilize you.

### 2. Bombing Run
Press **B** over the target square.

### 3. Missile Threat
Heat-seeking missiles spawn behind you.

### 4. Flares
Press **F** to deploy flares that spoof missiles.

### 5. FBW Comparison
Toggle **C** to see how unstable the aircraft becomes without FBW.

---

# Running the Simulation

[**Full code here**](fbw_code.md)

Requirements:

- Python 3  
- pygame  
- numpy  

---

# Summary

This simulation demonstrates:

- How an FBW attitude controller stabilizes an unstable aircraft  
- How feedback loops map to real physics  
- How bombs, missiles, and flares can be integrated into an aerospace simulation  
- How control engineering, avionics, and gameplay logic come together in one project  

---

