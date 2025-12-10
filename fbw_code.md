```python

import math
import random
import pygame
import numpy as np

# ============================================================
# Global parameters
# ============================================================
DT = 0.01
FPS = 60
FORWARD_SPEED = 800.0

TOPVIEW_SPAN = 8000.0   # map shows x,y in [-8000, +8000]

TARGET_HALF_SIZE = 500
# Dynamics coefficients
PITCH_DAMP = 1.5
PITCH_GAIN = 10.0
ROLL_DAMP = 1.5
ROLL_GAIN = 12.0

TURN_GAIN = -0.5

# Random disturbances
PITCH_DISTURB = 2
ROLL_DISTURB = 2

# FBW control gains
KP_THETA = 2.5
KD_THETA = 1.2
KP_PHI = 2.5
KD_PHI = 1.2

MAX_DEFLECTION = 1.0

# Gravity
G = 9.81 * 50
BOMB_DROP_ALT = 2000

# Target corridor and airspace
TARGET_X = 8000.0
TARGET_Y = 0.0
TARGET_LENGTH = 1000.0
TARGET_WIDTH = 1000.0

AIRSPACE_RADIUS = 6000.0  # kept, but not used anymore in rectangle logic

# Missile parameters
MISSILE_SPEED = 2000.0          # faster than aircraft so it can catch you
MISSILE_KILL_RADIUS = 250.0
MISSILE_BACK_DISTANCE = 8000.0

# Flare parameters
FLARE_LIFETIME = 5.0      # seconds
FLARE_SPAWN_SPACING = 200.0
FLARE_MISSILE_RADIUS = 2500.0
FLARE_SUCCESS_PROB = 1.0  # 100% effective now


class Aircraft:
    def __init__(self):
        self.theta = 0.0
        self.q = 0.0
        self.phi = 0.0
        self.p = 0.0
        self.h = 2000.0

        # Start bottom-center of the top view, pointing up
        self.x = random.uniform(TOPVIEW_SPAN * -0.9, TOPVIEW_SPAN * 0.9)
        # self.x = 0
        self.y = -TOPVIEW_SPAN    # bottom edge
        self.psi = 0.0            # 0 rad = forward (up)

    def step(self, delta_e, delta_a, dt):
        # --- pitch dynamics with disturbances ---
        theta_dot = self.q
        q_dot = -PITCH_DAMP * self.q + PITCH_GAIN * delta_e
        q_dot += PITCH_DISTURB * random.gauss(0.0, 1.0)

        # --- roll dynamics with disturbances ---
        phi_dot = self.p
        p_dot = -ROLL_DAMP * self.p + ROLL_GAIN * delta_a
        p_dot += ROLL_DISTURB * random.gauss(0.0, 1.0)

        # integrate attitude
        self.theta += theta_dot * dt
        self.q += q_dot * dt
        self.phi += phi_dot * dt
        self.p += p_dot * dt

        # altitude
        h_dot = FORWARD_SPEED * math.sin(self.theta)
        self.h += h_dot * dt
        if self.h < 0.0:
            self.h = 0.0

        # --- heading & ground track ---
        psi_dot = -TURN_GAIN * self.phi    # sign chosen so right roll → right turn
        self.psi += psi_dot * dt

        x_dot = FORWARD_SPEED * math.sin(self.psi)
        y_dot = FORWARD_SPEED * math.cos(self.psi)
        self.x += x_dot * dt
        self.y += y_dot * dt

    def get_sensors(self):
        return {
            "theta": self.theta,
            "q": self.q,
            "phi": self.phi,
            "p": self.p,
            "h": self.h,
            "x": self.x,
            "y": self.y,
            "psi": self.psi,
        }


class FBWController:
    def __init__(self):
        self.theta_cmd = 0.0
        self.phi_cmd = 0.0

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
        self.phi_cmd = max(-max_cmd, min(max_cmd, self.phi_cmd))

    def compute_controls(self, sensors):
        theta = sensors["theta"]
        q = sensors["q"]
        phi = sensors["phi"]
        p = sensors["p"]

        e_theta = self.theta_cmd - theta
        delta_e = KP_THETA * e_theta - KD_THETA * q

        e_phi = self.phi_cmd - phi
        delta_a = KP_PHI * e_phi - KD_PHI * p

        delta_e = max(-MAX_DEFLECTION, min(MAX_DEFLECTION, delta_e))
        delta_a = max(-MAX_DEFLECTION, min(MAX_DEFLECTION, delta_a))

        return delta_e, delta_a


class Bomb:
    def __init__(self, x, y, h, psi):
        # Start where the aircraft was when B was pressed
        self.x = x
        self.y = y

        # Effective altitude: use real h, but cap so impact isn't too slow
        self.h = min(h, BOMB_DROP_ALT)

        # Bomb inherits the aircraft's horizontal velocity at release
        self.vx = FORWARD_SPEED * math.sin(psi)
        self.vy = FORWARD_SPEED * math.cos(psi)

        # Start with zero vertical speed, then accelerate downward
        self.vz = 0.0

        self.active = True
        self.impact_pos = None
        self.reported = False  # used to only report impact once

    def step(self, dt):
        if not self.active:
            return

        # Gravity acts downward
        self.vz -= G * dt

        # Update position
        self.x += self.vx * dt
        self.y += self.vy * dt
        self.h += self.vz * dt

        # Check for ground impact
        if self.h <= 0.0:
            self.h = 0.0
            self.active = False
            self.impact_pos = (self.x, self.y)


class Missile:
    def __init__(self, x, y, h, target_aircraft):
        self.x = x
        self.y = y
        self.h = h
        self.target = target_aircraft
        self.active = True

    def distance_to_target(self):
        dx = self.target.x - self.x
        dy = self.target.y - self.y
        dz = self.target.h - self.h
        return math.sqrt(dx * dx + dy * dy + dz * dz)

    def step(self, dt):
        if not self.active:
            return None

        dx = self.target.x - self.x
        dy = self.target.y - self.y
        dz = self.target.h - self.h
        dist = math.sqrt(dx * dx + dy * dy + dz * dz) + 1e-6

        ux = dx / dist
        uy = dy / dist
        uz = dz / dist

        self.x += MISSILE_SPEED * ux * dt
        self.y += MISSILE_SPEED * uy * dt
        self.h += MISSILE_SPEED * uz * dt

        if dist < MISSILE_KILL_RADIUS:
            self.active = False
            return "hit"
        return None


class Flare:
    def __init__(self, x, y, h, spawn_time):
        self.x = x
        self.y = y
        self.h = h
        self.spawn_time = spawn_time
        self.active = True

    def update(self, current_time):
        if current_time - self.spawn_time > FLARE_LIFETIME:
            self.active = False


class FBWSim:
    def __init__(self):
        pygame.init()
        self.width = 1100
        self.height = 650
        self.screen = pygame.display.set_mode((self.width, self.height))
        pygame.display.set_caption("B-2 Style FBW Mission")
        self.fbw_mode_text = "FBW: ON (attitude hold)"
        self.track_points = []
        self.clock = pygame.time.Clock()
        self.font = pygame.font.SysFont("consolas", 18)

        self.aircraft = Aircraft()
        self.controller = FBWController()
        self.running = True

        self.sim_time = 0.0

        # Bomb state (supports 2 bombs)
        self.bombs = []
        self.max_bombs = 2
        self.bombs_remaining = self.max_bombs
        self.bomb_hit_target = False
        self.bomb_miss = False
        self.bomb_result_time = None
        self.BOMB_MSG_DURATION = 8.0
        self.bomb_away_time = None
        self.BOMB_AWAY_DURATION = 3.0

        # Enemy airspace and missile
        self.entered_airspace = False
        self.next_missile_spawn_time = None  # for continuous spawns

        self.missile = None
        self.missile_launched = False
        self.missile_hit = False

        # Flares
        self.flares = []
        self.flare_message_time = None
        self.FLARE_MSG_DURATION = 3.0
        self.missile_spoofed = False
        self.missile_spoofed_time = None
        self.SPOOF_MSG_DURATION = 5.0
        self.flare_used = False  # track if player ever used flares

        # FBW mode
        self.fbw_on = True
        self.fbw_mode_text = "FBW: ON (attitude hold)"

        # Game-over state when missile hits and no flares were used
        self.game_over = False
        self.game_over_time = None
        self.GAME_OVER_DURATION = 3.0

    def in_target_corridor(self, x, y):
        # HIT if impact is inside the central 2 km x 2 km square
        return (-TARGET_HALF_SIZE <= x <= TARGET_HALF_SIZE and
                -TARGET_HALF_SIZE <= y <= TARGET_HALF_SIZE)

    def in_enemy_airspace(self, x, y):
        # Make entire 8k x 8k world box "enemy airspace"
        # so missile timer starts right away
        return (-TOPVIEW_SPAN <= x <= TOPVIEW_SPAN and
                -TOPVIEW_SPAN <= y <= TOPVIEW_SPAN)

    def draw_horizon(self, theta, phi):
        center_x = self.width // 2
        center_y = self.height // 2
        size = 300

        pitch_scale = 10.0
        y_offset = theta * pitch_scale

        half = size // 2
        pts = np.array([[-half, 0.0], [half, 0.0]])

        c = math.cos(phi)
        s = math.sin(phi)
        R = np.array([[c, -s], [s, c]])
        pts_rot = pts @ R.T

        sky_color = (80, 120, 200)
        ground_color = (100, 80, 60)
        self.screen.fill(sky_color)
        pygame.draw.rect(
            self.screen,
            ground_color,
            pygame.Rect(0, center_y + y_offset, self.width, self.height)
        )

        pts_screen = []
        for x, y in pts_rot:
            sx = int(center_x + x)
            sy = int(center_y + y + y_offset)
            pts_screen.append((sx, sy))

        pygame.draw.line(self.screen, (255, 255, 255), pts_screen[0], pts_screen[1], 3)

        pygame.draw.line(
            self.screen,
            (255, 255, 0),
            (center_x - 20, center_y),
            (center_x + 20, center_y),
            3,
        )
        pygame.draw.line(
            self.screen,
            (255, 255, 0),
            (center_x, center_y - 10),
            (center_x, center_y + 10),
            2,
        )

    def draw_ground_map(self, sensors):
        map_width = 350
        map_height = 200
        margin = 20
        x0 = margin
        y0 = self.height - map_height - margin

        # Outline box
        pygame.draw.rect(self.screen, (30, 30, 30), (x0, y0, map_width, map_height), 1)

        # World region: x,y in [-TOPVIEW_SPAN, +TOPVIEW_SPAN]
        span = TOPVIEW_SPAN

        def world_to_map(wx, wy):
            # Normalized to [-1, +1]
            nx = wx / span
            ny = wy / span

            # Clamp to keep dot inside
            nx = max(-1.0, min(1.0, nx))
            ny = max(-1.0, min(1.0, ny))

            # Map:
            #   nx = -1 -> left   , +1 -> right
            #   ny = -1 -> bottom , +1 -> top
            sx = x0 + int((nx + 1.0) * 0.5 * map_width)
            sy = y0 + map_height - int((ny + 1.0) * 0.5 * map_height)
            return sx, sy

        # --- draw target square (1k x 1k) centered at (0,0) ---
        tl = world_to_map(-TARGET_HALF_SIZE, +TARGET_HALF_SIZE)
        br = world_to_map(+TARGET_HALF_SIZE, -TARGET_HALF_SIZE)
        rect_w = br[0] - tl[0]
        rect_h = br[1] - tl[1]
        pygame.draw.rect(self.screen, (200, 50, 50), (tl[0], tl[1], rect_w, rect_h), 2)

        # --- draw aircraft ground track (tracer) ---
        if len(self.track_points) > 1:
            prev_screen = None
            for wx, wy in self.track_points:
                sx, sy = world_to_map(wx, wy)
                if prev_screen is not None:
                    pygame.draw.line(self.screen, (0, 180, 0), prev_screen, (sx, sy), 1)
                prev_screen = (sx, sy)

        # --- draw aircraft ---
        ax, ay = world_to_map(sensors["x"], sensors["y"])
        pygame.draw.circle(self.screen, (0, 255, 0), (ax, ay), 4)

        # --- draw bomb impacts (all bombs) ---
        for bomb in self.bombs:
            if bomb.impact_pos is not None:
                bx, by = world_to_map(bomb.impact_pos[0], bomb.impact_pos[1])
                pygame.draw.circle(self.screen, (255, 200, 0), (bx, by), 4)

        text = self.font.render("Top view (target + track)", True, (255, 255, 255))
        self.screen.blit(text, (x0, y0 - 20))

    def draw_rear_view(self, sensors):
        rv_width = 350
        rv_height = 200
        margin = 20
        x0 = self.width - rv_width - margin
        y0 = self.height - rv_height - margin

        pygame.draw.rect(self.screen, (10, 10, 10), (x0, y0, rv_width, rv_height))
        pygame.draw.rect(self.screen, (200, 200, 200), (x0, y0, rv_width, rv_height), 1)

        cx = x0 + rv_width // 2
        cy = y0 + rv_height - 30

        pygame.draw.polygon(
            self.screen,
            (0, 255, 0),
            [(cx, cy),
             (cx - 10, cy + 20),
             (cx + 10, cy + 20)]
        )

        label = self.font.render("Rear view (missile + flares)", True, (255, 255, 255))
        self.screen.blit(label, (x0 + 5, y0 + 5))

        # Missile
        if self.missile and (self.missile.active or self.missile_hit or self.missile_spoofed):
            dx = self.missile.x - sensors["x"]
            dy = self.missile.y - sensors["y"]
            dz = self.missile.h - sensors["h"]

            psi = sensors["psi"]
            c = math.cos(-psi)
            s = math.sin(-psi)
            fx = c * dx - s * dy
            fy = s * dx + c * dy

            rx = -fx
            ry = dz

            scale_x = 0.01
            scale_y = 0.01

            mx = cx + int(rx * scale_x)
            my = cy - int(ry * scale_y)

            mx = max(x0 + 10, min(x0 + rv_width - 10, mx))
            my = max(y0 + 20, min(y0 + rv_height - 10, my))

            if self.missile.active:
                color = (255, 50, 50)
            elif self.missile_hit:
                color = (255, 150, 150)
            else:
                color = (180, 180, 255)
            pygame.draw.circle(self.screen, color, (mx, my), 5)

        # Flares
        for flare in self.flares:
            if not flare.active:
                continue
            dx = flare.x - sensors["x"]
            dy = flare.y - sensors["y"]
            dz = flare.h - sensors["h"]

            psi = sensors["psi"]
            c = math.cos(-psi)
            s = math.sin(-psi)
            fx = c * dx - s * dy
            fy = s * dx + c * dy

            rx = -fx
            ry = dz

            scale_x = 0.01
            scale_y = 0.01

            fx_screen = cx + int(rx * scale_x)
            fy_screen = cy - int(ry * scale_y)

            fx_screen = max(x0 + 10, min(x0 + rv_width - 10, fx_screen))
            fy_screen = max(y0 + 20, min(y0 + rv_height - 10, fy_screen))

            pygame.draw.circle(self.screen, (255, 140, 0), (fx_screen, fy_screen), 4)

    def draw_hud_text(self, sensors):
        theta_deg = math.degrees(sensors["theta"])
        phi_deg = math.degrees(sensors["phi"])
        psi_deg = math.degrees(sensors["psi"])
        h = sensors["h"]
        x = sensors["x"]
        y = sensors["y"]

        dist_to_target = math.sqrt((x - TARGET_X) ** 2 + (y - TARGET_Y) ** 2)

        lines = [
            f"Yaw  : {psi_deg:6.1f} deg",
            f"Pitch: {theta_deg:6.1f} deg   Cmd: {math.degrees(self.controller.theta_cmd):6.1f} deg",
            f"Roll : {phi_deg:6.1f} deg   Cmd: {math.degrees(self.controller.phi_cmd):6.1f} deg",
            f"Alt  : {h:7.1f} m",
            f"Range to target: {dist_to_target:7.1f} m",
            self.fbw_mode_text,
            f"Bombs (B): {self.bombs_remaining}/{self.max_bombs}",
            "Flares (F): unlimited",
            "Controls: Arrows = pitch/roll, B = bomb, F = flares, C = toggle FBW, ESC = quit",
        ]

        if self.in_target_corridor(x, y) and self.bombs_remaining > 0 and not self.game_over:
            lines.append(f"Over target corridor: press B to drop bomb ({self.bombs_remaining}/{self.max_bombs})")

        if self.bomb_away_time is not None and self.sim_time - self.bomb_away_time < self.BOMB_AWAY_DURATION:
            lines.append("Bomb away")

        if self.bomb_result_time is not None and self.sim_time - self.bomb_result_time < self.BOMB_MSG_DURATION:
            if self.bomb_hit_target:
                lines.append("Bomb impact: TARGET HIT")
            elif self.bomb_miss:
                lines.append("Bomb impact: MISS")

        if self.missile_launched and self.missile and self.missile.active and not self.missile_hit:
            lines.append("Missile inbound")

        if self.missile_hit:
            lines.append("Missile impact: AIRCRAFT HIT")

        if self.flare_message_time is not None and self.sim_time - self.flare_message_time < self.FLARE_MSG_DURATION:
            lines.append("Flares deployed")

        if self.missile_spoofed and self.missile_spoofed_time is not None:
            if self.sim_time - self.missile_spoofed_time < self.SPOOF_MSG_DURATION:
                lines.append("Missile spoofed by flares")

        if self.game_over:
            lines.append("GAME OVER: Missile hit you")
            

        x0 = 20
        y0 = 20
        for line in lines:
            surf = self.font.render(line, True, (255, 255, 255))
            self.screen.blit(surf, (x0, y0))
            y0 += 22

    def spawn_missile_if_time(self, sensors):
        x = sensors["x"]
        y = sensors["y"]

        # First time entering enemy airspace starts the missile logic
        if not self.entered_airspace:
            if self.in_enemy_airspace(x, y):
                self.entered_airspace = True
                self.next_missile_spawn_time = self.sim_time + random.uniform(1.0, 4.0)
            return

        # If we're "out" of airspace (shouldn't really happen with whole box), do nothing
        if not self.in_enemy_airspace(x, y):
            return
        if self.missile_hit:
            return

        # If no next spawn time set, schedule one
        if self.next_missile_spawn_time is None:
            self.next_missile_spawn_time = self.sim_time + random.uniform(1.0, 4.0)

        # Wait until it's time and there's no active missile
        if self.sim_time >= self.next_missile_spawn_time:
            if not (self.missile and self.missile.active):
                back_dist = MISSILE_BACK_DISTANCE
                psi = sensors["psi"]

                mx = x - back_dist * math.sin(psi)
                my = y - back_dist * math.cos(psi)
                mh = max(sensors["h"], 300.0)

                self.missile = Missile(mx, my, mh, self.aircraft)
                self.missile_launched = True
                self.missile_spoofed = False
                self.missile_spoofed_time = None

                # Schedule the next missile spawn
                self.next_missile_spawn_time = self.sim_time + random.uniform(1.0, 4.0)

    def update_bomb_logic(self):
        # Step all active bombs
        for bomb in self.bombs:
            if bomb.active:
                bomb.step(DT)

        # Process impacts once
        for bomb in self.bombs:
            if (not bomb.active) and bomb.impact_pos is not None and not bomb.reported:
                bx, by = bomb.impact_pos
                if self.in_target_corridor(bx, by):
                    self.bomb_hit_target = True
                    self.bomb_miss = False
                else:
                    self.bomb_hit_target = False
                    self.bomb_miss = True
                self.bomb_result_time = self.sim_time
                bomb.reported = True

        if self.bomb_result_time is not None:
            if self.sim_time - self.bomb_result_time > self.BOMB_MSG_DURATION:
                self.bomb_result_time = None

    def update_missile_logic(self):
        if self.missile and self.missile.active:
            result = self.missile.step(DT)
            if result == "hit":
                self.missile_hit = True
                # If hit and no flares were ever used, trigger game over
                if not self.game_over:
                    self.game_over = True
                    self.game_over_time = self.sim_time

    def update_flares(self):
        # Age out flares
        for flare in self.flares:
            flare.update(self.sim_time)
        # Keep only active
        self.flares = [f for f in self.flares if f.active]

        # Check missile spoofing
        if self.missile and self.missile.active and self.flares:
            for flare in self.flares:
                dx = self.missile.x - flare.x
                dy = self.missile.y - flare.y
                dz = self.missile.h - flare.h
                dist = math.sqrt(dx * dx + dy * dy + dz * dz)
                if dist < FLARE_MISSILE_RADIUS:
                    # 100% success
                    if random.random() < FLARE_SUCCESS_PROB:
                        self.missile.active = False
                        self.missile_spoofed = True
                        self.missile_spoofed_time = self.sim_time
                    break

    def deploy_flares(self, sensors):
        # Spawn a short string of flares behind aircraft
        psi = sensors["psi"]
        back_dir_x = -math.cos(psi)
        back_dir_y = -math.sin(psi)

        for k in range(3):
            fx = sensors["x"] + back_dir_x * (100.0 + k * FLARE_SPAWN_SPACING)
            fy = sensors["y"] + back_dir_y * (100.0 + k * FLARE_SPAWN_SPACING)
            fh = sensors["h"]
            self.flares.append(Flare(fx, fy, fh, self.sim_time))

        self.flare_message_time = self.sim_time
        self.flare_used = True

    def run(self):
        while self.running:
            self.sim_time += DT

            # 1) Handle events (key presses, quit, etc.)
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
                        if self.fbw_on:
                            self.fbw_mode_text = "FBW: ON (attitude hold)"
                        else:
                            self.fbw_mode_text = "FBW: OFF (direct surfaces)"
                    if event.key == pygame.K_f:
                        sensors = self.aircraft.get_sensors()
                        self.deploy_flares(sensors)

            # 2) Poll held keys
            keys = pygame.key.get_pressed()
            if keys[pygame.K_ESCAPE]:
                self.running = False

            # If game over, just show frozen scene + message, then exit after some time
            if self.game_over:
                sensors = self.aircraft.get_sensors()
                self.draw_horizon(sensors["theta"], sensors["phi"])
                self.draw_hud_text(sensors)
                self.draw_ground_map(sensors)
                self.draw_rear_view(sensors)

                pygame.display.flip()
                self.clock.tick(FPS)

                if self.game_over_time is not None:
                    if self.sim_time - self.game_over_time > self.GAME_OVER_DURATION:
                        self.running = False
                continue

            # 3) Compute control inputs (FBW ON vs OFF)
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

            # 4) Step aircraft physics
            self.aircraft.step(delta_e, delta_a, DT)
            sensors = self.aircraft.get_sensors()

            # 5) Record position for tracer
            self.track_points.append((sensors["x"], sensors["y"]))
            if len(self.track_points) > 2000:
                self.track_points.pop(0)

            # 6) Update game logic
            self.update_bomb_logic()
            self.spawn_missile_if_time(sensors)
            self.update_missile_logic()
            self.update_flares()

            # 7) Draw everything
            self.draw_horizon(sensors["theta"], sensors["phi"])
            self.draw_hud_text(sensors)
            self.draw_ground_map(sensors)
            self.draw_rear_view(sensors)

            pygame.display.flip()
            self.clock.tick(FPS)

        pygame.quit()


if __name__ == "__main__":
    sim = FBWSim()
    sim.run()

```
