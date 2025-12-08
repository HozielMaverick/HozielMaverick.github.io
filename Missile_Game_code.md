```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
import random

# ============================================================================#
# Constants that do not change between engagements
# ============================================================================#
tmax = 75
dt = 0.001
animation_interval = 5  # milliseconds
kill_dist = 1.0

# Axes limits
AX_MIN = 0.0
AX_MAX_XY = 50000.0  # 50 km by 50 km in XY
AX_MAX_Z = 50000.0   # allow Z to go up to 50 km

# Enemy airspace limits for escape condition (XY only)
AIRSPACE_MIN = 0.0
AIRSPACE_MAX = 50000.0  # escaping only when X or Y leaves [0, 50000]

# Spawn regions
XY_MIN = 0.0
XY_MAX = 50000.0
Z_TARGET_MIN = 25000.0
Z_TARGET_MAX = 35000.0

# Base position range (10 km to 40 km)
BASE_MIN_XY = 10000.0
BASE_MAX_XY = 40000.0

# Bomb properties
BOMB_COLOR = "yellow"  # bomb color

# Missile launch time is based on radar detection
missile_launch_time = None

# How long we keep showing the result before starting a new engagement (in animation steps)
post_intercept_frames = 1000

# Score counters (persist across engagements)
missile_neutralized_target = 0
aircraft_escaped = 0
payload_delivered = 0  # TOTAL successful deliveries across completed engagements

# Aircraft type definitions
# (speeds + curve_min/curve_max from your old code, plus max_turn_rate)
AIRCRAFT_TYPES = {
    "F22": {
        "color": "cyan",
        "flare_trigger_distance": 100.0,
        "flare_success_prob": 0.50,
        "speed_min": 2300,
        "speed_max": 2400,
        "curve_min": 3,
        "curve_max": 8,
        "max_turn_rate": 1.0,   # rad/s
        "radar_detection_prob": 0.05,
    },
    "A10": {
        "color": "orange",
        "flare_trigger_distance": 100.0,
        "flare_success_prob": 0.70,
        "speed_min": 2000,
        "speed_max": 2100,
        "curve_min": 5,
        "curve_max": 10,
        "max_turn_rate": 0.7,
        "radar_detection_prob": 0.50,
    },
    "B2": {
        "color": "blue",
        "flare_trigger_distance": 100.0,
        "flare_success_prob": 0.50,
        "speed_min": 1700,
        "speed_max": 1800,
        "curve_min": 7,
        "curve_max": 12,
        "max_turn_rate": 0.4,
        "radar_detection_prob": 0.01,
    },
    "X15": {
        "color": "magenta",
        "flare_trigger_distance": 100.0,
        "flare_success_prob": 0.05,
        "speed_min": 3200,
        "speed_max": 3300,
        "curve_min": 10,
        "curve_max": 20,
        "max_turn_rate": 0.10,
        "radar_detection_prob": 0.50,
    },
}

# Globals that will be reinitialized for each engagement
targ_vel = None
miss_vel = None
times = None
n_points = None
target_states = None
missile_states = None
bomb_states = None
target_heading = None

missile_launched = False
missile_frozen = False  # missile stops when flares succeed

intercept_time = None
intercept_index = None  # index of first time we enter kill_dist
escape_index = None     # index of first time aircraft leaves airspace (XY only)
crash_index = None      # index of first time aircraft goes below 500 m
outcome_index = None    # frame index where outcome occurs
outcome_type = None     # "missile", "aircraft", "crash", or None

# Start locations (re-randomized each engagement)
missile_start_loc = None
aircraft_start_loc = None

# Radar and flare state
current_aircraft_name = None
current_aircraft_params = None
radar_detected = False
last_radar_check_second = -1
radar_detection_index = None  # frame index where radar detects the aircraft

flares_deployed = False
missile_tracking_flare = False
flare_position = None
flare_effective_index = None  # time index when flares succeed
flare_deploy_index = None     # time index when flares are deployed (success or fail)

# Bomb / base / explosion state
base_position = None
bomb_dropped = False
bomb_active = False
bomb_drop_index = None
explosion_index = None
payload_recorded = False          # per-engagement flag: payload hit target
payload_counted_in_score = False  # per-engagement flag: already added to payload_delivered

# Post-bomb OLD-LOGIC flight profile
climb_rate_curve = None
post_Straight_time = None
post_curve_time1 = None
post_curve_time2 = None
post_Straight_time2 = None
post_turn_angle1 = None
post_turn_angle2 = None
post_targ_vel = None
post_profile_initialized = False
post_drop_alt = None
post_psi = None  # heading angle used in old logic

# Spawn/escape sides (spawn still used, escape not needed any more for heading)
spawn_side = None
escape_side = None
escape_point = None  # not used for post-bomb heading now

# Animation control
current_index = 0
frame_step = 1  # will be updated after we know n_points

# We will reuse these scatter handles across engagements
start_aircraft_scatter = None
start_missile_scatter = None
base_marker = None

# Legend handle so we can update it dynamically
legend_handle = None


# ============================================================================#
# Helper: rotate vector "current" toward "desired" by at most max_angle
# ============================================================================#
def rotate_towards(current, desired, max_angle):
    cur = current / (np.linalg.norm(current) + 1e-9)
    des = desired / (np.linalg.norm(desired) + 1e-9)

    dot = np.clip(cur[0]*des[0] + cur[1]*des[1], -1.0, 1.0)
    angle = np.arccos(dot)

    if angle < 1e-6:
        return des.copy()

    if angle <= max_angle:
        return des.copy()

    s = max_angle / angle
    new_vec = (1.0 - s) * cur + s * des
    new_vec[2] = 0.0
    new_vec /= (np.linalg.norm(new_vec) + 1e-9)
    return new_vec


# ============================================================================#
# CREATE 3D PLOT (fixed axes)
# ============================================================================#
fig = plt.figure(figsize=(14, 10))
ax = fig.add_subplot(111, projection="3d")

ax.set_xlim(AX_MIN, AX_MAX_XY)
ax.set_ylim(AX_MIN, AX_MAX_XY)
ax.set_zlim(AX_MIN, AX_MAX_Z)
ax.set_box_aspect([AX_MAX_XY, AX_MAX_XY, AX_MAX_Z])

ax.set_xlabel("X (m)")
ax.set_ylabel("Y (m)")
ax.set_zlabel("Z (m)")
ax.set_title("3D Missile-Aircraft Pursuit Simulation - Trace And Chase")
ax.grid(True)
ax.view_init(elev=20, azim=45)

# Preallocate arrays
times = np.arange(0, tmax, dt)
n_points = len(times)

target_states = np.zeros((n_points, 3))
missile_states = np.zeros((n_points, 3))
bomb_states = np.full((n_points, 3), np.nan)

# Artists
target_point, = ax.plot([], [], [], "o", markersize=10, label="Aircraft")
target_trail, = ax.plot([], [], [], "-", linewidth=2, alpha=0.5)

missile_point, = ax.plot([], [], [], "ro", markersize=8, label="Heat Seeking Missile")
missile_trail, = ax.plot([], [], [], "r-", linewidth=1.5, alpha=0.5, label="Missile Trail")

bomb_point, = ax.plot([], [], [], "o", markersize=8, color=BOMB_COLOR, label="Bomb")
bomb_trail, = ax.plot([], [], [], "-", linewidth=1.5, color=BOMB_COLOR, alpha=0.7, label="Bomb Trail")

base_marker = ax.scatter([], [], [], marker="s", s=80, color="brown", label="Enemy Base")

time_text = ax.text2D(0.02, 0.95, "", transform=ax.transAxes, fontsize=12)
speed_text = ax.text2D(0.02, 0.90, "", transform=ax.transAxes, fontsize=10)
distance_text = ax.text2D(0.02, 0.85, "", transform=ax.transAxes, fontsize=10)
aircraft_text = ax.text2D(0.02, 0.80, "", transform=ax.transAxes, fontsize=10)

score_text = ax.text2D(0.90, 0.02, "", transform=ax.transAxes, fontsize=12)

# Radar indicator bottom-left
radar_dot = ax.text2D(0.02, 0.05, "●", transform=ax.transAxes, fontsize=18, color="green")
radar_label = ax.text2D(0.05, 0.05, "RADAR LOCK", transform=ax.transAxes, fontsize=12, color="black")

# Markers
intercept_marker, = ax.plot([], [], [], "kX", markersize=10, label="Intercept")
crash_marker, = ax.plot([], [], [], "kP", markersize=10, label="Crash")
flare_marker, = ax.plot([], [], [], marker="*", linestyle="None", color="green",
                        markersize=15, label="Successful Flares, Aircraft is Safe")
explosion_marker, = ax.plot([], [], [], marker="o", linestyle="None", markersize=20,
                            color="gold", markeredgecolor="red", label="Explosion")


# ============================================================================#
# FUNCTION TO GENERATE A NEW ENGAGEMENT
# ============================================================================#
def generate_new_engagement():
    global targ_vel, miss_vel
    global target_states, missile_states, bomb_states
    global missile_launched, missile_frozen, intercept_time, intercept_index, escape_index
    global crash_index, outcome_index, outcome_type
    global current_index, frame_step
    global missile_start_loc, aircraft_start_loc
    global start_aircraft_scatter, start_missile_scatter
    global missile_launch_time
    global current_aircraft_name, current_aircraft_params
    global radar_detected, last_radar_check_second, radar_detection_index
    global flares_deployed, missile_tracking_flare, flare_position
    global flare_effective_index, flare_deploy_index
    global target_point, target_trail
    global legend_handle
    global flare_marker, radar_dot
    global intercept_marker, crash_marker
    global base_position, base_marker
    global bomb_dropped, bomb_active, bomb_drop_index, explosion_index
    global payload_recorded, payload_counted_in_score
    global target_heading
    global spawn_side, escape_side, escape_point
    global climb_rate_curve, post_Straight_time, post_curve_time1, post_curve_time2
    global post_Straight_time2, post_turn_angle1, post_turn_angle2
    global post_targ_vel, post_profile_initialized, post_drop_alt, post_psi

    # Reset radar UI
    radar_dot.set_color("green")

    # Pick aircraft type
    type_of_aircraft = random.randint(1, 4)
    if type_of_aircraft == 1:
        current_aircraft_name = "B2"
    elif type_of_aircraft == 2:
        current_aircraft_name = "F22"
    elif type_of_aircraft == 3:
        current_aircraft_name = "A10"
    else:
        current_aircraft_name = "X15"

    current_aircraft_params = AIRCRAFT_TYPES[current_aircraft_name]

    # Reset radar/flare/missile state
    radar_detected = False
    radar_detection_index = None
    last_radar_check_second = -1
    missile_launch_time = None

    flares_deployed = False
    missile_tracking_flare = False
    flare_position = None
    flare_effective_index = None
    flare_deploy_index = None

    missile_launched = False
    missile_frozen = False

    intercept_index = None
    intercept_time = None
    escape_index = None
    crash_index = None
    outcome_index = None
    outcome_type = None

    # Reset bomb & explosion state for THIS engagement
    bomb_states[:] = np.nan
    bomb_dropped = False
    bomb_active = False
    bomb_drop_index = None
    explosion_index = None
    payload_recorded = False          # did we hit near the base?
    payload_counted_in_score = False  # has that hit been added to global tally?

    # Reset post-bomb profile
    climb_rate_curve = random.choice([-0.001, 0.001])
    post_Straight_time = None
    post_curve_time1 = None
    post_curve_time2 = None
    post_Straight_time2 = None
    post_turn_angle1 = None
    post_turn_angle2 = None
    post_targ_vel = None
    post_profile_initialized = False
    post_drop_alt = None
    post_psi = None

    # Clear markers
    flare_marker.set_data([], [])
    flare_marker.set_3d_properties([])
    crash_marker.set_data([], [])
    crash_marker.set_3d_properties([])
    intercept_marker.set_data([], [])
    intercept_marker.set_3d_properties([])
    explosion_marker.set_data([], [])
    explosion_marker.set_3d_properties([])

    # Random missile start location on ground
    missile_start_loc = np.array([
        random.uniform(XY_MIN, XY_MAX),
        random.uniform(XY_MIN, XY_MAX),
        0.0,
    ])

    # Choose random spawn side for aircraft
    spawn_side = random.choice(["left", "right", "bottom", "top"])
    if spawn_side == "left":
        ax_x = XY_MIN
        ax_y = random.uniform(XY_MIN, XY_MAX)
    elif spawn_side == "right":
        ax_x = XY_MAX
        ax_y = random.uniform(XY_MIN, XY_MAX)
    elif spawn_side == "bottom":
        ax_x = random.uniform(XY_MIN, XY_MAX)
        ax_y = XY_MIN
    else:  # top
        ax_x = random.uniform(XY_MIN, XY_MAX)
        ax_y = XY_MAX

    aircraft_start_loc = np.array([
        ax_x,
        ax_y,
        random.uniform(Z_TARGET_MIN, Z_TARGET_MAX),
    ])

    # Random base position in [10 km, 40 km]
    base_position = np.array([
        random.uniform(BASE_MIN_XY, BASE_MAX_XY),
        random.uniform(BASE_MIN_XY, BASE_MAX_XY),
        0.0,
    ])

    # Update base marker
    base_marker._offsets3d = (
        np.array([base_position[0]]),
        np.array([base_position[1]]),
        np.array([base_position[2]]),
    )

    # Randomize speeds (ingress speed)
    speed_min = current_aircraft_params["speed_min"]
    speed_max = current_aircraft_params["speed_max"]
    targ_vel = random.randint(speed_min, speed_max)
    miss_vel = random.randint(3700, 3800)

    # Initialize states
    target_states[:] = 0.0
    missile_states[:] = 0.0
    target_states[0] = aircraft_start_loc
    missile_states[0] = missile_start_loc

    # Initial heading: toward base
    to_base = np.array([
        base_position[0] - aircraft_start_loc[0],
        base_position[1] - aircraft_start_loc[1],
        0.0,
    ])
    norm_tb = np.linalg.norm(to_base[:2])
    if norm_tb > 0:
        target_heading = to_base / norm_tb
    else:
        target_heading = np.array([1.0, 0.0, 0.0])

    max_turn_rate = current_aircraft_params["max_turn_rate"]

    print("====================================================")
    print(f"New engagement generated with {n_points} points.")
    print(f"Aircraft type = {current_aircraft_name}")
    print(f"Spawn side = {spawn_side}")
    print(f"Target speed (ingress) = {targ_vel} m/s")
    print(f"Missile speed = {miss_vel} m/s")
    print(f"Target start = {aircraft_start_loc}")
    print(f"Missile start = {missile_start_loc}")
    print(f"Enemy base at = {base_position}")

    # Generate trajectories
    for i in range(1, n_points):
        t = times[i]

        # ==================== AIRCRAFT DYNAMICS =====================

        if bomb_dropped and post_profile_initialized and (bomb_drop_index is not None) and (i > bomb_drop_index):
            # ---------------- POST-BOMB: OLD LOGIC (FULL TRAJECTORY) ----------------
            t_rel = (i - bomb_drop_index) * dt

            # Time boundaries relative to bomb drop
            t1 = post_Straight_time
            t2 = t1 + post_curve_time1
            t3 = t2 + post_Straight_time2
            t4 = t3 + post_curve_time2

            # Heading rate psi_dot
            if t_rel <= t1:
                psi_dot = 0.0
            elif t_rel <= t2:
                psi_dot = post_turn_angle1 / post_curve_time1
            elif t_rel <= t3:
                psi_dot = 0.0
            elif t_rel <= t4:
                psi_dot = post_turn_angle2 / post_curve_time2
            else:
                psi_dot = 0.0

            post_psi += psi_dot * dt

            dx = np.cos(post_psi)
            dy = np.sin(post_psi)

            # XY motion with post-bomb speed
            target_states[i, 0] = target_states[i - 1, 0] + post_targ_vel * dt * dx
            target_states[i, 1] = target_states[i - 1, 1] + post_targ_vel * dt * dy

            # Vertical modulation as in old code (two cosine bumps)
            # First turn contribution
            if t_rel <= t1:
                tau1 = 0.0
            elif t_rel <= t2:
                tau1 = t_rel - t1
            else:
                tau1 = post_curve_time1

            # Second turn contribution
            if t_rel <= t3:
                tau2 = 0.0
            elif t_rel <= t4:
                tau2 = t_rel - t3
            else:
                tau2 = post_curve_time2

            z_offset1 = climb_rate_curve * (post_targ_vel ** 2) * (
                1.0 - np.cos(np.pi * tau1 / post_curve_time1)
            )
            z_offset2 = climb_rate_curve * (post_targ_vel ** 2) * (
                1.0 - np.cos(np.pi * tau2 / post_curve_time2)
            )

            target_states[i, 2] = post_drop_alt + z_offset1 + z_offset2

            # Clamp altitude
            if target_states[i, 2] < 0.0:
                target_states[i, 2] = 0.0
            elif target_states[i, 2] > AX_MAX_Z:
                target_states[i, 2] = AX_MAX_Z

        else:
            # ---------------- PRE-BOMB: FLY TOWARD BASE (YOUR CURRENT LOGIC) -------
            vec = np.array([
                base_position[0] - target_states[i - 1, 0],
                base_position[1] - target_states[i - 1, 1],
                0.0,
            ])

            norm_vec = np.linalg.norm(vec[:2])
            if norm_vec > 0:
                desired_dir = vec / norm_vec
            else:
                desired_dir = target_heading.copy()

            max_angle_step = max_turn_rate * dt
            target_heading = rotate_towards(target_heading, desired_dir, max_angle_step)

            target_states[i, 0] = target_states[i - 1, 0] + targ_vel * dt * target_heading[0]
            target_states[i, 1] = target_states[i - 1, 1] + targ_vel * dt * target_heading[1]
            target_states[i, 2] = target_states[i - 1, 2]  # constant altitude pre-bomb

        # ======================= BOMB DYNAMICS ================================
        if not bomb_dropped:
            # Check proximity to base to drop bomb
            if (
                abs(target_states[i, 0] - base_position[0]) <= 50.0
                and abs(target_states[i, 1] - base_position[1]) <= 50.0
            ):
                bomb_dropped = True
                bomb_active = True
                bomb_drop_index = i
                bomb_states[i] = target_states[i]
                print(
                    f"{current_aircraft_name} drops bomb at t = {t:.2f}s, "
                    f"pos = {target_states[i]}"
                )

                # ------------ INITIALIZE POST-BOMB OLD-LOGIC PROFILE -------------
                post_drop_alt = target_states[i, 2]

                # Timing based on aircraft agility
                Straight_time = random.randint(1, 5)
                curve_min = current_aircraft_params["curve_min"]
                curve_max = current_aircraft_params["curve_max"]
                curve_time1 = random.randint(curve_min, curve_max)
                mult = random.choice([0.5, 1.1])
                curve_time2 = curve_time1 * mult
                Straight_time2 = random.randint(1, 5)

                post_Straight_time = Straight_time
                post_curve_time1 = curve_time1
                post_curve_time2 = curve_time2
                post_Straight_time2 = Straight_time2

                # Post-bomb speed
                speed_min = current_aircraft_params["speed_min"]
                speed_max = current_aircraft_params["speed_max"]
                post_targ_vel = random.randint(speed_min, speed_max)

                # Turn angles like old code
                sign = random.choice([-1, 1])
                post_turn_angle1 = sign * -np.pi * 4 / 3
                sign2 = random.uniform(-1, 1)
                while abs(sign2) < 0.4:
                    sign2 = random.uniform(-1, 1)
                post_turn_angle2 = sign2 * post_turn_angle1

                # Initial heading angle psi from current heading
                post_psi = np.arctan2(target_heading[1], target_heading[0])

                post_profile_initialized = True

                print(f"POST-BOMB PROFILE:")
                print(f"  Straight_time = {Straight_time:.1f} s")
                print(f"  Curve_time1   = {curve_time1:.1f} s")
                print(f"  Straight_time2= {Straight_time2:.1f} s")
                print(f"  Curve_time2   = {curve_time2:.1f} s")
                print(f"  post_targ_vel = {post_targ_vel} m/s")
                print(f"  climb_rate_curve = {climb_rate_curve}")
        else:
            # Bomb already dropped: propagate bomb vertically
            if bomb_active:
                bomb_states[i, 0] = bomb_states[i - 1, 0]
                bomb_states[i, 1] = bomb_states[i - 1, 1]
                bomb_states[i, 2] = bomb_states[i - 1, 2] - 2.0 * miss_vel * dt

                if bomb_states[i, 2] <= 0.0:
                    bomb_states[i, 2] = 0.0
                    bomb_active = False
                    explosion_index = i
                    print(
                        f"Bomb hits ground at t = {t:.2f}s, "
                        f"pos = {bomb_states[i]}"
                    )
                    # Only mark this engagement as having delivered payload.
                    dist_to_base = np.linalg.norm(
                        bomb_states[i, :2] - base_position[:2]
                    )
                    if (
                        bomb_dropped
                        and (bomb_drop_index is not None)
                        and (i >= bomb_drop_index)
                        and (dist_to_base <= 200.0)
                        and (not payload_recorded)
                    ):
                        payload_recorded = True
                        print("Payload successfully delivered on target!")
            else:
                bomb_states[i] = bomb_states[i - 1]

        if not bomb_dropped and np.isnan(bomb_states[i, 0]):
            bomb_states[i] = np.array([np.nan, np.nan, np.nan])

        # ====================== RADAR + MISSILE DYNAMICS ======================
        if not radar_detected:
            current_second = int(t)
            if current_second != last_radar_check_second and current_second >= 1:
                last_radar_check_second = current_second
                p_det = current_aircraft_params["radar_detection_prob"]
                if random.random() < p_det:
                    radar_detected = True
                    radar_detection_index = i
                    missile_launch_time = t
                    missile_launched = True
                    print(
                        f"Radar detected {current_aircraft_name} "
                        f"at t = {t:.2f}s, missile launched"
                    )

        if not missile_launched:
            missile_states[i] = missile_states[i - 1]
            continue

        if missile_frozen:
            missile_states[i] = missile_states[i - 1]
            continue

        if missile_tracking_flare and flare_position is not None:
            guidance_target = flare_position
        else:
            guidance_target = target_states[i]

        direction = guidance_target - missile_states[i - 1]
        distance_to_guidance_target = np.linalg.norm(direction)

        if distance_to_guidance_target > 0:
            unitvec = direction / distance_to_guidance_target
            dr = unitvec * miss_vel * dt
            missile_states[i] = missile_states[i - 1] + dr
        else:
            missile_states[i] = missile_states[i - 1]

        # Flares
        if (not flares_deployed) and missile_launched:
            distance_to_aircraft = np.linalg.norm(
                target_states[i] - missile_states[i]
            )
            if distance_to_aircraft <= current_aircraft_params["flare_trigger_distance"]:
                flares_deployed = True
                flare_deploy_index = i
                print(
                    f"{current_aircraft_name} deploys flares at t = {t:.2f}s, "
                    f"distance = {distance_to_aircraft:.2f} m"
                )

                if random.random() < current_aircraft_params["flare_success_prob"]:
                    missile_tracking_flare = True
                    missile_frozen = True
                    flare_position = missile_states[i].copy()
                    flare_effective_index = i
                    print("Flares succeeded, missile is now spoofed and frozen in place")
                else:
                    missile_tracking_flare = False
                    flare_position = None
                    flare_effective_index = None
                    print("Flares failed, missile continues tracking aircraft")

    # ====================== POST-PROCESS OUTCOMES ============================
    distances = np.linalg.norm(target_states - missile_states, axis=1)
    hit_indices = np.where(distances <= kill_dist)[0]

    if hit_indices.size > 0:
        intercept_index = int(hit_indices[0])
        intercept_time = times[intercept_index]
        print(
            f"Hit: first entry into kill_dist at t = {intercept_time:.2f}s, "
            f"d = {distances[intercept_index]:.2f} m"
        )
    else:
        intercept_index = None
        intercept_time = None

    # Escape (any side of airspace)
    escape_mask = (
        (target_states[:, 0] < AIRSPACE_MIN)
        | (target_states[:, 0] > AIRSPACE_MAX)
        | (target_states[:, 1] < AIRSPACE_MIN)
        | (target_states[:, 1] > AIRSPACE_MAX)
    )
    escape_indices = np.where(escape_mask)[0]
    if escape_indices.size > 0:
        escape_index = int(escape_indices[0])
        print(
            "Escape: aircraft leaves XY airspace at t = "
            f"{times[escape_index]:.2f}s, "
            f"pos = {target_states[escape_index]}"
        )
    else:
        escape_index = None

    # Crash (Z < 500 m)
    crash_mask = target_states[:, 2] < 500.0
    crash_indices = np.where(crash_mask)[0]
    if crash_indices.size > 0:
        crash_index = int(crash_indices[0])
        print(
            "Crash: aircraft altitude below 500 m at t = "
            f"{times[crash_index]:.2f}s, "
            f"pos = {target_states[crash_index]}"
        )
    else:
        crash_index = None

    # Outcome logic – write directly to globals outcome_index / outcome_type
    outcome_index = None
    outcome_type = None

    has_flares_success = missile_tracking_flare and (flare_effective_index is not None)
    has_bomb = bomb_dropped and (bomb_drop_index is not None)
    has_explosion = has_bomb and (explosion_index is not None)

    if has_flares_success:
        candidates = []
        if escape_index is not None:
            candidates.append((escape_index, "aircraft"))
        if crash_index is not None:
            candidates.append((crash_index, "crash"))
        if candidates:
            outcome_index, outcome_type = min(candidates, key=lambda x: x[0])

    else:
        if not has_bomb:
            candidates = []
            if intercept_index is not None:
                candidates.append((intercept_index, "missile"))
            if crash_index is not None:
                candidates.append((crash_index, "crash"))
            if escape_index is not None:
                candidates.append((escape_index, "aircraft"))
            if candidates:
                outcome_index, outcome_type = min(candidates, key=lambda x: x[0])
        else:
            if intercept_index is None:
                candidates = []
                if escape_index is not None:
                    candidates.append((escape_index, "aircraft"))
                if crash_index is not None:
                    candidates.append((crash_index, "crash"))
                if candidates:
                    outcome_index, outcome_type = min(candidates, key=lambda x: x[0])
            else:
                if (bomb_drop_index is None) or (intercept_index <= bomb_drop_index):
                    candidates = [(intercept_index, "missile")]
                    if crash_index is not None:
                        candidates.append((crash_index, "crash"))
                    if escape_index is not None:
                        candidates.append((escape_index, "aircraft"))
                    outcome_index, outcome_type = min(candidates, key=lambda x: x[0])
                else:
                    if has_explosion and (intercept_index < explosion_index) and payload_recorded:
                        outcome_index = explosion_index
                        outcome_type = "missile"
                    else:
                        candidates = [(intercept_index, "missile")]
                        if crash_index is not None:
                            candidates.append((crash_index, "crash"))
                        if escape_index is not None:
                            candidates.append((escape_index, "aircraft"))
                        outcome_index, outcome_type = min(candidates, key=lambda x: x[0])

    if outcome_type == "missile":
        print("Outcome for this engagement: MISSILE NEUTRALIZED TARGET")
    elif outcome_type == "aircraft":
        print("Outcome for this engagement: AIRCRAFT ESCAPED")
    elif outcome_type == "crash":
        print("Outcome for this engagement: AIRCRAFT CRASHED (MISSILE WIN)")
    else:
        print("Outcome for this engagement: NONE")

    # Frame step
    global frame_step
    frame_step = max(1, n_points // 500)
    current_index = 0

    # Update start markers and legend
    aircraft_color = current_aircraft_params["color"]
    target_point.set_color(aircraft_color)
    target_trail.set_color(aircraft_color)
    target_point.set_label(current_aircraft_name)

    global start_aircraft_scatter, start_missile_scatter, legend_handle

    if start_aircraft_scatter is None:
        start_aircraft_scatter = ax.scatter(
            target_states[0, 0],
            target_states[0, 1],
            target_states[0, 2],
            c=aircraft_color,
            s=100,
            marker="s",
            label="Aircraft Start",
        )
    else:
        start_aircraft_scatter._offsets3d = (
            np.array([target_states[0, 0]]),
            np.array([target_states[0, 1]]),
            np.array([target_states[0, 2]]),
        )
        start_aircraft_scatter.set_color(aircraft_color)

    if start_missile_scatter is None:
        start_missile_scatter = ax.scatter(
            missile_start_loc[0],
            missile_start_loc[1],
            missile_start_loc[2],
            c="black",
            s=100,
            marker="^",
            label="Missile Start",
        )
    else:
        start_missile_scatter._offsets3d = (
            np.array([missile_start_loc[0]]),
            np.array([missile_start_loc[1]]),
            np.array([missile_start_loc[2]]),
        )

    if legend_handle is not None:
        legend_handle.remove()
    legend_handle = ax.legend()


# Initialize first engagement
generate_new_engagement()


# ============================================================================#
# ANIMATION FUNCTIONS
# ============================================================================#
def init():
    target_point.set_data([], [])
    target_point.set_3d_properties([])

    target_trail.set_data([], [])
    target_trail.set_3d_properties([])

    missile_point.set_data([], [])
    missile_point.set_3d_properties([])

    missile_trail.set_data([], [])
    missile_trail.set_3d_properties([])

    bomb_point.set_data([], [])
    bomb_point.set_3d_properties([])

    bomb_trail.set_data([], [])
    bomb_trail.set_3d_properties([])

    time_text.set_text("")
    speed_text.set_text("")
    distance_text.set_text("")
    score_text.set_text("")
    aircraft_text.set_text("")

    intercept_marker.set_data([], [])
    intercept_marker.set_3d_properties([])

    crash_marker.set_data([], [])
    crash_marker.set_3d_properties([])

    flare_marker.set_data([], [])
    flare_marker.set_3d_properties([])

    explosion_marker.set_data([], [])
    explosion_marker.set_3d_properties([])

    return (
        target_point,
        target_trail,
        missile_point,
        missile_trail,
        bomb_point,
        bomb_trail,
        time_text,
        speed_text,
        distance_text,
        aircraft_text,
        score_text,
        intercept_marker,
        crash_marker,
        flare_marker,
        explosion_marker,
    )


def update(frame):
    global current_index
    global missile_neutralized_target, aircraft_escaped, payload_delivered
    global payload_recorded, payload_counted_in_score

    # Radar color
    if radar_detection_index is not None and current_index >= radar_detection_index:
        radar_dot.set_color("red")
    else:
        if (current_index // 5) % 2 == 0:
            radar_dot.set_color("green")
        else:
            radar_dot.set_color("lime")

    # If bomb successfully delivered and explosion time reached,
    # increment payload_delivered exactly once at that moment.
    if (
        payload_recorded
        and not payload_counted_in_score
        and explosion_index is not None
        and current_index >= explosion_index
    ):
        payload_delivered += 1
        payload_counted_in_score = True
        print(f"Payload counter incremented to {payload_delivered} at t = {times[explosion_index]:.2f}s")

    # Engagement reset
    if outcome_index is not None:
        if current_index >= min(n_points - 1, outcome_index + post_intercept_frames):
            # Update scores for *finished* engagement
            if outcome_type in ("missile", "crash"):
                missile_neutralized_target += 1
            elif outcome_type == "aircraft":
                aircraft_escaped += 1

            print(
                f"Scores -> Missile: {missile_neutralized_target}, "
                f"Aircraft: {aircraft_escaped}, "
                f"Payloads: {payload_delivered}"
            )

            generate_new_engagement()
    else:
        if current_index >= n_points - 1:
            generate_new_engagement()

    i = current_index
    if outcome_index is not None and i >= outcome_index:
        idx = outcome_index
    else:
        idx = i

    # Aircraft
    target_point.set_data([target_states[idx, 0]], [target_states[idx, 1]])
    target_point.set_3d_properties([target_states[idx, 2]])

    target_trail.set_data(target_states[: idx + 1, 0], target_states[: idx + 1, 1])
    target_trail.set_3d_properties(target_states[: idx + 1, 2])

    # Missile
    missile_point.set_data([missile_states[idx, 0]], [missile_states[idx, 1]])
    missile_point.set_3d_properties([missile_states[idx, 2]])

    missile_trail.set_data(missile_states[: idx + 1, 0], missile_states[: idx + 1, 1])
    missile_trail.set_3d_properties(missile_states[: idx + 1, 2])

    # Bomb + tracer
    if not np.isnan(bomb_states[idx, 0]):
        bomb_point.set_data([bomb_states[idx, 0]], [bomb_states[idx, 1]])
        bomb_point.set_3d_properties([bomb_states[idx, 2]])

        indices = np.arange(n_points)
        valid_bomb = (~np.isnan(bomb_states[:, 0])) & (indices < idx)
        bomb_trail.set_data(bomb_states[valid_bomb, 0], bomb_states[valid_bomb, 1])
        bomb_trail.set_3d_properties(bomb_states[valid_bomb, 2])
    else:
        bomb_point.set_data([], [])
        bomb_point.set_3d_properties([])
        bomb_trail.set_data([], [])
        bomb_trail.set_3d_properties([])

    # Explosion marker
    if explosion_index is not None and idx >= explosion_index:
        explosion_marker.set_data(
            [bomb_states[explosion_index, 0]],
            [bomb_states[explosion_index, 1]],
        )
        explosion_marker.set_3d_properties([0.0])
    else:
        explosion_marker.set_data([], [])
        explosion_marker.set_3d_properties([])

    # Speed / distance
    if idx > 0:
        dx = target_states[idx, 0] - target_states[idx - 1, 0]
        dy = target_states[idx, 1] - target_states[idx - 1, 1]
        dz = target_states[idx, 2] - target_states[idx - 1, 2]
        speed = np.sqrt(dx**2 + dy**2 + dz**2) / dt
    else:
        speed = targ_vel

    distance = np.linalg.norm(target_states[idx] - missile_states[idx])

    time_text.set_text(f"Time = {times[idx]:.2f} s")
    speed_text.set_text(f"Target Speed = {speed:.1f} m/s")
    distance_text.set_text(f"Distance = {distance:.1f} m")
    aircraft_text.set_text(f"Aircraft: {current_aircraft_name}")

    # Scoreboard shows totals over completed past engagements
    score_text.set_text(
        f"Missile Hit Target: {missile_neutralized_target}\n"
        f"Aircraft Escaped: {aircraft_escaped}\n"
        f"Payloads Delivered: {payload_delivered}"
    )

    # Intercept marker
    if intercept_index is not None and idx >= intercept_index:
        intercept_marker.set_data(
            [target_states[intercept_index, 0]],
            [target_states[intercept_index, 1]],
        )
        intercept_marker.set_3d_properties([target_states[intercept_index, 2]])
    else:
        intercept_marker.set_data([], [])
        intercept_marker.set_3d_properties([])

    # Crash marker
    if outcome_type == "crash" and crash_index is not None and idx >= crash_index:
        crash_marker.set_data(
            [target_states[crash_index, 0]],
            [target_states[crash_index, 1]],
        )
        crash_marker.set_3d_properties([target_states[crash_index, 2]])
    else:
        crash_marker.set_data([], [])
        crash_marker.set_3d_properties([])

    # Flare marker
    if (
        missile_tracking_flare
        and flare_position is not None
        and flare_effective_index is not None
        and idx >= flare_effective_index
    ):
        flare_marker.set_data([flare_position[0]], [flare_position[1]])
        flare_marker.set_3d_properties([flare_position[2]])
    else:
        flare_marker.set_data([], [])
        flare_marker.set_3d_properties([])

    current_index = min(i + frame_step, n_points - 1)

    return (
        target_point,
        target_trail,
        missile_point,
        missile_trail,
        bomb_point,
        bomb_trail,
        time_text,
        speed_text,
        distance_text,
        aircraft_text,
        score_text,
        intercept_marker,
        crash_marker,
        flare_marker,
        explosion_marker,
    )


# ============================================================================#
# RUN ANIMATION
# ============================================================================#
anim = FuncAnimation(
    fig,
    update,
    init_func=init,
    blit=False,
    interval=animation_interval,
)

print("Showing animation...")
plt.show()

```
