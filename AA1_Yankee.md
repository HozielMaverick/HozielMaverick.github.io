# AA-1 Yankee Electric Propulsion Conversion Project

## Overview
This project explored the full engineering workflow required to replace the AA-1 Yankee’s original Lycoming O-235-C2C gasoline engine with a modern electric propulsion system. The objective was to understand performance impacts, battery requirements, and the redesign of avionics and electrical architecture needed to safely integrate an all-electric powerplant into a certified light aircraft.

---

## My Role
As the electrical engineering lead on the team, I was responsible for:

1. **Propulsion component selection**  
   – Motor, inverter, battery modules, and contactor architecture  

2. **Performance modelling**  
   – Energy budgeting for all flight phases (takeoff, climb, cruise, descent, reserve) + analysis of flight range and fligh time.

3. **Avionics + electrical systems integration**  
   – HV/LV architecture, logic-enable paths, CAN-bus communication, safety interlocks

---

## Technical Details

### 1. Propulsion Component Selection
The first task was selecting an electric motor to replace the 115-HP (≈85 kW) Lycoming O-235-C2C engine [1].  
To ensure sufficient performance margin and avoid under-powering the aircraft, we selected the **MGM Compro BLDCin 200 kW** motor [5], which provides:

- 200 kW peak power  
- 150 kW continuous power  
- An RPM range compatible with the AA-1 propeller, eliminating the need for a gearbox  

Although this motor exceeds the power of the original Lycoming engine, excess thrust can be safely managed through **inverter power limiting**, preventing overspeed of flight control surfaces.

---

### 2. Battery Mass Allowance
The mechanical team determined that replacing the piston engine and associated systems resulted in a new empty weight of **1150 lb**, not including batteries. With a maximum takeoff weight (MTOW) of **1500 lb**, and reserving 350 lb for passengers and baggage [2], the allowable battery mass is:

m_battery = 1500 lb − 1150 lb = 350 lb ≈ 160 kg

MGM Compro advertises up to **400 Wh/kg** energy density [4], but a conservative **300 Wh/kg** was used to account for non-cell mass. This yields a pack energy of:

E_pack = 160 kg × 300 Wh/kg = 48 kWh

Only **80% usable energy** is available in aviation applications for safety, giving:

E_usable = 0.8 × 48 = 38.4 kWh

---

### 3. Performance Modeling

#### **3.1 Takeoff + Climb Energy**

Potential energy:

E_climb = (m g h) / (η × 3.6×10^6)

Using:  
m = 1500 lb = 681 kg  
h = 8000 ft  
η = 0.76  

E_climb = 5.95 kWh

Takeoff energy (70 mph):

E_to = (1/2) m V_to² = 0.09 kWh

E_climb+to = 6.04 kWh

---

#### **3.2 Reserve Energy**

E_reserve = 60 kW × (10/60) = 10 kWh

---

#### **3.3 Cruise Energy + Duration**

E_cruise = E_usable − E_climb+to − E_reserve  
E_cruise = 38.4 − 6.04 − 10 = 22.36 kWh

t_cruise = E_cruise / 60 = 0.373 h = 22 min

R_cruise = 117 kt × 0.373 = 43.5 NM

---

### 4. Climb Range and Duration

t_climb = 8000 / 710 = 11 min

![Climb Geometry](/images/aero_climb_angle_and_speed.jpg)

Angle:

α = arcsin(3.607 / 42.47) = 4.86°

Horizontal speed:

Vx_climb = 95 × cos(4.86°) = 94.66 mph = 82.26 kt

Range:

R_climb = 82.26 × 0.1833 = 15.1 NM

---

### 5. Descent Range and Duration

Approach speed range: 80–90 mph → use 85 mph midpoint.

![Descent Geometry](/images/aero_descent_angle_and_speed.jpg)

Descent angle:

β = arcsin(50 / 750) = 3.82°

Horizontal:

Vx_descent = 85 × cos(3.82°) = 84.81 mph = 73.70 kt

Vertical:

Vy_descent = 85 × sin(3.82°) = 498 ft/min

t_descent = 8000 / 498 = 16 min

Range:

R_descent = 73.70 × 0.2666 = 19.65 NM

Verification:  
Interpolated AC 700-028 value = 19.22 NM

![Descent Verification](/images/aero_descent_verif.jpg)

---

### 6. Total Mission Performance

R_total = R_climb + R_cruise + R_descent  
R_total = 15.1 + 43.5 + 19.5 = **78.1 NM**

t_total = t_climb + t_cruise + t_descent  
t_total = 11 + 22 + 16 = **49 min**

These results maintain a **10 kWh energy reserve**, giving a conservative estimate of electric mission performance.

---

## Avionics & Electrical Integration

### 1. Removal of Lycoming-Dependent Electrical Systems

![Electrical Systems Removed](/images/aero_remove_electricals.jpg)

The following systems were removed entirely:  
- Alternator and field circuits  
- Voltage regulator  
- Oil-pressure switch  
- Overvoltage protection relay  
- Starter solenoid and starter motor  
- Magneto switch  
- Fuel-pump, starter-flasher, and oil-temperature indicators  

Because these systems relied on mechanical engine operation, they were eliminated for the electric-propulsion architecture.  

---

### 2. High-Voltage (HV) Architecture

**AAMCi V2-800 Contactor Module**  
- Connects/disconnects the HV pack  
- Performs pre-charge  
- Manages main contactor logic  
- Provides HV interlocks and surge protection

**BMS-16Ai Battery Management System**  
- Monitors cell voltages, temperatures, SOC, SOH  
- Performs active balancing  
- Authorizes HV activation  
- Communicates with inverter and contactor via CAN 

**HBCi 400-800 High-Voltage Inverter**  
- Converts HV DC to 3-phase AC  
- Provides torque commands, RPM feedback, temperature data  
- Exchanges CAN messages  
- Accepts analog throttle input  

**HV Master Switch (Honeywell TL-Series)**  
- Manual arming switch  
- Interrupts only 12-V enable lines  
- Required because MGM does not supply a cockpit arming switch  

---

### 3. Reworked Electrical System Architecture

![Electric System Diagram](/images/ESD.jpg)

Circuit breakers are **not shown** for diagram clarity, but all components on both the HV and LV buses are properly fused or breaker-protected.

The system is divided into **three coordinated subsystems**:
---

### 3.1 Low-Voltage (LV) Bus
- Powered by the 12-V battery  
- Feeds avionics and controls  
- Sends throttle analog signal to inverter  
- Receives charge via DC-DC converter  

---

### 3.2 High-Voltage (HV) Bus
- Feeds inverter, contactor, propulsion motor  
- Powers DC-DC converter  
- Internal interlocks enforce safe activation order  

---

### 3.3 Logic & Control Layer (CAN Bus)
- Connects BMS, inverter, and avionics  
- Transmits SOC/SOH, temps, RPM, torque, faults  
- Uses twisted-pair wiring for EMI reduction  

---

### 4. Cockpit Display Redesign

To maintain familiar pilot workload and ergonomics, the cockpit was redesigned to preserve original gauge functionality where possible. The original panel from [3] was marked with:

- **Red** → instruments removed entirely
- **Yellow** → instruments replaced with electric equivalents

![Cockpit Display Revisions](/images/aero_CD.jpg)  
![Component Conversion List](/images/aero_list.jpg)

Examples: - The tachometer becomes a digital RPM indicator driven from the inverter. - Oil temperature and pressure gauges are replaced with SOC/SOH from the BMS and motor temperature from the inverter.

---

### 5. Reassigned Controls (Starred Components)
Some cockpit controls were repurposed to serve new functions under electric propulsion:

**Suction Gauge → AHRS**  
Uses a standalone AHRS module (GSU 25, MGL AHRS-1, etc.) to provide attitude data.

**Fuel Tank Selector → Manual Coolant Override**  

![Manual Coolant Override](/images/aero_MCO.jpg)

**Carb Heat Control → Aviation-Grade USB-C Charger**  
The Carb Heat Control becomes unused in an electric aircraft. Its panel location is repurposed for a **shielded, aviation-grade USB-C power outlet**. This provides clean, RF-filtered power for portable avionics such as EFB tablets, handheld GPS receivers, or communication devices, improving cockpit usability without increasing pilot workload.

---

### 6. Throttle Potentiometer System

To preserve the pilot’s familiar throttle-push control, the mechanical throttle knob now drives a **potentiometer-based voltage divider**. This produces an analog voltage:

V_out = 5 V × (R_top / (R_top + R_bottom))

This signal is read by the inverter’s analog-input channel and mapped to torque or power demand. 

![Throttle Potentiometer System](/images/aero_TC.jpg)

This architecture also enables software-level **power limiting**, ensuring the electric motor cannot exceed speeds that would exceed the aircraft’s structural control-surface limits.

---

## References
[1] Grumman Owners & Pilots Association, AA-1/AA-1A/AA-1B/AA-1C POH Data.  
[2] AA-1A Pilot Operating Handbook.  
[3] AA-1C Pilot Operating Handbook.  
[4] MGM Compro Battery Systems, https://mgm-compro.com/ess-batteries/  
[5] MGM Compro BLDCin 200 kW Motor, https://mgm-compro.com/electric-motor/mgm-bldcin-200kw/  
[6] MGM Compro Propulsion Systems Overview, https://mgm-compro.com/eps-propulsion-systems/  
[7] MGM Compro HBCi 400/800 Inverter, https://mgm-compro.com/controllers-inverters-escs/hbci-400800/  
[8] MGM Compro AAMCi Series Contactor, https://mgm-compro.com/electronics-contactors/aamci-series-v2-800/  
[9] Climb Rate Data, Manufacturer Performance Tables.  
[10] Hepperle, M., Electric Flight, DLR.  
[11] Transport Canada Advisory Circular AC 700-028.  
