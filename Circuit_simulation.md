# Circuit Simulation Framework – SPICE Style MNA Solver

This project is a SPICE style circuit simulation engine built across four major assignments. It parses netlists, builds Modified Nodal Analysis (MNA) matrices, and solves circuits in DC, transient, and frequency domain settings. Over the four stages it progresses from linear DC operating point, to nonlinear diode operating point with Newton Raphson, to transient analysis with Backward Euler, and finally to a Harmonic Balance solver that can be compared directly to Backward Euler steady state results.

---

## Features Overview

- Netlist based circuit description (R, C, L, I, V, diodes)  
- Automated Modified Nodal Analysis (MNA) matrix construction  
- Linear DC operating point solving  
- Nonlinear DC operating point with diodes using Newton Raphson  
- Transient simulation with Backward Euler (.tran)  
- Harmonic Balance solver for periodic steady state  
- Direct comparison of HB results with BE steady state waveforms  
- MATLAB implementation with modular functions for parsing, stamping, and solving  

---

# 1. Linear DC Operating Point Solver (.op, Dev1)

The first assignment implements a DC operating point solver for purely linear circuits.

Key aspects:

- Parses a simple SPICE style netlist  
- Builds the MNA matrices for resistors, independent current sources, and voltage sources  
- Solves the linear system for node voltages and selected branch currents  
- Provides the base infrastructure for later nonlinear and time domain work  

This stage focuses on getting the MNA formulation correct and matching DC results to known solutions or reference tools.

**DC Netlist**
![DC Node Voltages Netlist](/images/dc_netlist.jpg)

**DC MNA Matrices**
![DC MNA Matrices](/images/dc_mna.jpg)

---

# 2. Nonlinear DC Operating Point with Diodes (.op, Dev2)

The second assignment extends the DC solver to handle nonlinear devices, specifically diodes.

New capabilities:

- Adds diode elements to the netlist format  
- Formulates the nonlinear MNA equations that include the diode I–V relationship  
- Implements Newton Raphson iteration at the DC operating point  
- Computes the Jacobian and residual for the diode currents and voltages  
- Achieves convergence to a consistent DC operating point even with nonlinear elements  

At this stage the simulator becomes a nonlinear DC analysis tool that can handle rectifiers, clipping networks, and other diode based circuits.

**DC Non Linear Netlist**
![DC non Linear Netlist](/images/dc_nl_netlist.jpg)


**DC Non Linear MNA Matrices**
![DC non Linear MNA Matrices](/images/dc_nl_mna.jpg)

---

# 3. Transient Analysis with Backward Euler (.tran, Dev3)

The third assignment introduces time domain simulation and the .tran style analysis.

Main features:

- Implements Backward Euler as an implicit integration method  
- Converts capacitors and inductors into companion models at each timestep  
- Rebuilds and solves the MNA system on every time step  
- Supports circuits that include both linear elements and diodes  
- Produces time domain waveforms for node voltages and branch currents  

This stage turns the framework into a dynamic circuit simulator that can model charging and discharging, diode switching behavior, and other time varying phenomena.

**Transient Analysis Netlist**
![Transient Circuit](/images/trancirc.jpg)

**Transiant Analysis Output**
![Transient Output](/images/tranoutput.jpg)

---

# 4. Harmonic Balance vs Backward Euler (Dev4)

The fourth assignment focuses on periodic steady state analysis and compares two approaches:

- Time domain Backward Euler  
- Frequency domain Harmonic Balance (HB)

Additions in this stage:

- Builds Fourier transform matrices (Gamma and inverse Gamma)  
- Constructs the diagonal admittance matrix Y(jkω₀) for linear components in the frequency domain  
- Forms the nonlinear algebraic system for HB and solves it with Newton Raphson in the frequency domain  
- Recovers time domain periodic waveforms by inverse transform  
- Compares HB steady state results with the long time Backward Euler solution for the same circuit  

This shows how a frequency domain method like Harmonic Balance can reach the periodic steady state solution directly, while Backward Euler approaches the same steady state through time stepping.

**Steady State Analysis Circuit, f=1kHz, Is=1e-14, Vt=25mv**
![Circuit](/images/hbcirc.jpg)

**HB VS BE Output**
![Harmonic Balance vs Backward Euleur Results](/images/hbbe.jpg)

---

# Architectural Summary

### Netlist Parsing

- Reads SPICE style netlist files  
- Maps component names and node labels to internal indices  
- Stores element data in structured arrays suitable for stamping  

### MNA Assembly

- Stamps resistors, sources, capacitors, inductors, and diodes into the MNA system  
- Handles augmented rows and columns for voltage sources  
- Builds the correct system for DC, transient, or HB depending on the analysis type  

### Solver Layer

- Linear solve for DC operating point (Dev1)  
- Newton Raphson for nonlinear DC operating point with diodes (Dev2)  
- Time stepping loop with Backward Euler for .tran simulations (Dev3)  
- Frequency domain Newton Raphson for Harmonic Balance steady state (Dev4)  

### Outputs

- DC node voltages and branch currents  
- Transient waveforms in the time domain  
- Periodic steady state waveforms from HB  
- Harmonic content and comparison between HB and BE  

---

# Skills Demonstrated

- Circuit analysis and simulation using MNA  
- Nonlinear numerical methods such as Newton Raphson  
- Implicit time integration with Backward Euler  
- Harmonic Balance formulation in the frequency domain  
- MATLAB programming for solvers, matrix assembly, and plotting  
- Systematic extension of a codebase from linear DC to nonlinear and then to transient and HB analysis  

---

# Access to Full Code

[**Dev1**](Dev1.md)

[**Dev2**](Dev2.md)

[**Dev3**](Dev3.md)

[**Dev4**](Dev4.md)
