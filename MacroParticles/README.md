# MacroParticles: semi-solid VdW-like particle simulation

This project implements an interactive Van-der-Waals-like particle physics simulation. It allows to create new particles, interact with them and visualize the strain between them. Additionally you can enable the rendering of the dynamic plots of total kinetic energy and average/max strain over time.

![An example of a particle system][img/img1.png]
![Strain graph and energy/strain plots][img/img2.png]

Graphs:

[red]: log of total kinetic energy
[green]: average relative strain
[orange]: max relative strain

### Controls

- Left click with `Ctrl` pressed creates a new particle. Particles can be placed on top of each other -- the interaction potential was designed to handle this situation
- Pressing `Shift` creates a temporary particle. It follow the cursor as long `Shift` is pressed and interacts with the rest of the particles as if it was just a normal particle
- `Backspace` removes the last added particle and `Del` deletes them all
- Pressing `F` freezes/unfreezes the simulation

### Building

The project uses a simple graphics library for Windows (found in ../_Graphics). To build the project you can either use the included Code::Blocks project file (*.cbp) or edit the compiler paths in `build.bat` script and build the project using the provided `Makefile`.
