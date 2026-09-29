# Model and validation notes

These notes describe the implemented educational model in [`src/main.js`](src/main.js). The equations below are source descriptions, not a claim that this implementation agrees with a laboratory apparatus. Material values are nominal constants in the app. **All displayed transit times, forces, fields, and currents are simulation estimates.** There is no independent physical benchmark or uncertainty estimate in this repository yet.

The Vite production build passed during this update, and the interactive flows were reviewed in a browser. Those checks establish that the app can load and run; they do not check the numerical physics. No automated physics invariant, convergence study, or comparison with physical measurements has been completed for this update. The checklist below remains proposed work.

## Geometry and shared assumptions

The magnet is a uniformly magnetized cylinder, 12 mm in diameter and 25 mm long. Its selectable material changes the modeled remanence and density; the default N42 setting uses 1.30 T and 7,500 kg/m³. Its field calculation uses the ideal-cylinder solution implemented in `magnetField`, with Bulirsch's complete elliptic integral routine. That field is exact for the *assumed uniform-cylinder model*, not for every real magnet.

The coil has a 30 mm bore, 0.5 mm wire, and a selectable 50–400 turns distributed across an approximately 12 mm long winding. The source samples the turns across the winding and calculates the magnet-to-coil flux linkage by reciprocity using a loop vector potential and numerical quadrature (`fluxLinkage`, `coilLoops`). For roughly equal turn fluxes, $\Lambda\approx N\Phi$; more generally $\Lambda=\sum_i\Phi_i$. A position table is interpolated for the live derivative. The induced electromotive force approximates $\mathcal E=-d\Lambda/dt$, including translation and modeled rotation.

Wire resistance is $R_{\rm coil}=\rho\ell/A$, using the selected wire's fixed resistivity. Air-core inductance uses the empirical multilayer-coil expression printed in Wheeler's 1928 paper, with the radius, length, and winding depth converted to inches and the result converted from microhenries to henries. That paper attributes this multilayer expression to L. A. Hazeltine. A 1 Ω meter resistance is added in series. The circuit step solves

$$
L\,\frac{dI}{dt}=V_{\rm battery}+\mathcal E-(R_{\rm coil}+1\,\Omega)I.
$$

The meter and main current readout show total coil current. The induced-current trace uses $I-V_{\rm battery}/(R_{\rm coil}+1\,\Omega)$ as its displayed component. The magnet's electromagnetic force uses the flux-linkage position derivative and coil current. Swing mode drives the magnet; dragging uses a spring-like pointer target, so neither is a free-motion laboratory experiment.

The optional iron core is a 27 mm diameter, 40 mm long rod. It is modeled as a uniformly magnetized, laminated material with an effective susceptibility derived from an approximate demagnetizing factor and capped at 2 T. The plotted linkage includes the core's magnet-driven contribution under a battery-off reference condition. When the battery is on, current changes the core response, so Fig. 3 is a reference curve rather than a full plot of the operating core state. Core hysteresis, detailed domain behavior, and core eddy currents are omitted. Coil winding resistance and the tabulated tube resistivities do not vary with temperature. The model does not include winding capacitance, frequency-dependent skin effects, manufacturing tolerances, or external magnetic sources.

## Tube drop

The tube is 30 cm long with a 14 mm inside diameter and adjustable 0.5–3 mm wall thickness. Its wall is represented by **120 coaxial rings**, each 2.5 mm tall (`buildTube`). Each conductive ring has resistance based on the selected fixed resistivity and an estimated self-inductance; pairs of rings are linked through a coaxial-loop mutual-inductance calculation. The magnet's flux through each ring comes from the same ideal-cylinder magnet model. The update solves the coupled ring-current and vertical magnet-motion equations with gravity, electromagnetic drag, and the modeled inductances. Plastic has infinite resistance in the model, so it carries no eddy current and has no magnetic drag.

The magnet moves only along the tube axis. The tube is treated as nonmagnetic and uniform, with no seams, temperature change, radial motion, tilt, air resistance, or impact rebound. The terminal-speed readout is a **quasistatic estimate at mid-tube** from the model's local flux slopes and resistances; it is not a measured terminal velocity. Playback uses a fixed 1/60 s physics tick and eight 1/480 s tube substeps per tick. The tube matrix is assembled for that substep. Resistive heating is calculated as $\sum_k I_k^2R_k$ and does not feed back into material properties.

Transit timing begins when the magnet's leading edge crosses the tube's top plane and ends when its trailing edge clears the bottom plane. The displayed duration therefore depends on the magnet length as well as the 30 cm tube length. A saved A/B trial is a record of two model runs under their selected settings, not a comparison with measured data.

## What the field drawing includes

In **tube mode**, the lines show the magnet's field alone and the strength map is hidden. The tube's colored bands and arrows separately encode the computed ring currents. The source deliberately omits the rings' field from the field-line drawing.

In **coil mode**, the app carries pretraced magnet-only lines when the coil contribution is small. Above the source's current/core thresholds it traces a combined field from the magnet, coil current, and optional modeled core magnetization (`Bworld`, `updateFieldLines`). Thus a magnet-only trace at low current is a visualization simplification, while the combined trace still inherits the core and coil approximations above. The colors encode $|B|$ on a logarithmic scale. Line spacing and count are chosen for legibility; they are not a separate field measurement.

## Sources for methods

- Norman Derby and Stanislaw Olbert, “[Cylindrical Magnets and Ideal Solenoids](https://arxiv.org/abs/0909.3880),” *American Journal of Physics* **78**, 229 (2010), [doi:10.1119/1.3256157](https://doi.org/10.1119/1.3256157). This is the source identified in `magnetField` for the ideal-cylinder field and `cel` evaluation. Their paper also presents its own tube-drop simulation and measured comparison; those results have **not** been reproduced as a validation of this app.
- Harold A. Wheeler, “[Simple Inductance Formulas for Radio Coils](https://doi.org/10.1109/JRPROC.1928.221309),” *Proceedings of the IRE* **16**, 1398–1400 (1928). The multilayer expression used by `buildCoil` appears as Eq. (1) in this paper and is attributed there to L. A. Hazeltine. It is an empirical air-core approximation, not a calibrated inductance for this app's winding or iron-core option.

## Physics benchmark checklist — proposed, not yet run

These checks should be automated before describing the numbers as quantitatively validated. Define tolerances and record the app version, geometry, and time step with every result.

1. **Magnet field:** compare the calculated on-axis field with the closed-form uniformly magnetized cylinder expression; check radial field symmetry and convergence near the cylinder ends against an independent implementation.
2. **Flux and Faraday sign:** for an aligned coil, check position symmetry of $\Lambda$, odd symmetry of its position derivative, zero induced EMF when all geometry is stationary, and polarity reversal when the magnet is flipped. Test rotation separately from translation.
3. **Circuit:** with a fixed magnet and an air core, compare battery-only current after a voltage step against the analytic series-$RL$ response. Check the displayed total current against the induced-component trace with battery voltage both zero and nonzero.
4. **Plastic control:** verify all ring currents, magnetic drag, and resistive heat are zero within numerical tolerance, and compare the timed fall with a gravity-only calculation using the same entry and exit planes.
5. **Conductive energy accounting:** over an unforced drop, compare gravitational energy lost with kinetic energy gained, $\int\sum_k I_k^2R_k\,dt$, and change in magnetic energy. State any omitted terms before setting a pass tolerance.
6. **Convergence:** repeat representative coil and tube cases with finer flux quadrature, smaller ring height, and smaller simulation steps. Confirm the tube solve's matrix uses the same substep length as the motion integrator. Report changes in peak current, drag, and transit time.
7. **Physical comparison:** if measured transit times or field maps are available, use independently measured magnet moment/remanence, dimensions, wall thickness, conductivity at measured temperature, and timing positions. Compare across more than one tube material and state uncertainty rather than fitting one case and calling it validation.
