# Induction bench

An interactive browser app for exploring **electromagnetic induction**. Move a cylindrical magnet relative to a coil to see flux linkage, induced voltage, current, force, and magnetic field. Drop the same magnet through a tube to compare eddy current braking in conducting materials with a plastic control.

The app covers these two induction experiments. It is not a survey of all electromagnetism. Its numbers are **simulation estimates**, not measurements or independently validated predictions; see [Model and validation](MODEL.md).

## Public demo

The GitHub Pages demo is paused. To use the app, follow the local instructions below.

## Run locally

From this directory:

```sh
npm ci
npm run dev
```

Open the local address printed by Vite. To check a production build, run `npm run build`; `npm run preview` serves the built files locally. The visualization uses WebGL through Three.js, so use a browser with WebGL enabled. No account or API key is needed for the simulation.

## Learn, then explore

The five-step guided path starts with a still magnet, then asks you to move it and observe the sign of the induced voltage. It next explores speed or coil turns, reverses the magnet's poles, and compares conducting and insulating tubes. You can hide the guide and use the full bench controls at any point.

In the coil experiment, drag the magnet or use the labeled position slider. A focused range slider works with the keyboard's arrow keys. Adjust the swing, coil turns, angle, wire, and magnet material to see which changes affect the response. The **Battery** control adds a separate source of current; its Electromagnet experiment explains how a battery drives coil poles. When it is on, distinguish total coil current from the induced component shown on the time plot. The **Soft iron** option is a simplified core model.

In **Tube drop**, choose a tube material and wall thickness, then drop the magnet. Use pause, step, slower playback, and replay to inspect the brief fall. The trial cards retain two results for comparison, including the material, thickness, and modeled transit time. Switch the plot between speed and magnetic drag; both traces align by progress through the tube. Select **Review A** or **Review B**, then move the Review position slider to inspect the saved ring-current bands and readouts. Plastic is the nonconducting control in this model.

## Reading the display

- **Flux linkage** $\Lambda$ is the sum of magnetic flux through the coil turns, approximately $N\Phi$ when the turns have similar flux. Changes in the magnet-driven linkage induce voltage: $\mathcal E\approx-d\Lambda/dt$. Motion, magnet rotation, or coil rotation can change it. A field can be present while the magnet is still and the induced voltage is zero. With both battery and iron core on, Fig. 3 shows a battery-off reference curve because the core response depends on current.
- **Coil mode field lines** show the modeled magnet field at low coil current. When the coil current or core contribution becomes large enough, the app traces the modeled combined magnet, coil, and optional core field. The strength map uses the same scope. Field lines are selected traces, not individual particles or a count of flux quanta.
- **Tube mode field lines** show the magnet's modeled field only. The strength map is hidden in this mode. Colored bands and arrows on the tube represent the separately calculated induced ring currents. The drawn lines do not include the field created by those eddy currents.
- **Current:** with a battery connected, the main coil current includes the battery and induction response; the oscilloscope's induced-current channel separates the response relative to the steady battery baseline. In tube mode, the meter samples one ring below the magnet, while the other readouts summarize the largest ring current, magnetic drag, and resistive heating.
- **Transit time** starts when the leading end of the magnet enters the tube and stops when its trailing end clears the bottom. It is a result of the model and selected geometry, not a timing measurement from a physical tube.

## Model and validation

[MODEL.md](MODEL.md) records the implemented geometry, equations, approximations, display scope, references, and a benchmark checklist. The paper references identify methods used in the source; they do not establish that this particular implementation or its numerical outputs have been independently validated.
