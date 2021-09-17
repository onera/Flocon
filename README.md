# Flocon

Flocon is an open-source benchmark library aimed at gathering flow control problems.

These benchmarks are built to be relevant for various topics such as model reduction/identification, actuators/sensors placement, control synthesis and its implementation.
For that purpose, each benchmark can be modified (e.g. number, position and nature of actuators/sensors) to highlight some specific point.

The library is implemented in python but the core of each benchmark may rely on another language, e.g. for simulation. Specific requirements are indicated in the associated README files.

## Available benchmarks

The list of benchmarks will be completed progressively.

|   Name     | Description                                        |  Simulator |
|-----------:|:---------------------------------------------------|:----------:|
|   `step`   | The system is an amplifier and the control goal is disturbance rejection. See [this page](step).|  FreeFem++ |
|   `ks`     | Non-parallel Kuramoto-Sivashinsky equation. See [this page](ks).                                |  FreeFem++ |

## Authors

The Flocon library is developped by Pierre Vuillemin (pierre.vuillemin@onera.fr) and Charles Poussot-Vassal (charles.poussot-vassal@onera.fr) within an internal ONERA project.
