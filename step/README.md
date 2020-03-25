# Step use-case

This use case comes from the article [1] and is schematised in the Figure below.

![Default configuration](./static/step_default.png)

## Requirements

This core simulator of this benchmark is implemented [FreeFem++](https://freefem.org/).

## Getting started

The benchmark can be loaded as follows,

```python
import Flocon
# Load the step benchmark
step = Flocon.get_bench('step')
# By default, the configuration is the same as in [1]
step.plot_config()
```
Note that this creates a (sqlite) database in which simulations results (in particular baseflow) are stored.



## Adjusting the congifuration

### Actuators
### Sensors

## References

[1] Hervé, A., Sipp, D., Schmid, P. J., & Samuelides, M. (2012). A physics-based approach to flow control using system identification. Journal of Fluid Mechanics, 702, 26-58.
