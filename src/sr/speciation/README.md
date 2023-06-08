# SRangesBirthDeathModel

The `SRangesBirthDeathModel` class is a variant of the fossilized birth-death model under budding (asymmetric) speciation with stratigraphic ranges.

## Class Description

The `SRangesBirthDeathModel` class extends the `SABirthDeathModel` class and represents a birth-death model with stratigraphic ranges.

### Methods

- `q(double t, double c1, double c2)`: Calculates the q-value for a given time, c1, and c2.
- `log_q(double t, double c1, double c2)`: Calculates the logarithm of the q-value for a given time, c1, and c2.
- `q_tilde(double t, double c1, double c2)`: Calculates the q-tilde value for a given time, c1, and c2.
- `log_q_tilde(double t, double c1, double c2)`: Calculates the logarithm of the q-tilde value for a given time, c1, and c2.
- `log_lambda_times_int_limits_p(double tOld, double tYoung, double c1, double c2)`: Calculates the logarithm of lambda times the integral limits p for given times tOld, tYoung, c1, and c2.
- `findAncestralRangeLastNode(Node node)`: Finds the last node of the ancestral range for a given node.
- `calculateLogP()`: Calculates the log probability for the birth-death model with stratigraphic ranges.