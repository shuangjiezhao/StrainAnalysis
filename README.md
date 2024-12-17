# Strain Analysis Toolkit

## Overview
The Strain Analysis Toolkit is a powerful Python package designed for in-depth analysis and visualization of structural deformations in 2D materials, such as graphene and Covalent Organic Frameworks (COFs). It provides a comprehensive set of tools to quantify and understand atomic-scale strain, enabling researchers to gain valuable insights into the mechanical behavior of these materials.

## Key Features
- Atomic layer separation analysis
- Interlayer distance calculation
- Detailed strain mapping in both in-plane and out-of-plane directions
- Graphene and COF structure manipulation
- Customizable visualization options

## Installation

### Prerequisites
- Python 3.8 or higher
- NumPy
- Matplotlib
- ASE (Atomic Simulation Environment)

### Install via pip
```bash
pip install git+https://github.com/shuangjiezhao/strain-analysis.git
```
## Usage
```python
### Example
pythonCopyfrom strain_analysis import StrainAnalysis

### Initialize strain analysis
sa = StrainAnalysis()

### Separate layers from original and optimized structures
layers_before = sa.layer_sep(original_structure, True)
layers_after = sa.layer_sep(optimized_structure, True)

# Generate flakes for strain analysis
flake_before = sa.flake_generator(layers_before[0], 2, 2, 1)
flake_after = sa.flake_generator(layers_after[0], 2, 2, 1)

# Prepare data for strain analysis
sa.prepare_G(flake_before, flake_after, reference_distance)

# Visualize out-of-plane strain
fig, ax = sa.out_plane(color=['blue', 'red'])
plt.show()

# Visualize in-plane strain
fig, ax = sa.in_plane(color=['green', 'yellow'])
plt.show()

# Add COF frame on top of graphene
sa.add_COF(fig, ax, 'green', 'cof_frame', 'your/path/')
```

## Contributing
We welcome contributions to this project! If you have any ideas, bug fixes, or feature requests, please feel free to submit a pull request. Before starting, please review the contributing guidelines.

## License
This project is licensed under the MIT License.
