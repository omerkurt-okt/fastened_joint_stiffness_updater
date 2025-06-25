# Fastened Joint Stiffness Updater

A comprehensive finite element analysis utility for calculating and updating CBUSH element stiffness values using the Huth formula for fastened joints in Nastran FE models.

![Python](https://img.shields.io/badge/python-3.7+-blue.svg)
![License](https://img.shields.io/badge/license-Apache%202.0-green.svg)
![PyQt5](https://img.shields.io/badge/GUI-PyQt5-orange.svg)

## Overview

The CBUSH Huth Stiffness Updater is an advanced engineering tool designed for structural analysts working with fastened joints in finite element models. It automates the calculation of fastener stiffness values based on the empirical Huth formula, which accounts for the flexibility of both the fastener and connected plates.

### Key Capabilities

- **Huth Formula Implementation**: Accurate calculation of fastener stiffness based on plate properties and fastener characteristics
- **Material Orientation Support**: Advanced handling of composite materials with proper coordinate transformations
- **Multi-Plate Connection Analysis**: Automatic identification and analysis of complex multi-plate joints
- **Airbus Method Support**: Implementation of Airbus-specific fastener modeling approach
- **Parallel Processing**: Multi-core processing for large models with thousands of fasteners
- **Classical Lamination Theory**: Full CLT implementation for composite plates (PCOMP, PCOMPG, PCOMPP)
- **Fastener Database**: Comprehensive fastener specification management
- **SET Management**: Import/export of fastener groups from Nastran SET cards

## Features

### Technical Features

- Support for multiple Nastran element types (CBUSH, CQUAD4, CTRIA3, CQUAD8, CTRIA6)
- Material property extraction (MAT1, MAT8, PSHELL, PCOMP, PCOMPG, PCOMPP)
- RBE3 connection handling for CBUSH-to-plate connections
- Coordinate system transformations for material orientation
- Multiple connection analysis methods (Normal, HyperMesh, Two-Plate Enforcement)
- Mix mode support for existing Airbus/Normal connections

### User Interface

- Intuitive PyQt5-based graphical interface
- Multiple tabs for different functionalities:
  - Calculation tab for main processing
  - Fastener Groups management
  - Fastener Specifications database
  - Comprehensive help documentation
- Real-time progress tracking and logging
- Import/Export capabilities for configurations

## Usage
### Basic Workflow

- Load Model: Select your Nastran BDF file containing CBUSH elements
- Import Fastener Sets (Optional): Load fastener group definitions from SET cards
- Configure Parameters: Set default fastener properties and calculation methods
- Process Model: Run the calculation with desired number of parallel workers
- Export Results: Save the updated BDF file with new stiffness values

### Fastener Group Management
The tool supports automatic detection of fastener groups from SET naming conventions:

Pattern: DESCRIPTION1__DESCRIPTION2__DIAMETER__SPECIFICATION__ROWNUMBER__CPNUMBER__ENUM

Example: "Frame200__Duct__6.35__ELS438__2__2__0"

### Connection Processing Modes

- Convert to Airbus Method: Creates additional CBUSH spanning outer plates with only axial stiffness
- Convert to Normal Method: Removes Airbus CBUSH elements and updates regular CBUSHes
- Mix Mode: Preserves existing connection types (default)

### Material Orientation Options

- Directional Modulus: Uses transformed Ex for K2, Ey for K3 based on material orientation
- Effective Modulus: Uses geometric mean E_eff = √(E1×E2) for both K2 and K3

## Technical Details
### Huth Formula Implementation
The tool implements the empirical Huth formula for shear stiffness:

K_shear = 1 / C_total

Where:
C_total = [(t₁ + t₂)/(2d)]^a × [b₁×C₁ + b₂×C₂]
C₁ = 1/(t₁×E₁) + 1/(2×t₁×E_f)
C₂ = 1/(t₂×E₂) + 1/(2×t₂×E_f)

### Material Parameters
![image](https://github.com/user-attachments/assets/7eddd8f8-af41-44f7-beb0-c4444acfcc68)

### Coordinate Systems
The tool handles multiple coordinate systems:
- Element local coordinate systems
- Material coordinate systems (theta_mcid)
- Fastener coordinate systems

## Acknowledgments
- pyNastran library for BDF file handling
- PyQt5 for the graphical user interface
- TF-X Structural Analysis Team for validation and feedback
