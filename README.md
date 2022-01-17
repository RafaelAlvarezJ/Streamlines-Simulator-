# Streamlines-Simulator

In Petroleum Engineering, specifically Reservoir Engineering, a simulator is a valuable computational tool for forecasting fluid production in a reservoir. It requires input information such as rock and fluid properties.

Simulators are often built using finite differences or finite element method; they return great results, but computational cost sometimes makes them really slow. A streamlines simulator has a big advantage over them because all transport calculations are IMPES (which are simple, and require no solver like GMRES), therefore simulations are a lot faster.

This project is a streamlines simulator I've been developing. It considers two wells: an injector and a producer, located in opposites corners of the domain, constant injection and production rate, and a homogeneus and isotropic reservoir (constant porosity and permeability), more input data can be found in "DatosPollock". It is not finished yet cause I want to optimize it using numba, but it's already giving interesting and fast results!
