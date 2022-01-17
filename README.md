# Streamlines-Simulator

In Petroleum Engineering, specifically in Reservoir Engineering, a simulator is a valuable computational tool for forecasting fluid production in a reservoir. It requires input information such as rock and fluid properties, and also total simulation time, injection/production rate and so on.

Simulators are often built using finite differences or finite element method; they return great results, but computational cost sometimes makes them really slow. A streamlines simulator has a big advantage over them because all transport calculations are IMPES (which are simple, and require no solver like GMRES), therefore simulations are a lot faster.

This project is a streamlines simulator I've been developing. It is not finished yet but it's already giving interesting and fast results!
