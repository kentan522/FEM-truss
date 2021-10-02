# FEM-truss
Finite element analysis on the collapse mechanism of a truss using object-oriented programming on MATLAB. 
Based on Imperial College London's MEng Civil Engineering 2nd Year Computational Methods 2020/2021 coursework. 

The TRUSS.m file represents the truss class that can be used to generate an instance of a truss object.

The TRUSS_script.m file represents an example script that can be used to conduct the analysis on a defined truss object.

## Explaination of methods of TRUSS.m file

Some syntax:

𝜆<sub>tot</sub> = load proportionality factor

𝐹<sub>tot</sub> = vector of total forces for each member

𝑈<sub>tot</sub> = vector of total displacement for each member

| Methods | Description |
| ------------- | ------------- |
| Class constructor  | To initialise the truss object based on its given parameters |
| Assembly | To assemble the global stiffness matrix for the global truss structure  |
| Solver| To partition and solve the matrix, obtaining displacements at each dofs as well as the reaction forces at the restrained dofs  |
| Post processor| To obtain the 𝜆<sub>tot</sub> , 𝐹<sub>tot</sub> , 𝑈<sub>tot</sub> as well as the critical member at which failure will occur locally |
| Plotting | To plot the deformation state after each event of a member failure, as well as the axial forces within each member |
| Updater | To replace each failed member by its critical squash/buckling load resolved in the direction of its nodal dofs, stored in the global nodal force vector  |
| Matrix singularity check | To check for the stability of the KFF sub-matrix before attempting to solve the global stiffness matrix  |

## Explaination of parts within TRUSS_script.m file
Part 1: Formulation of truss geometry by specifying nodal coordinates, dofs and elemental connection.

Part 2: Materially Non-Linear Analysis (failure by plastic collapse only).

Part 3: Geometrically Non-Linear Analysis (failure by plastic collapse or elastic buckling).

Part 4: Varying the radii of each members to conduct analysis on truss geometries with members of varying radii.
