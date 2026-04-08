Architecture and Numerical Coupling of the Purkinje-Myocardium System
1. Code Overview and Interface Design
The current implementation utilizes a highly modular, interface-driven C++ architecture built on top of OpenFOAM. It cleanly separates the definition of the physiological domains (1D Purkinje network, 3D myocardium, ECG) from the physics of how they exchange charge, and from the numerical scheme used to advance them in time.

The core abstract endpoints, tissueCouplingEndpoint and networkCouplingEndpoint, allow the domains to remain agnostic to one another. They strictly expose necessary fields (mesh, voltages, source terms, terminal nodes) and accept feedback currents.

2. Current Implementation of the Domains
The Purkinje Network (purkinjeNetworkModel): Acts as a 1D graph-backed domain. It effectively manages its internal state, reading network topology and Purkinje-Ventricular Junction (PVJ) mappings. It calculates the global volume of PVJ spheres across all parallel processors to ensure grid-independent coupling. It correctly subtracts applied coupling currents from its internal source to account for charge leaving the network.

The Myocardium (MyocardiumDomain): Acts as the primary 3D continuous domain. It exposes a volumetric sourceField where external coupling currents are deposited before its continuous PDE (e.g., Monodomain or Bidomain) is solved.

The PVJ Mapper (PVJMapper): Encapsulates the parallel and geometric complexities. It locates the 3D tissue cells associated with 1D terminal nodes, manages parallel reductions when gathering the 3D voltage (V 
m
​	
 ), and converts discrete 1D currents into continuous volumetric 3D sources using the pre-calculated sphere volumes to avoid numerical singularities.

3. The Coupler and Directionality
The physics of the connection are governed by the PVJResistanceCoupler. It evaluates the standard monodomain-monodomain resistive coupling term at the PVJs:

I 
PM
​	
 = 
R 
pvj
​	
 
V 
p
​	
 −V 
m
​	
 
​	
 
Memory Management: The coupler maintains internal buffers (terminalCurrentBuffer_, terminalSourceBuffer_) to prevent continuous memory reallocation during the tight time-stepping loop.

Coupling Directionality: The coupler explicitly handles the CouplingMode.

Unidirectional (networkToMyocardium): The coupler computes the forward current to drive the tissue but forces the network-side current buffers to 0.0. This isolates the 1D solver from the 3D tissue, allowing the Purkinje network to act as a strict Dirichlet-like driver.

Bidirectional: The computed difference is applied to both domains, allowing for physiological retrograde conduction (e.g., ectopic beats in the tissue propagating back up the Purkinje network).

4. Current Integration Schema
The system time integration is managed by staggeredElectrophysicsAdvanceScheme, which implements a Weakly Coupled Segregated Scheme (an explicit staggered/Gauss-Seidel update):

Prepare Secondary: Coupler reads V 
p
(n)
​	
  and V 
m
(n)
​	
 , calculates I 
PM
(n)
​	
 .

Advance Upstream: Purkinje network advances to V 
p
(n+1)
​	
  using the computed current.

Prepare Primary: Coupler re-evaluates the coupling using the latest V 
p
(n+1)
​	
  and the stale V 
m
(n)
​	
 .

Advance Primary: Myocardium advances to V 
m
(n+1)
​	
  using the updated volumetric source.

Advance Downstream: The ECG solves passively using the finalized tissue voltage.

This approach is highly computationally efficient per time step and successfully allows the Purkinje network to temporally lead the myocardium.

5. Main Problem with the Staggered Schema
While the explicit staggered approach is perfectly stable for unidirectional coupling, it poses a severe risk for bidirectional coupling.

The resistive coupling term represents a stiff numerical spring. In bidirectional mode, the massive current sink of the 3D tissue can artificially drain the 1D Purkinje nodes in a single explicit time step. If the resistance R 
pvj
​	
  is low (strong coupling), the explicit time lag between the domain updates will likely result in non-physical oscillations or catastrophic solver divergence unless the time step (Δt) is forced to be impractically small.

6. Possible Solutions Using Existing Pointers
To achieve stability in bidirectional mode without abandoning the segregated solver architecture, the current codebase provides a direct path forward via the existing pimpleControl* pimplePtr in the advance scheme signature.

Solution A: Strongly Coupled Segregated Scheme (PIMPLE Loop)

By wrapping the upstream and primary advance sequence in a while (pimplePtr->loop()) block, you can convert the explicit scheme into an implicit one.

Implementation: The segregated solvers run iteratively within a single time step. The Purkinje network and Myocardium take turns solving and updating the coupling term I 
PM
​	
  until both voltage fields converge below a specified tolerance.

Advantage: Strictly resolves the bidirectional source-sink mismatch, offering the stability of an implicit scheme while leveraging standard OpenFOAM fvMatrix machinery.

Solution B: Semi-Implicit Source Treatment

Currently, the coupling term is treated as an explicit right-hand-side vector (deposited into sourceField).

Implementation: Instead of evaluating the entire fraction explicitly, the coupler could pass the linear coefficient  
R 
pvj
​	
 
1
​	
  directly to the domains. The solvers can then add this coefficient to the main diagonal of their respective matrix systems (fvMatrix in OpenFOAM), treating the V 
p
​	
  or V 
m
​	
  variable implicitly, while keeping the neighbor's voltage on the right-hand side.

Advantage: Drastically improves the stability of the segregated scheme without necessarily requiring the multiple iterations of a full PIMPLE loop.

(Note: A true Monolithic/Block-Coupled approach, where both domains are solved simultaneously in a single massive matrix, would offer maximum stability but is generally avoided in standard OpenFOAM as it requires complex third-party block-matrix extensions).