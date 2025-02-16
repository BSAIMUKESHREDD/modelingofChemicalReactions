model parth_21624
  import Modelica.Math;
  
  // Constants
  constant Real R(unit="J/(mol*K)") = 8.314;  // Universal gas constant
  constant Real Ea_ads(unit="J/mol") = 50000;  // Activation energy for adsorption
  constant Real Ea_des(unit="J/mol") = 80000;  // Activation energy for desorption
  constant Real k_ads(unit="1/s") = 1.0e3;   // Adsorption rate constant
  constant Real k_des(unit="1/s") = 5.0e4;   // Desorption rate constant
  constant Real k_react(unit="1/s") = 1.5e3; // Surface reaction rate constant
  constant Real deltaH(unit="J/mol") = -30000; // Reaction enthalpy change
  
  // Pore Properties (Knudsen Diffusion)
  constant Real pore_radius(unit="m") = 5e-9;   // Pore radius
  constant Real porosity = 0.4;       // Porosity (dimensionless)
  constant Real tortuosity = 2.5;     // Tortuosity (dimensionless)
  constant Real M_Cl2(unit="kg/mol") = 70.906e-3;   // Molar mass of Cl2
  
  // Material Properties
  constant Real rho_SnO2(unit="kg/m^3") = 6900;  // Density of SnO2
  constant Real Cp_SnO2(unit="J/(kg*K)") = 643;    // Specific heat of SnO2
  
  // Initial Conditions (molar concentrations in mol/m^3)
  Real Cl2_conc(start=20e-6*41.6, unit="mol/m^3");   // Chlorine concentration
  Real O2_conc(start=0.21*41.6, unit="mol/m^3");     // Oxygen concentration
  Real N2_conc(start=0.79*41.6, unit="mol/m^3");     // Nitrogen concentration
  
  Real O2_ads(start=1e-6, unit="mol/m^2");      // Adsorbed O2
  Real O_minus_ads(start=0, unit="mol/m^2");    // Adsorbed O- species
  Real Cl_ads(start=0, unit="mol/m^2");         // Adsorbed chlorine
  Real SnCl4_conc(start=0, unit="mol/m^3");     // Tin chloride formation
  
  Real total_sites(unit="") = 1e-6;      // Total available sites on SnO2 surface
  Real free_sites(start=1e-6, unit="");  // Free sites available for adsorption
  
  Real rate_ads(unit="mol/(m^2*s)");                // Adsorption rate
  Real rate_des(unit="mol/(m^2*s)");                // Desorption rate
  Real rate_react(unit="mol/(m^2*s)");              // Reaction rate
  Real Dk(unit="m^2/s");                      // Knudsen diffusivity
  Real Deff(unit="m^2/s");                    // Effective diffusivity
  Real T(start=300, unit="K");            // Temperature
  
  // Knudsen Diffusion
  equation
  Dk = (2 * pore_radius / 3) * sqrt((8 * R * T) / (Modelica.Constants.pi * M_Cl2));
  Deff = porosity * Dk / tortuosity;

  // Adsorption Dynamics
  rate_ads = k_ads * Deff * Cl2_conc * free_sites;
  rate_des = k_des * Cl_ads;
  
  // Surface Reaction
  rate_react = k_react * Cl_ads * O_minus_ads;
  
  // Mass Balance Equations
  der(Cl2_conc) = -rate_ads + rate_des;
  der(O2_conc) = -rate_react;
  der(O2_ads) = rate_react;
  der(O_minus_ads) = rate_ads - rate_react;
  der(Cl_ads) = rate_ads - rate_react - rate_des;
  der(SnCl4_conc) = rate_react;
  der(free_sites) = -rate_ads + rate_des + rate_react;
  der(N2_conc) = 0;
  
  // Energy Balance
  der(T) = (-deltaH * (rate_ads - rate_des)) / (rho_SnO2 * Cp_SnO2);
  
  annotation (experiment(StartTime=0, StopTime=1000, Tolerance=1e-6));
end parth_21624;
