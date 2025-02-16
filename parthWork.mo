 model parthWork
  import Modelica.Math;

  // Constants
  constant Real R = 8.314;    // Universal gas constant (J/mol·K)
  constant Real Ea_ads = 50000;  // Activation energy for adsorption (J/mol)
  constant Real Ea_des = 80000;  // Activation energy for desorption (J/mol)
  constant Real k_ads = 1.0e3;   // Adsorption rate constant (1/s)
  constant Real k_des = 5.0e4;   // Desorption rate constant (1/s)
  constant Real deltaH = -30000; // Reaction enthalpy change (J/mol)

  // Pore Properties (Knudsen Diffusion)
  constant Real pore_radius = 5e-9;   // Pore radius (m)
  constant Real porosity = 0.4;       // Porosity (dimensionless)
  constant Real tortuosity = 2.5;     // Tortuosity (dimensionless)
  constant Real M_Cl2 = 70.906;      // Molar mass of Cl2 (g/mol)

  // Material Properties
  constant Real rho_SnO2 = 6900;  // Density of SnO2 (kg/m³)
  constant Real Cp_SnO2 = 643;    // Specific heat of SnO2 (J/kg·K)

  // Variables
  Real Cl2_conc(start=20e-6);   // Chlorine concentration (ppm)
  Real O2_ads(start=1e-6);      // Adsorbed O2 (mol/m²)
  Real O_minus_ads(start=0);    // Adsorbed O- species (mol/m²)
  Real Cl_ads(start=0);         // Adsorbed chlorine (mol/m²)
  Real T(start=300);            // Temperature (K)
  Real rate_ads;                // Adsorption rate (mol/m²·s)
  Real rate_des;                // Desorption rate (mol/m²·s)
  Real Dk;                      // Knudsen diffusivity (m²/s)
  Real Deff;                    // Effective diffusivity (m²/s)

// Equations
equation
  Dk = (2 * pore_radius / 3) * sqrt((8 * R * T) / (Modelica.Constants.pi * M_Cl2));
  Deff = porosity * Dk / tortuosity;

  // Adsorption considering O2 and O-
  rate_ads = k_ads * Deff * Cl2_conc * O2_ads * (1 - Cl_ads);
  rate_des = k_des * Deff * Cl_ads;

  // O2 to O- conversion reaction
  der(O_minus_ads) = k_ads * O2_ads - rate_ads;
  der(O2_ads) = -k_ads * O2_ads + rate_des;

  // Mass balance equations for Cl2
  der(Cl2_conc) = -rate_ads + rate_des;
  der(Cl_ads) = rate_ads - rate_des;

  // Energy balance
  der(T) = (-deltaH * (rate_ads - rate_des)) / (rho_SnO2 * Cp_SnO2);

  annotation (experiment(StartTime=0, StopTime=1000, Tolerance=1e-6));
end parthWork;
