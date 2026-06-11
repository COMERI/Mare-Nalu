Simulation:
  name: NaluSim

linear_solvers:

  - name: solve_scalar
    type: tpetra
    method: gmres
    preconditioner: sgs
    tolerance: 1e-5
    max_iterations: 50
    kspace: 50
    output_level: 0

  - name: solve_cont
    type: tpetra
    method: gmres
    preconditioner: muelu
    tolerance: 1e-5
    max_iterations: 50
    kspace: 50
    output_level: 0
    muelu_xml_file_name: ../../xml/matches_ml_default.xml

realms:

  - name: realm_2
    mesh: ../../mesh/horseshoe.g
    use_edges: no

    equation_systems:
      name: theEqSys
      max_iterations: 1

      solver_system_specification:
        velocity: solve_scalar
        pressure: solve_cont
        enthalpy: solve_scalar
        turbulent_ke: solve_scalar

      systems:

        - LowMachEOM:
            name: myLowMach
            max_iterations: 1
            convergence_tolerance: 1e-5

        - TurbKineticEnergy:
            name: myTke
            max_iterations: 1
            convergence_tolerance: 1.e-5

        - Enthalpy:
            name: myEnth
            max_iterations: 1
            convergence_tolerance: 1e-5

    material_properties:

      target_name: block_5

      constant_specification:
       universal_gas_constant: 8314.4621
       reference_pressure: 101325.0

      reference_quantities:
        - species_name: FakeAir
          mw: 28.0
          mass_fraction: 1.0

      specifications:
 
        - name: density
          type: ideal_gas

        - name: viscosity
          type: polynomial
          coefficient_declaration:
           - species_name: FakeAir
             coefficients: [1.7894e-5, 273.11, 110.56]

        - name: specific_heat
          type: polynomial
          coefficient_declaration:
           - species_name: FakeAir
             low_coefficients: [3.298677000E+00, 1.408240400E-03, -3.963222000E-06, 
                                5.641515000E-09, -2.444854000E-12,-1.020899900E+03]
             high_coefficients: [3.298677000E+00, 1.408240400E-03, -3.963222000E-06, 
                                 5.641515000E-09, -2.444854000E-12,-1.020899900E+03]

    initial_conditions:
      - constant: ic_1
        target_name: block_5
        value:
          pressure: 0
          velocity: [0,0,0]  
          temperature: 300.0 
          turbulent_ke: 1.0e-6
  
    boundary_conditions:

    - inflow_boundary_condition: bc_inflow
      target_name: surface_1
      inflow_user_data:
        velocity: [0,0,-10.0]
        temperature: 300.0
        turbulent_ke: 1.0e-3

    - open_boundary_condition: bc_side
      target_name: surface_2
      open_user_data:
        velocity: [0,0,0]
        pressure: 0.0
        temperature: 300.0
        turbulent_ke: 1.0e-6

    - wall_boundary_condition: bc_nocht
      target_name: surface_3
      wall_user_data:
        velocity: [0,0,0]
        temperature: 301.0
        use_wall_function_projected: yes
        use_neumann_condition: yes
        projected_distance: 1.0e-3

    - wall_boundary_condition: bc_cht
      target_name: surface_4
      wall_user_data:
        velocity: [0,0,0]
        temperature: 350.0
        use_wall_function_projected: yes
        use_neumann_condition: yes
        projected_distance: 1.0e-3

    solution_options:
      name: myOptions
      turbulence_model: ksgs
    
      options:
          
        - peclet_function_form:
            velocity: tanh
            enthalpy: tanh
            turbulent_ke: tanh

        - peclet_function_tanh_transition:
            velocity: 200.0
            enthalpy: 2.0
            turbulent_ke: 2.0

        - peclet_function_tanh_width:
            velocity: 20.0
            enthalpy: 4.0
            turbulent_ke: 4.0
   
        - laminar_prandtl:
            enthalpy: 1.0
            turbulent_ke: 1.0

        - turbulent_prandtl:
            enthalpy: 1.0
            turbulent_ke: 0.9

        - source_terms:
            continuity: density_time_derivative

        - limiter:
            pressure: no
            velocity: yes
            enthalpy: yes 
            turbulent_ke: yes 

    post_processing:
    
    - type: surface
      physics: surface_force_and_moment_wall_function_projected
      output_file_name: elemPipeLOWp.dat
      frequency: 1
      parameters: [0,0,0]
      target_name: [surface_3, surface_4]

    output:
      output_data_base_name: output/fluidsLOWp.e
      output_frequency: 5
      output_node_set: no
      output_variables:
       - velocity
       - pressure
       - temperature
       - temperature_bc
       - turbulent_ke
       - turbulent_viscosity
       - enthalpy
       - assembled_wall_area
       - yplus

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_time: 0.01
      time_step: 0.0001
      time_stepping_type: fixed
      time_step_count: 0
      nonlinear_iterations: 1 

      realms:
        - realm_2
