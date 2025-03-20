!< PRISM, general parameters.
module adam_prism_parameters
    !< PRISM, general parameters.
    
    use PENF
    
    implicit none
    save
    private
    public :: mu0, eps0     
    
    real(R8P), parameter :: mu0      = 1.256637e-6_R8P          !< vacuum magnetic permeability
    real(R8P), parameter :: eps0     = 8.854187e-12_R8P         !< vacuum dielectric constant
    real(R8P), parameter :: e_charge = -1.6e-19_R8P             !< Elettron charge value
    real(R8P), parameter :: e_mass   = 9.11e-31_R8P             !< Elettron mass value
    real(R8P), parameter :: q_over_m = e_charge/e_mass          !< q/m value for elettron
    real(R8P), parameter :: k_B      = 1.380649e-23_R8P         !< Boltzmann coefficient
    real(R8P), parameter :: T_e      = 300000_R8P               !< Plasma Temperature
    real(R8P), parameter :: vth      = sqrt(k_B * T_e / e_mass) !< Thermal Velocity value
    real(R8P), parameter :: pi       = 3.141592653589793_R8P    !< Pi value

endmodule adam_prism_parameters