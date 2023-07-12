from CoolProp.CoolProp import PropsSI
import CoolProp
from ambiance import Atmosphere
import numpy as np

def std_ISA():
    """
    Pegar os dados padroes da atmosfera ISA
    """
    std_atmosphere = {
        "rho_isa" : {"value" : 1.225, "unit" : "Kg/m^3"},
        "gamma" : {"value" : 1.4},
        "R_isa" : {"value" : 287.053, "unit" : "m^2/s^2K"},
        "g_isa" : {"value" : 9.81, "unit" : "m/s^2"},
        "p0" : {"value" : 101325, "unit" : "Pa"},
        "T0" : {"value" : 288.15, "unit" : "K"},
        "mi0" : {"value" : 1.789*10**-5, "unit" : "kg/ms"},
        "cp_isa" : {"value" : 1006, "unit" : "J/kgK"}
    }
    return std_atmosphere

def height_properties(height : float = None):
    """
    Calcular as propriedades atmosféricas em uma determinada altitude
    """
    properties = {
        "rho_h" : {"value" : Atmosphere(height).density, "unit" : "kg/m^3"},
        "v_sound_h" : {"value" : Atmosphere(height).speed_of_sound, "unit" : "m/s"},
        "T_h" : {"value" : Atmosphere(height).temperature, "unit" : "K"},
        "P_h" : {"value" : Atmosphere(height).pressure, "unit" : "Pa"},
        "g_h" : {"value" : Atmosphere(height).grav_accel, "unit" : "m/s^2"},
        "mu_h" : {"value" : Atmosphere(height).dynamic_viscosity, "unit" : "kg/ms"},
        "cp" : {"value" : PropsSI('C', 'T', Atmosphere(height).temperature, 'P', Atmosphere(height).pressure, 'air'), "unit" : "J/kg/K"}
    }
    return properties

def P_estag(gamma : float, T_isen: float, T_real : float, P_real : float):
    '''
    calculo da pressao de estagnacao
    args:
        gamma : constante do ar
        T_isen : temperatura isentropica
        T_real : temperatura real
        P_reeal : pressao real
    '''
    return P_real * (T_isen / T_real) ** (gamma/(gamma - 1))

def T_estag(gamma : float, mach : float, T : float):
    '''
    calculo da temperatura de estagnacao
    args:
        gamma : constante do ar
        T : temperatura
        mach : velocidade em mach
    '''
    return T * (1 + ((gamma - 1)/2) * mach ** 2)
    

def P_real(gamma : float, P_total : float, T_total : float, T_real : float):
    '''
    calculo da pressao real
    args:
        gamma: constante do ar
        P_total : pressao total
        T_total : temperatura total
        T_real : temperatura real
    '''
    return P_total*(T_real/T_total)**(gamma/(gamma-1))

def get_air_cp_gamma(T : float):
    """
    The function `get_air_cp_gamma` takes a temperature value as input and returns the specific heat
    capacity (cp) and specific heat ratio (gamma) of air at that temperature.
    
    :param T: The parameter T represents the temperature of the air
    :return: the specific heat capacity (cp) and the specific heat ratio (gamma) of air at a given
    temperature (T).
    """
    vectorT = [250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1000]
    vector_cp = [1003, 1005, 1008, 1013, 1020, 1029, 1040, 1051, 1063, 1075, 1087, 1099, 1121, 1142]
    vector_gamma = [1.401, 1.4, 1.398, 1.395, 1.391, 1.387, 1.381, 1.376, 1.37, 1.364, 1.359, 1.354, 1.344, 1.336]

    cp = np.interp(T, vectorT, vector_cp)
    gamma = np.interp(T, vectorT, vector_gamma)
    
    return cp, gamma   

if __name__ == "__main__":
    print(height_properties(1000))
    


