import numpy

def D_ram(m_dot_air : float, v_in : float):
    """
    The function calculates the ram drag force based on the mass flow rate of air and the velocity of
    the air.
    
    :param m_dot_air: The parameter "m_dot_air" represents the mass flow rate of air. It is a float
    value
    :type m_dot_air: float
    :param v_in: The parameter "v_in" represents the velocity of the air entering a system
    :type v_in: float
    :return: the product of m_dot_air and v_in, which is the value of F.
    """
    F = m_dot_air * v_in
    return F

def emp_stat(m_dot_air : float, m_dot_fuel : float, v_jet : float):
    """
    The function calculates the static thrust force produced by an engine based on the mass flow rates of air
    and fuel and the velocity of the jet.
    
    :param m_dot_air: The parameter m_dot_air represents the mass flow rate of air
    :type m_dot_air: float
    :param m_dot_fuel: The parameter m_dot_fuel represents the mass flow rate of the fuel
    :type m_dot_fuel: float
    :param v_jet: The parameter `v_jet` represents the velocity of the jet or exhaust gases
    :type v_jet: float
    :return: the value of the variable F, which is the product of (m_dot_air + m_dot_fuel) and v_jet.
    """
    F = (m_dot_air + m_dot_fuel) * v_jet
    return F

def emp(m_dot_air : float, m_dot_fuel : float, v_jet : float, v_in : float):
    """
    The function calculates the net thrust force of an engine based on the mass flow rates of air and
    fuel, the velocity of the jet, and the velocity of the incoming air.
    
    :param m_dot_air: The parameter `m_dot_air` represents the mass flow rate of air
    :type m_dot_air: float
    :param m_dot_fuel: The parameter "m_dot_fuel" represents the mass flow rate of fuel
    :type m_dot_fuel: float
    :param v_jet: The parameter `v_jet` represents the velocity of the jet exhaust
    :type v_jet: float
    :param v_in: The parameter `v_in` represents the velocity of the air entering the engine
    :type v_in: float
    :return: the value of F, which represents the net force produced by the engine.
    """
    stat = emp_stat(m_dot_air= m_dot_air,
                    m_dot_fuel= m_dot_fuel, 
                    v_jet= v_jet)
    ram = D_ram(m_dot_air= m_dot_air,
                v_in= v_in)
    F = stat - ram
    return F

def ene_cin(m_dot_air : float, m_dot_fuel : float, v_jet : float, v_in : float):
    """
    The function calculates the energy produced by a jet engine based on the mass flow rates of air and
    fuel, the velocity of the jet, and the velocity of the incoming air.
    
    :param m_dot_air: The parameter "m_dot_air" represents the mass flow rate of air
    :type m_dot_air: float
    :param m_dot_fuel: The parameter m_dot_fuel represents the mass flow rate of the fuel
    :type m_dot_fuel: float
    :param v_jet: The parameter "v_jet" represents the velocity of the jet stream, which is the speed at
    which the air-fuel mixture is expelled from the engine
    :type v_jet: float
    :param v_in: The parameter `v_in` represents the velocity of the air/fuel mixture before it enters
    the combustion chamber
    :type v_in: float
    :return: the value of E, which is the energy produced by the engine.
    """
    E = ((m_dot_air + m_dot_fuel) * (v_jet ** 2) - m_dot_air * (v_in ** 2))/2
    return E

def pow_aircraft(m_dot_air : float, m_dot_fuel : float, v_jet : float, v_in : float):
    """
    The function calculates the power output of an aircraft engine based on the mass flow rates of air
    and fuel, the jet velocity, and the inlet velocity.
    
    :param m_dot_air: The mass flow rate of air (in kg/s)
    :type m_dot_air: float
    :param m_dot_fuel: The parameter m_dot_fuel represents the mass flow rate of fuel
    :type m_dot_fuel: float
    :param v_jet: The parameter "v_jet" represents the velocity of the jet exhaust
    :type v_jet: float
    :param v_in: The parameter `v_in` represents the velocity of the incoming air
    :type v_in: float
    :return: the power output of the aircraft.
    """
    p = emp(m_dot_air= m_dot_air,
            m_dot_fuel= m_dot_fuel,
            v_jet= v_jet,
            v_in= v_in)
    p = v_in*p
    return p

def prop_ef(m_dot_air : float, m_dot_fuel : float, v_jet : float, v_in : float, LCV : float):
    """
    The function calculates the propulsive efficiency of an aircraft based on the given parameters.
    
    :param m_dot_air: The mass flow rate of air (in kg/s)
    :type m_dot_air: float
    :param m_dot_fuel: The parameter `m_dot_fuel` represents the mass flow rate of fuel
    :type m_dot_fuel: float
    :param v_jet: The parameter "v_jet" represents the velocity of the jet exhaust
    :type v_jet: float
    :param v_in: The parameter "v_in" represents the velocity of the incoming air or fluid
    :type v_in: float
    :param LCV: LCV stands for Lower Calorific Value, which is the amount of heat energy released per
    unit mass of fuel when it is completely burned. It is typically measured in units of energy per unit
    mass, such as joules per kilogram (J/kg) or British thermal units per pound (BT
    :type LCV: float
    :return: the propulsive efficiency (np) of an aircraft.
    """

    p = pow_aircraft(m_dot_air= m_dot_air,
                     m_dot_fuel= m_dot_fuel,
                     v_jet= v_jet,
                     v_in= v_in)
    
    E = ene_cin(m_dot_air= m_dot_air,
                m_dot_fuel= m_dot_fuel, 
                v_jet= v_jet,
                v_in= v_in)
    np = p/E
    return np

def term_ef(m_dot_air : float, m_dot_fuel : float, v_jet : float, v_in : float, LCV : float):
    """
    The function calculates the thermal efficiency of a jet engine based on the mass flow rates of air
    and fuel, jet velocity, inlet velocity, and the lower calorific value of the fuel.
    
    :param m_dot_air: The mass flow rate of air (in kg/s)
    :type m_dot_air: float
    :param m_dot_fuel: The parameter `m_dot_fuel` represents the mass flow rate of the fuel in kg/s
    :type m_dot_fuel: float
    :param v_jet: The parameter "v_jet" represents the velocity of the jet stream, which is the speed at
    which the air-fuel mixture is expelled from the engine
    :type v_jet: float
    :param v_in: The parameter "v_in" represents the velocity of the air/fuel mixture entering the
    combustion chamber
    :type v_in: float
    :param LCV: LCV stands for Lower Calorific Value, which is the amount of heat released per unit mass
    of fuel when it is completely burned. It is typically measured in units of energy per unit mass,
    such as joules per kilogram (J/kg) or megajoules per kilogram (
    :type LCV: float
    :return: the thermal efficiency (nth) of a jet engine.
    """
    E = ene_cin(m_dot_air= m_dot_air, 
                m_dot_fuel= m_dot_fuel, 
                v_jet= v_jet, 
                v_in= v_in)
    
    nth = E/(m_dot_fuel * LCV)
    return nth

def tot_ef(m_dot_air : float, m_dot_fuel : float, v_jet : float, v_in : float, LCV : float):
    """
    The function calculates the total efficiency of a jet engine based on the mass flow rates of air and
    fuel, jet velocity, inlet velocity, and lower calorific value.
    
    :param m_dot_air: The mass flow rate of air (in kg/s)
    :type m_dot_air: float
    :param m_dot_fuel: The parameter "m_dot_fuel" represents the mass flow rate of the fuel
    :type m_dot_fuel: float
    :param v_jet: The parameter "v_jet" represents the velocity of the jet stream
    :type v_jet: float
    :param v_in: The parameter `v_in` represents the velocity of the air/fuel mixture entering the
    combustion chamber
    :type v_in: float
    :param LCV: The LCV stands for Lower Calorific Value. It is the amount of heat released per unit
    mass of fuel when it is completely burned. It is typically measured in units of energy per unit
    mass, such as joules per kilogram (J/kg) or British thermal units per pound (BT
    :type LCV: float
    :return: The function `tot_ef` returns the total efficiency `nt` of a propulsion system.
    """
    np = prop_ef(m_dot_air= m_dot_air,
                 m_dot_fuel= m_dot_fuel,
                 v_jet= v_jet,
                 v_in= v_in)
    nth = term_ef(m_dot_air= m_dot_air,
                  m_dot_fuel= m_dot_fuel,
                  v_jet= v_jet,
                  v_in= v_in)
    
    nt = np* nth
    return nt

def get_density(P : float, T : float, R : float = 287.1):
    """
    The function calculates the density of a gas given its pressure, temperature, and gas constant.
    
    :param P: Pressure (in pascals)
    :type P: float
    :param T: The parameter T represents the temperature in Kelvin
    :type T: float
    :param R: The parameter R represents the ideal gas constant, which has a value of 287.1 J/(kg·K)
    :type R: float
    :param M: The parameter M represents the molar mass of the gas
    :type M: float
    :return: the density (rho) calculated using the given formula.
    """
    rho = P /(R * T)
    return rho