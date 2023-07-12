
import atmosphere 
import json
from CoolProp.CoolProp import PropsSI
import re
import thermo

# The `Compressor` class is a generator object for a turbofan compressor, which calculates various
# properties such as ideal temperature, real temperature, work, and work per kilogram of fuel.
class Compressor():
    """
    The function initializes an object with given standard data, atmospheric data, and optional parameters, and
    sets attributes of the object based on the combined data.
    
    :param std_data: std_data is a dictionary containing standard data for the object being initialized.
    It could include attributes like name, age, gender, etc
    :type std_data: dict
    :param atm_data: The `atm_data` parameter is a dictionary that contains data related to an specific atmsphere.
    It could include information such as density, R.
    :type atm_data: dict
    :param ef_comp: The `ef_comp` parameter is a float that represents the efficiency component
    :type ef_comp: float
    :param rc: The parameter `rc` stands for "compressive ratio" and is a float value. It represents
    the coeficient of fluid compression
    :type rc: float
    """
    def __init__(self,
                std_data : dict,
                atm_data : dict, 
                ef_comp : float = None, 
                rc : float = None, 
                **kwargs) -> None:
        
        self.ef_comp = ef_comp
        self.rc = rc
        self.results = {}
        combined_start_data = {**std_data, **atm_data, **kwargs}

        for key, value in combined_start_data.items():
            if type(value) == dict:
                value = value['value']
            setattr(self, key, value)

    def get_ideal_temp(self):
        """
        The function calculates the ideal output temperature based on the input temperature and other
        parameters.
        :return: the value of `self.T_saida`, which is the calculated ideal temperature of the output.
        """
        self.T_saida = self.T_entrada * self.rc ** ((self.gamma - 1) / self.gamma)
        return self.T_saida
    
    def get_real_temp(self):
        """
        The function calculates the real output temperature based on the ideal output temperature, input
        temperature, and compressor efficiency.
        :return: the value of `self.T_saida`, which is the calculated real temperature of the output.
        """
        self.T_saida = self.T_entrada + ((self.T_saida - self.T_entrada) / self.ef_comp)
        return self.T_saida
    
    def get_work(self):
        """
        The function calculates the work used for compression.
        :return: the value of the variable `self.W`, which represents the work used for compression.
        """
        self.W = self.m_dot_air * self.cp * (self.T_saida - self.T_entrada)
        self.results['W'] = self.W
        return self.W
    
    def get_work_kg(self):
        """
        The function calculates the energy gain per kilogram of fuel.
        :return: the energy gain per kilogram of fuel, which is calculated using the specific heat capacity
        (cp) and the temperature difference between the outlet (T_saida) and the inlet (T_entrada).
        """
        self.W_kg = self.cp * (self.T_saida - self.T_entrada)
        self.results['W_kg'] = self.W_kg
        return self.W_kg
    
    def step(self, **kwargs):
        """
        The `step` function takes in keyword arguments, sets the corresponding attributes of the object,
        performs some calculations, and returns a dictionary of results.
        :return: a dictionary containing the following key-value pairs:
        - "P_entrada" : self.P_entrada
        - "P_saida" : self.P_saida
        - "T_entrada" : self.T_entrada
        - "T_saida" : self.T_saida
        """
        for key, value in kwargs.items():
            if type(value) == dict:
                value = value['value']
            setattr(self, key, value)

        self.P_saida = self.P_entrada * self.rc
        self.get_ideal_temp()
        self.get_real_temp()
        self.get_work_kg()
        self.get_work()

        results = {"P_entrada" : self.P_entrada, 
                   "P_saida" : self.P_saida,
                   "T_entrada" : self.T_entrada,
                   "T_saida" : self.T_saida}
        self.results.update(results)

        return self.results

# The Turbine class is a generator object for a turbine component, with attributes and methods for
# initializing and setting turbine properties.
class Turbine():
    """
    The `__init__` function initializes an object with provided data and sets attributes based on the
    data.
    
    :param std_data: std_data is a dictionary containing standard data. It is used to initialize
    attributes of the class object
    :type std_data: dict
    :param atm_data: The `atm_data` parameter is a dictionary that contains data related to the
    atmosphere. It is used to initialize the `self` object in the `__init__` method of a class
    :type atm_data: dict
    :param ef_turb: The `ef_turb` parameter is a float that represents the efficiency of the turbine
    :type ef_turb: float
    :param rc: The `rc` parameter is a float that represents the value of the compression rate attribute
    :type rc: float
    """
    def __init__(self,
                std_data : dict,
                atm_data : dict, 
                ef_turb : float = None, 
                rc : float = None, 
                **kwargs) -> None:
        
        self.ef_turb = ef_turb
        self.rc = rc
        self.results = {}
        combined_start_data = {**std_data, **atm_data, **kwargs}

        for key, value in combined_start_data.items():
            if type(value) == dict:
                value = value['value']
            setattr(self, key, value)
    
    def get_temp_from_work(self):
        """
        The function calculates the temperature and pressure at the outlet based on the input parameters.
        :return: The function is not returning any value.
        """

        pattern = r"\d+"
        matches = re.findall(pattern, self.id)
        id_number = int(matches[0]) if matches else None

        self.T_saida = self.T_entrada - (self.W_requerido / (self.ef_turb*(self.m_dot_air + self.m_dot_f)*self.cp))
        self.T_iso = self.T_entrada - (self.delta_T_compressor[id_number]/self.ef_turb)

        self.P_saida = self.P_entrada / ((self.T_entrada/self.T_iso)**(self.gamma/(self.gamma-1)))
        return 
        
   
    def get_iso_temp(self):
        """
        The function calculates the ideal temperature based on the input values.
        :return: the value of `self.T_iso_saida`, which is the calculated ideal temperature.
        """
        self.T_iso_saida = self.T_entrada / ((self.P_entrada / self.P_saida)**((self.gamma - 1)/self.gamma))
        return self.T_iso_saida
    
    def get_real_temp(self):
        """
        The function calculates the real output temperature based on the input temperature, isolation
        temperature, and turbine efficiency.
        :return: the value of the variable `self.T_saida`, which represents the calculated real
        temperature of the output.
        """
        self.T_saida = self.T_entrada - (self.T_entrada - self.T_iso_saida)/self.ef_turb
        return self.T_saida
    

    def get_work(self):
        """
        The function calculates the work generated by a turbine.
        :return: The variable `W` is being returned.
        """
        if not hasattr(self, 'm_dot_air'):
            print('To calculate the total work W is required a m_dot_air in turbine call')
            self.m_dot_air = 1
        self.W = self.m_dot_air * self.cp * (self.T_saida - self.T_entrada)
        self.results['W'] = self.W
        return self.W
    
    def get_expand_ratio(self):
        """
        The function calculates the expansion ratio based on the input and output temperatures.
        :return: The method is returning the expand ratio, which is calculated as `self.expand_ratio =
        cte**ratio`.
        """
        ratio = self.T_entrada/self.T_saida
        cte = self.gamma/(self.gamma - 1)
        self.expand_ratio = cte**ratio
        self.results['expand ratio'] = self.expand_ratio
        return self.expand_ratio
    
    def get_work_kg(self):
        """
        The function calculates the specific work of a turbine based on the specific heat capacity,
        temperature difference, and efficiency.
        :return: the value of the variable `self.W_kg`.
        """
        self.W_kg = self.cp * (self.T_saida - self.T_entrada)
        self.results['W_kg'] = self.W_kg
        return self.W_kg
    
    def get_pressure(self):
        """
        The function calculates the output pressure of a turbine based on the input pressure and the
        expansion ratio.
        :return: the value of `self.P_saida`, which is the calculated output pressure of the turbine.
        """
        if not hasattr(self, "expand_ratio"):
            self.get_expand_ratio()
        self.P_saida = self.P_entrada/self.expand_ratio
        return self.P_saida
    
    def get_mach_jato(self):
        """
        The function calculates the Mach number of the turbine's exit velocity.
        :return: the value of the variable `self.M`, which represents the Mach number of the turbine's exit
        velocity.
        """
        a = 2/(self.gamma - 1)
        b = (self.P_entrada / self.P_h)**((self.gamma - 1)/self.gamma)
        self.M = a*(b - 1)**0.5
        self.results['M_saida'] = self.M
        return self.M
    
    def get_jet_vel(self):
        """
        The function calculates the jet velocity and propulsive efficiency based on given parameters.
        :return: the value of `self.v_jet`.
        """

        a = 1 - (self.P_saida / self.P_entrada) ** ((self.gamma - 1)/self.gamma)
        v = 2 * self.cp * self.T_entrada * a
        self.v_jet = v ** 0.5
        FG = self.v_jet
        FN = self.v_jet - self.aircraft_vel
        self.eta_prop = (2 * self.aircraft_vel)/(self.aircraft_vel + self.v_jet)
        
        results = {'v_jet' : self.v_jet, 'eta_prop' : self.eta_prop}
        self.results.update(results)
        return self.v_jet

    
    def step(self, **kwargs):
        """
        The `step` function takes in keyword arguments, sets the corresponding attributes of the object, and
        then calculates and returns various results.
        :return: a dictionary containing the values of "P_entrada", "P_saida", "T_entrada", and "T_saida".
        """
        for key, value in kwargs.items():
            if type(value) == dict:
                value = value['value']
            setattr(self, key, value)
        
        if self.rc is not None:
            self.P_saida = self.P_entrada / self.rc
            self.get_iso_temp()
            self.get_real_temp()
            self.get_expand_ratio()
            self.get_pressure()
        else:
            self.get_temp_from_work()
        
        self.get_work_kg()
        self.get_work()
        self.get_jet_vel()

        results = {"P_entrada" : self.P_entrada, 
                   "P_saida" : self.P_saida,
                   "T_entrada" : self.T_entrada,
                   "T_saida" : self.T_saida}
        self.results.update(results)
        
        return self.results
    
#TODO rever o calculo para acrescentar as propriedades de estagnacao
# The `Diffuser` class initializes an object with standard and atmospheric data, and sets attributes
# based on the combined data.
class Diffuser():
    """
    The `__init__` function initializes an object with provided data and sets attributes based on the
    data.
    
    :param std_data: std_data is a dictionary containing standard data. It is used to initialize
    attributes of the class object
    :type std_data: dict
    :param atm_data: The `atm_data` parameter is a dictionary that contains data related to the
    atmosphere. It is used to initialize the `self` object in the `__init__` method of a class
    :type atm_data: dict
    """
    def __init__(self,
                std_data : dict,
                atm_data : dict, 
                **kwargs) -> None:
        
        self.results = {}
        combined_start_data = {**std_data, **atm_data, **kwargs}
        for key, value in combined_start_data.items():
            if type(value) == dict:
                value = value['value']
            setattr(self, key, value)
    
    def get_temperature(self):
        """
        The function calculates the temperature at the exit of a flow based on the input temperature, Mach
        number, and specific heat ratio.
        :return: the value of `self.T_saida`, which is the calculated temperature.
        """
        self.T_saida = float(self.T_entrada * (1 + (self.gamma - 1) * self.M**2/2))
        return self.T_saida

    def get_pressure_subsonic(self):
        """
        The function calculates the pressure at the exit of a subsonic flow using the given inputs.
        :return: the value of `self.P_saida`.
        """
        self.P_saida = float((self.T_saida/self.T_entrada)**(self.gamma/(self.gamma - 1))*self.P_entrada)
        return self.P_saida
    
    def step(self, **kwargs):
        """
        The function updates the attributes of an object based on the provided keyword arguments, retrieves
        temperature and pressure values, and returns a dictionary of the updated results.
        :return: the updated results dictionary.
        """

        for key, value in kwargs.items():
            if type(value) == dict:
                value = value['value']
            setattr(self, key, value)

        self.get_temperature()
        self.get_pressure_subsonic()

        results = {"P_entrada" : self.P_entrada, 
                   "P_saida" : self.P_saida,
                   "T_entrada" : self.T_entrada,
                   "T_saida" : self.T_saida}
        self.results.update(results)   

        return self.results 

class Combustion_Chamber():
    '''
    The `Combustion_Chamber` class represents a combustion chamber. It has the following attributes:
    - `std_data`: A dictionary containing standard data.
    - `atm_data`: A dictionary containing atmospheric data.
    - `eta`: An optional float representing the combustion efficiency.
    - `pi`: An optional float representing the combustion chamber's resistance coefficient.
    The `__init__` method initializes the `Combustion_Chamber` object. It takes in the `std_data` and `atm_data` dictionaries, as well as the optional `eta` and `pi` values. Additional keyword arguments can also be passed.
    These additional keyword arguments (**kwargs) can include any other properties that you want to add to the `Combustion_Chamber` object, such as the fuel type (`fuel_type`), inlet temperature (`T_entrada`), inlet pressure (`P_entrada`), air mass flow rate (`m_dot_air`), among others.
    The method combines the `std_data`, `atm_data`, and additional keyword arguments into the `combined_start_data` dictionary. It then iterates over the items in this dictionary and sets each key-value pair as an attribute of the object. If the value is a dictionary, it extracts the 'value' key from it.
    '''

    def __init__(self,
                std_data : dict,
                atm_data : dict, 
                eta : float = None, 
                pi : float = None, 
                **kwargs) -> None:
        
        self.eta = eta
        self.pi = pi
        self.results = {}
        combined_start_data = {**std_data, **atm_data, **kwargs}

        for key, value in combined_start_data.items():
            if type(value) == dict:
                value = value['value']
            setattr(self, key, value)

    def get_pressure(self):
        """
        The function calculates the output pressure based on the input pressure and a constant factor.
        :return: the value of `self.P_saida`.
        """

        self.P_saida = self.P_entrada * self.pi
        return self.P_saida
    
    def get_mf_from_T_saida(self):
        """
        The function calculates the fuel mass flow rate of a fluid based on its properties and the temperature
        difference between the inlet and outlet.
        :return: the value of `self.m_f`.
        """
        if hasattr(self, 'fuel_type'):
            try:
                self.cp_fuel = PropsSI('C', 'T', self.T_saida, 'P', self.P_saida, self.fuel_type)
            except:
                raise ValueError(f'Your fluid {self.fuel_type} are not present in CoolProp ref, please add the Cp value (J/kg) with the argument cp_fuel')

        self.m_dot_f = self.m_dot_air*self.cp*(self.T_saida - self.T_entrada) / (self.LCV*self.eta - self.cp_fuel*self.T_saida)
        self.results['m_dot_fuel'] = self.m_dot_f
        return self.m_dot_f
    
    def get_T_saida_from_mf(self):
        return NotImplemented

    def step(self, **kwargs):
        """
        The `step` function updates the attributes of an object based on the provided keyword arguments,
        calculates pressure and temperature values, and returns a dictionary of the results.
        :return: a dictionary containing the values of "P_entrada", "P_saida", "T_entrada", and "T_saida".
        """
        for key, value in kwargs.items():
            if type(value) == dict:
                value = value['value']
            setattr(self, key, value)
        
        self.get_pressure()
        if hasattr(self, 'T_saida'):        
            self.get_mf_from_T_saida()
        else:
            self.get_T_saida_from_mf()
        
        results = {"P_entrada" : self.P_entrada, 
                   "P_saida" : self.P_saida,
                   "T_entrada" : self.T_entrada,
                   "T_saida" : self.T_saida}
        
        self.results.update(results)
        return self.results
    
class Fan():
    
    """
    The above code defines a Python class called "Fan" that represents a fan and provides methods to
    calculate outlet pressure, outlet temperature, and generated power.
    
    :param std_data: std_data is a dictionary containing standard data for the fan. It may include the
    following keys:
    :type std_data: dict
    :param atm_data: The `atm_data` parameter is a dictionary that contains atmospheric data. It could
    include values such as atmospheric pressure, atmospheric temperature, and specific heat ratio. These
    values are used in the calculations performed by the methods of the `Fan` class
    :type atm_data: dict
    :param eta: Efficiency of the fan. It is a float value
    :type eta: float
    :param pi: The parameter "pi" represents the constant value used to calculate the outlet pressure of
    the fan
    :type pi: float
    :param split_ratio: The split_ratio parameter represents the ratio of the mass flow rate of the
    bypass air to the total mass flow rate of the fan. It is used in the calculation of the generated
    power (W) by the fan
    :type split_ratio: float
    :param bypass: The "bypass" parameter is a float value that represents the bypass ratio of the fan.
    It is used in the calculation of the power generated by the fan (method get_W()). The bypass ratio
    is the ratio of the mass flow rate of the bypass air to the mass flow rate of the primary
    :type bypass: float
    """

    def __init__(self,
                std_data : dict,
                atm_data : dict, 
                eta : float = None, 
                pi : float = None, 
                split_ratio : float = None,
                bypass : float = None,
                **kwargs) -> None:
        
        self.eta = eta
        self.pi = pi
        self.split_ratio = split_ratio
        self.bypass = bypass
        self.results = {}
        combined_start_data = {**std_data, **atm_data, **kwargs}

        for key, value in combined_start_data.items():
            if type(value) == dict:
                value = value['value']
            setattr(self, key, value)

    def get_pressure(self):
        """
        The function calculates the output pressure based on the input pressure and a constant value.
        :return: the value of `self.P_saida`.
        """
        self.P_saida = self.pi*self.P_entrada
        return self.P_saida
    
    def get_temperature(self):
        """
        The function calculates the temperature at the outlet based on the input temperature, specific heat
        ratio, and efficiency.
        :return: the value of `self.T_saida`.
        """
        Tis = self.T_entrada*self.pi**((self.gamma-1)/self.gamma)
        self.T_saida = self.T_entrada + (Tis-self.T_entrada)/self.eta
        return self.T_saida
    
    def get_W(self):
        """
        The function calculates and returns the value of W using the given formula.
        :return: The method is returning the value of the variable `self.W`.
        """
        self.W = self.bypass * self.m_dot_air * self.cp * (self.T_saida-self.T_entrada)
        self.results['W'] = self.W
        return self.W
    
    def step(self, **kwargs):
        """
        The `step` function sets the attributes of an object based on the provided keyword arguments,
        retrieves pressure and temperature values, and returns a dictionary of the results.
        :return: a dictionary containing the values of "P_entrada", "P_saida", "T_entrada", and "T_saida".
        """

        for key, value in kwargs.items():
            if type(value) == dict:
                value = value['value']
            setattr(self, key, value)

        self.get_pressure()
        self.get_temperature()
        if self.bypass != None:
            self.get_W()

        results = {"P_entrada" : self.P_entrada, 
                   "P_saida" : self.P_saida,
                   "T_entrada" : self.T_entrada,
                   "T_saida" : self.T_saida}
        
        self.results.update(results)
        return self.results
    
# The Nozzle class calculates various parameters related to fluid flow based on given inputs.
class Nozzle():
    """
    The above function initializes the class with given parameters and calculates various fluid flow
    parameters based on the inputs.
    
    :param std_data: std_data is a dictionary containing standard data related to fluid flow. It may
    include the following keys:
    :type std_data: dict
    :param atm_data: The `atm_data` parameter is a dictionary that contains atmospheric data. It is used
    to provide additional information about the atmospheric conditions for the fluid flow calculations
    :type atm_data: dict
    :param eta_n: The parameter `eta_n` represents the isentropic efficiency of the system. It is used
    in the calculations to determine the pressure and temperature at the exit of the fluid flow
    :type eta_n: float
    :param pi_duct: `pi_duct` is a dimensionless parameter that represents the pressure loss in the
    duct. It is used to calculate the entrance pressure (`P_entrada`) by multiplying it with the initial
    pressure (`P_entrada`)
    :type pi_duct: float
    """
    def __init__(self,
                std_data : dict,
                atm_data : dict, 
                eta_n : float = None, 
                pi_duct : float = None, 
                **kwargs) -> None:

        self.eta_n = eta_n
        self.pi_duct = pi_duct
        self.results = {}
        combined_start_data = {**std_data, **atm_data, **kwargs}

        for key, value in combined_start_data.items():
            if type(value) == dict:
                value = value['value']
            setattr(self, key, value)
            

    def get_parameters(self):
        """
        The function calculates various parameters related to fluid flow based on given inputs.
        """
        self.P_entrada = self.P_entrada * (1 - self.pi_duct)
        self.P_crit = self.P_entrada * (1 - ((self.gamma - 1)/(self.gamma + 1)/self.eta_n))**(self.gamma/(self.gamma - 1))

        if self.P_crit >= self.p0:
            self.P_saida = self.P_crit
            self.T_saida = self.T_entrada * (1-self.eta_n*(1 - (self.P_saida/self.P_entrada)**((self.gamma - 1)/self.gamma)))
            self.v_jet = (2 * self.cp_fuel * (self.T_entrada - self.T_saida))**0.5

        else:
            self.P_saida = self.p0

        self.results['v_jet'] = self.v_jet
        self.T_saida = self.T_entrada * (1-self.eta_n*(1 - (self.P_saida/self.P_entrada)**((self.gamma - 1)/self.gamma)))
        cp_saida, _ = atmosphere.get_air_cp_gamma(T= self.T_saida)
        self.v_jet = (2 * cp_saida  * (self.T_entrada - self.T_saida))**0.5

    def step(self, **kwargs):
        """
        The `step` function updates the attributes of an object based on the provided keyword arguments,
        retrieves the parameters, updates the results dictionary, and returns the updated results.
        :return: the updated results dictionary.
        """

        for key, value in kwargs.items():
            if type(value) == dict:
                value = value['value']
            setattr(self, key, value)

        self.get_parameters()

        results = {"P_entrada" : self.P_entrada, 
                   "P_saida" : self.P_saida,
                   "T_entrada" : self.T_entrada,
                   "T_saida" : self.T_saida}
        
        self.results.update(results)
        return self.results
    
class Mixer():
    """
    The function initializes an object with given data and sets attributes based on the data.
    
    :param std_data: std_data is a dictionary containing standard data values
    :type std_data: dict
    :param atm_data: A dictionary containing atmospheric data such as temperature, pressure, humidity,
    etc
    :type atm_data: dict
    :param eta_m: The `eta_m` represents the isentropic efficiency of the system. It is used
    in the calculations to determine the pressure and temperature at the exit of the fluid flow
    :type eta_m: float
    :param pi: The `pi` is a dimensionless parameter that represents the pressure loss in the
    duct. It is used to calculate the entrance pressure (`P_entrada`) by multiplying it with the initial
    pressure (`P_entrada`)
    :type pi: float
    """
    def __init__(self,
                std_data : dict,
                atm_data : dict, 
                eta_m : float = None, 
                pi : float = None, 
                **kwargs) -> None:

        self.eta_m = eta_m
        self.pi = pi
        self.results = {}
        combined_start_data = {**std_data, **atm_data, **kwargs}

        for key, value in combined_start_data.items():
            if type(value) == dict:
                value = value['value']
            setattr(self, key, value)

    def get_mass_flow(self):
        """
        The function calculates the mass flow rate of air and updates the results dictionary with the value.
        :return: the value of the variable `self.m_dot_air`.
        """
        self.m_dot_air = self.m_dot_air*self.bypass*self.split_ratio
        self.results['m_dot_air'] = self.m_dot_air
        return self.m_dot_air
    
    def get_temperature(self):
        """
        The function calculates and returns the temperature of the output based on the given inputs.
        :return: the value of `self.T_saida`, which is the calculated temperature.
        """
        self.cp_bypass, _ = atmosphere.get_air_cp_gamma(self.T_bypass)
        self.cp_intern, _ = atmosphere.get_air_cp_gamma(self.T_entrada)
        self.T_saida = (self.m_dot_air * self.T_bypass * self.cp_bypass + ( self.m_dot_air + self.m_dot_f ) * self.T_entrada * self.cp_intern  ) / ( self.m_dot_air * self.cp_bypass + ( self.m_dot_air + self.m_dot_f ) * self.cp_intern )
        return self.T_saida

    def get_pressure(self):
        """
        The function calculates the output pressure based on the input pressure and a constant factor.
        :return: the value of `self.P_saida`.
        """
        self.P_saida = self.P_entrada*self.pi
        return self.P_saida

    def step(self, **kwargs):
        """
        The function "step" is used to retrieve the mass flow, temperature, and pressure.
        """

        for key, value in kwargs.items():
            if type(value) == dict:
                value = value['value']
            setattr(self, key, value)

        self.get_mass_flow()
        self.get_temperature()
        self.get_pressure()

        results = {"P_entrada" : self.P_entrada, 
                   "P_saida" : self.P_saida,
                   "T_entrada" : self.T_entrada,
                   "T_saida" : self.T_saida}
        
        self.results.update(results)
        return self.results


class Afterburner():
    """
    The  __init__  method in the  Afterburner  class serves as the constructor for objects of this class. It initializes the object's attributes and sets their initial values based on the provided arguments and keyword arguments. 
 
    The method takes in several arguments: 
    -  std_data  and  atm_data  are dictionaries containing standard data and atmospheric data, respectively. 
    -  pi  is an optional float argument representing a pressure ratio. 
    
    The method initializes the  pi  attribute with the provided value and creates an empty  results  dictionary attribute. 
    
    Next, it combines the  std_data ,  atm_data , and additional keyword arguments ( kwargs ) into a single dictionary called  combined_start_data  using the  **  operator. 
    
    Then, in a loop, it iterates over the key-value pairs in  combined_start_data . If the value is of type  dict , it extracts the  'value'  key from the dictionary and assigns it to the corresponding attribute. This allows for flexible assignment of attribute values, accommodating nested dictionaries if present. 
    
    Overall, the  __init__  method sets up the initial state of an  Afterburner  object by assigning attribute values based on the provided data and keyword arguments.
    """
    def __init__(self,
                std_data : dict,
                atm_data : dict, 
                eta_ab : float,
                pi : float = None, 
                **kwargs) -> None:

        self.pi = pi
        self.eta_ab = eta_ab
        self.results = {}
        combined_start_data = {**std_data, **atm_data, **kwargs}

        for key, value in combined_start_data.items():
            if type(value) == dict:
                value = value['value']
            setattr(self, key, value)

    def get_temperature(self):
        """
        The function calculates the temperature of the outlet based on the given inputs.
        :return: the value of `self.T_saida`.
        """
        cp, _ = atmosphere.get_air_cp_gamma(T= self.T_entrada)
        self.T_saida = (self.m_dot_fab * (self.eta_ab * self.LCV - cp * self.T_entrada) + (self.m_dot_air + self.m_dot_f) * cp * self.T_entrada)/(cp * (self.m_dot_air + self.m_dot_f - self.m_dot_fab))
        return self.T_saida
    
    def get_pressure(self):
        """
        The function calculates the output pressure based on the input pressure and a constant factor.
        :return: the value of `self.P_saida`.
        """
        self.P_saida = self.P_entrada * self.pi
        return self.P_saida
    
    def step(self, **kwargs):
        """
        The function updates the attributes of an object with the values provided in the kwargs dictionary,
        retrieves temperature and pressure values, updates the results dictionary, and returns the updated
        results.
        :return: the updated results dictionary, which contains the values of P_entrada, P_saida, T_entrada,
        and T_saida.
        """

        for key, value in kwargs.items():
            if type(value) == dict:
                value = value['value']
            setattr(self, key, value)       

        self.get_temperature()
        self.get_pressure()

        results = {"P_entrada" : self.P_entrada, 
                   "P_saida" : self.P_saida,
                   "T_entrada" : self.T_entrada,
                   "T_saida" : self.T_saida}
        
        self.results.update(results)
        return self.results
    

if __name__ == "__main__":
    std = atmosphere.std_ISA()
    atm = atmosphere.height_properties(height= 1000)

    args = {'T_entrada' : 257,
            'P_entrada' : 97000,
            'aircraft_vel' : 600}
    compressor = Compressor(ef_comp= 0.9, rc= 45, std_data= std, atm_data= atm, **args)
    print(compressor.step())