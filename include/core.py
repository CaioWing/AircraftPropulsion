import numpy as np
import atmosphere 
import json
import thermo

import structures

def truncate_dict_values(data):
    """
    The function `truncate_dict_values` recursively truncates float values in a dictionary to two
    decimal places.
    
    :param data: The `data` parameter is a dictionary that contains key-value pairs. The values can be
    of any data type, including nested dictionaries
    :return: a modified version of the input dictionary where all float values are rounded to 2 decimal
    places.
    """
    if isinstance(data, dict):
        return {key: truncate_dict_values(value) for key, value in data.items()}
    elif isinstance(data, float):
        return round(data, 3)
    else:
        return data
    
def check_requirements(list_requirements : list, attrb_dict : dict):
    """
    The function "check_requirements" checks if a list of required attributes is present in a dictionary
    and prints the attributes that are missing.
    
    :param list_requirements: A list of required attributes that need to be checked in the attrb_dict
    dictionary
    :type list_requirements: list
    :param attrb_dict: The `attrb_dict` parameter is a dictionary that contains attributes as keys and
    their corresponding values
    :type attrb_dict: dict
    """

    add_attrb = []
    for required_attrb in list_requirements:
        if required_attrb not in attrb_dict.keys():
            add_attrb.append(required_attrb)
    print(f'Os seguintes atributos nao foram identificados: {add_attrb}')

# The Motor class represents a motor system and allows for the addition of components such as
# compressors and turbines, solving the system, and saving the simulation data.
class Motor():

    """
    The above function is an initializer for a class that sets attributes based on input parameters
    and checks for required properties.
    
    :param height: The height parameter represents the height at which the properties of the
    atmosphere are being calculated. It is a float value representing the height in meters, defaults
    to 0
    :type height: float (optional)
    :param number_lines: The parameter "number_lines" is an integer that represents the number of
    lines, defaults to 1
    :type number_lines: int (optional)
    :param properties: The "properties" parameter is a dictionary that allows you to pass additional
    attributes and their corresponding values to the object being initialized. These attributes can
    be any valid Python object or data type
    :type properties: dict
    """

    components_lines = {}
    compressor_id = 0
    turbine_id = 0
    combustion_chamber_id = 0
    fan_id = 0
    virtual_component_id = 0
    nozzle_id = 0
    mixer_id = 0
    diffuser_id = 0
    afterburner_id = 0
    required_properties = ['LCV']


    def __init__(self, height : float = 0,
                 number_lines : int = 1,
                 properties : dict = {}) -> None:
        
        self.number_lines = number_lines
        self.std = atmosphere.std_ISA()
        self.atm = atmosphere.height_properties(height= height)
        self.adicional_properties = properties

        check_requirements(list_requirements= self.required_properties,
                            attrb_dict= properties)

        for propertie, value in properties.items():
            setattr(self, propertie, value)

        if 'aircraft_vel' in self.adicional_properties.keys():
            M = self.adicional_properties['aircraft_vel']/self.atm['v_sound_h']['value']
            self.adicional_properties['M'] = M
        

    def add_virtual_component(self, T : float = None, P : float = None):
        '''
        funcao para adicionar um componente virtual no motor,
        pode ser usado para acrescentar temperaturas em algumas
        partes do motor sem a necessidade de criar componentes 
        específicos
        '''
        name = f'virtual {self.virtual_component_id}'
        self.components_lines[name] = {'virtual_T' : T, 'virtual_P' : P}


    def add_compressor(self, ef_comp : float = None, rc : float = None, **kwargs):
        """
        The function `add_compressor` adds a compressor to a problem.
        
        :param ef_comp: The parameter "ef_comp" represents the efficiency of the compressor. It is a float
        value that indicates how effectively the compressor can compress the gas. Higher values indicate
        higher efficiency
        :type ef_comp: float
        :param rc: The parameter `rc` stands for the "compression ratio" of the compressor. It is a measure
        of how much the compressor is able to increase the pressure of the gas being compressed. It is
        defined as the ratio of the discharge pressure to the suction pressure
        :type rc: float
        """
        name = f'compressor {self.compressor_id}'
        compressor = structures.Compressor(std_data= self.std, 
                                atm_data= self.atm, 
                                ef_comp= ef_comp, 
                                rc= rc,
                                **kwargs,
                                **self.adicional_properties)
        self.components_lines[name] = compressor
        self.compressor_id += 1

    def add_turbine(self, ef_turb : float = None, rc : float = None, **kwargs):
        """
        The function `add_turbine` adds a turbine to the problem.
        
        :param ef_turb: The `ef_turb` parameter represents the efficiency of the turbine. It is a float
        value that indicates how effectively the turbine converts the energy of the wind into mechanical
        energy
        :type ef_turb: float
        :param rc: The parameter `rc` stands for compression rate. It is a dimensionless parameter that
        represents the increase os pressure given the exit of turbine.
        :type rc: float
        """

        name = f'turbina {self.turbine_id}'
        turbine = structures.Turbine(std_data= self.std, 
                            atm_data= self.atm, 
                            ef_turb= ef_turb, 
                            rc= rc, 
                            **kwargs,
                            **self.adicional_properties)
        self.components_lines[name] = turbine
        self.turbine_id += 1

    def add_combustion_chamber(self, eta : float, pi : float, **kwargs):
        """
        The function `add_combustion_chamber` adds a combustion chamber component to a dictionary of
        components.
        
        :param eta: Eta is the efficiency of the combustion chamber. It represents the ratio of the actual
        energy output of the combustion chamber to the theoretical maximum energy output. It is usually
        expressed as a decimal value between 0 and 1
        :type eta: float
        :param pi: The parameter "pi" represents the pressure ratio across the combustion chamber
        :type pi: float
        """
        name = f'combustion chamber {self.combustion_chamber_id}'
        chamber = structures.Combustion_Chamber(std_data= self.std,
                                                atm_data= self.atm,
                                                eta= eta,
                                                pi= pi,
                                                **self.adicional_properties,
                                                **kwargs)
        self.components_lines[name] = chamber
        self.combustion_chamber_id += 1
    
    def add_diffuser(self, **kwargs):
        """
        The  `add_diffuser`  method adds a diffuser component to a dictionary of components. It takes in keyword arguments to customize the diffuser's properties. 
        The  `name`  variable is generated using a formatted string that includes the diffuser's ID. 
        A  `diffuser`  object is created from the  `structures.Combustion_Chamber`  class, using the standard data ( `self.std` ), atmospheric data ( `self.atm` ), and additional properties ( `self.additional_properties` ) passed as keyword arguments. 
        The  `diffuser`  object is then added to the  `components_lines`  dictionary with the generated  `name`  as the key. 
        """
        name = f'diffuser {self.diffuser_id}'
        diffuser = structures.Diffuser(std_data= self.std,
                                                atm_data= self.atm,
                                                **self.adicional_properties,
                                                **kwargs)
        self.components_lines[name] = diffuser
        self.diffuser_id += 1

    def add_fan(self, eta : float, pi : float, **kwargs):
        """
        The function `add_fan` adds a fan component to a dictionary of components, with specified parameters
        such as efficiency, pressure ratio, split ratio, and bypass ratio.
        
        :param eta: The parameter "eta" represents the efficiency of the fan. It is a float value that
        indicates how effectively the fan converts input power into useful work. Higher values of eta
        indicate higher efficiency
        :type eta: float
        :param pi: The parameter "pi" represents the pressure ratio of the fan. It is a float value that
        indicates the ratio of the outlet pressure to the inlet pressure of the fan
        :type pi: float
        :param split_ratio: The `split_ratio` parameter is a float that represents the ratio of the airflow
        that is bypassed around the fan. It is used to calculate the airflow through the fan and the bypass
        airflow. The default value is 1, which means no bypass airflow, defaults to 1
        :type split_ratio: float (optional)
        :param bypass: The `bypass` parameter is an optional parameter that represents the bypass ratio of
        the fan. It is a float value that determines the ratio of the bypass air flow to the core air flow
        in a turbofan engine. If no value is provided for `bypass`, it defaults to `None`
        :type bypass: float
        """

        name = f'fan {self.fan_id}'
        fan = structures.Fan(std_data= self.std, 
                             atm_data= self.atm,
                            eta = eta, 
                            pi = pi, 
                            **self.adicional_properties,
                            **kwargs)
        
        self.components_lines[name] = fan
        self.fan_id += 1

    def add_nozzle(self, eta_n : float, pi_duct : float, **kwargs):
        """
        The `add_nozzle` function adds a nozzle component to a dictionary of components, with specified
        properties and additional keyword arguments.
        
        :param eta_n: The parameter `eta_n` represents the nozzle efficiency. It is a float value that
        indicates how effectively the nozzle converts the energy of the fluid into useful work
        :type eta_n: float
        :param pi_duct: The parameter `pi_duct` represents the pressure ratio across the duct. It is a float
        value that determines the pressure difference between the inlet and outlet of the duct
        :type pi_duct: float
        """
        name = f'nozzle {self.nozzle_id}'
        nozzle = structures.Nozzle(std_data= self.std, 
                            atm_data= self.atm,
                            eta_n = eta_n, 
                            pi_duct = pi_duct, 
                            **self.adicional_properties,
                            **kwargs)  
        self.components_lines[name] = nozzle
        self.fan_id += 1     

    def add_mixer(self, eta_m : float, pi : float, **kwargs):
        """
        The function `add_mixer` adds a mixer component to a dictionary of components, with specified
        parameters.
        
        :param eta_m: The parameter `eta_m` represents the efficiency of the mixer, which is a float value.
        It is used to calculate the performance of the mixer in terms of how well it mixes the components
        :type eta_m: float
        :param pi: The parameter `pi` represents the pressure ratio of the mixer. It is a float value that
        determines the ratio of the outlet pressure to the inlet pressure of the mixer
        :type pi: float
        """

        name =  f'mixer {self.mixer_id}'
        mixer = structures.Mixer(std_data= self.std,
                                 atm_data= self.atm,
                                 eta_m= eta_m,
                                 pi= pi,
                                **self.adicional_properties,
                                **kwargs)
        self.components_lines[name] = mixer
        self.mixer_id += 1   

    def add_afterburner(self, eta_ab : float, pi : float, **kwargs):
        name =  f'afterburner {self.afterburner_id}'
        afterburner = structures.Afterburner(std_data= self.std,
                                 atm_data= self.atm,
                                 pi= pi,
                                 eta_ab= eta_ab,
                                **self.adicional_properties,
                                **kwargs)
        self.components_lines[name] = afterburner
        self.mixer_id += 1   
    
    def get_motor_attrb(self, **kwargs):
        """
        The function calculates various attributes related to an aircraft's motor efficiency and power and
        updates the historical data with the results.
        """

        for key, value in kwargs.items():
            if type(value) == dict:
                value = value['value']
            setattr(self, key, value)

        self.aircraft_pow = thermo.pow_aircraft(m_dot_air= self.m_dot_air,
                                                m_dot_fuel= self.m_dot_f,
                                                v_jet= self.v_jet,
                                                v_in= self.aircraft_vel)

        self.aircraft_prop_ef = thermo.prop_ef(m_dot_air= self.m_dot_air,
                                               m_dot_fuel= self.m_dot_f,
                                               v_jet= self.v_jet,
                                               v_in= self.aircraft_vel,
                                               LCV= self.LCV)
        
        self.aircraft_therm_ef = thermo.term_ef(m_dot_air= self.m_dot_air,
                                                m_dot_fuel= self.m_dot_f,
                                                v_jet= self.v_jet,
                                                v_in= self.aircraft_vel,
                                                LCV= self.LCV)
        self.aicraft_total_ef = self.aircraft_therm_ef * self.aircraft_prop_ef

        self.F_gr = self.m_dot_air * self.v_jet
        self.F_net = self.m_dot_air * (self.v_jet - self.aircraft_vel)
        self.air_rho = thermo.get_density(P= self.P_final,
                                           T= self.T_final)
        self.A_saida = self.m_dot_air / self.air_rho / self.aircraft_vel

        results = {'aircraft results' : {'prop. efficiency' : self.aircraft_prop_ef,
                                         'therm. efficiency' : self.aircraft_therm_ef,
                                         'total efficiency' : self.aicraft_total_ef,
                                         'power' : self.aircraft_pow,
                                         'net thrust' : self.F_net,
                                         'gross thrust' : self.F_gr,
                                         'Area exit' : self.A_saida}}
        self.hist_data.update(results)
        



    def solve(self, initial_data: dict, save_simulation = True):
        """
        The `solve` function iterates through a dictionary of components and their properties, performing
        calculations and storing the results in a history data dictionary.
        
        :param initial_data: The `initial_data` parameter is a dictionary that contains the initial values
        for the simulation. It is used to set the starting values for the simulation components
        :type initial_data: dict
        """
        index = 0
        self.hist_data = {'start': initial_data}
        properties = {}
        last_id = 'start'
        for id, component in self.components_lines.items():  
            component.id = id
            if 'virtual' in id:
                data = {
                'T_saida' : component['virtual_T'],
                'P_saida' : component['virtual_P'],
            }
                if data['P_saida'] == None:
                    data['P_saida'] = self.hist_data[last_id]['P_saida']

                outputs = data
                
            else:
                if index > 0:
                    component.cp, component.gamma = atmosphere.get_air_cp_gamma(self.hist_data[last_id]['T_saida'])
                    data = {
                        'T_entrada' : self.hist_data[last_id]['T_saida'],
                        'P_entrada' : self.hist_data[last_id]['P_saida'],
                    }
                else:
                    data = self.hist_data['start']

                outputs = component.step(**data, **properties)

            properties['P_final'] = component.P_saida
            properties['T_final'] = component.T_saida

            if 'fan' in id:
                properties['T_bypass'] = component.T_saida
            
            if 'compressor' in id:
                if 'delta_T_compressor' not in properties.keys():
                    properties['delta_T_compressor'] = []
                properties['delta_T_compressor'].append(component.T_saida - component.T_entrada)

            if hasattr(component, 'W'):
                if 'W_requerido' not in properties.keys():
                    properties['W_requerido'] = 0
                properties['W_requerido'] += component.W
            
            if hasattr(component, 'v_jet'):
                properties['v_jet'] = component.v_jet

            if hasattr(component, 'm_dot_f'):
                properties['m_dot_f'] = component.m_dot_f

            self.hist_data[id] = outputs
            last_id = id
            index += 1
        print(properties)
        self.get_motor_attrb(**properties)


        


    def save_simulation(self, filename : str = 'simulation.json'):
        """
        The function `save_simulation` saves the truncated data dictionary to a JSON file with the
        specified filename.
        
        :param filename: The `filename` parameter is a string that specifies the name of the file where the
        simulation data will be saved. By default, it is set to 'simulation.json', defaults to
        simulation.json
        :type filename: str (optional)
        """
        truncated_data = truncate_dict_values(self.hist_data)
        with open(filename, 'w') as file:
            json.dump(truncated_data, file, indent= 4)


if __name__ == "__main__": 
    properties = {'aircraft_vel' : 238.4,
                  'm_dot_air' : 1,
                  'm_dot_fab' : 0.012,
                  'cp_fuel' : 1244,
                  'LCV' : 43*10**6,
                  'bypass' : 0.5,
                  'split_ratio' : 0.5}
    
    teste = Motor(height= 1000,
                  properties= properties)
    
    args = {'T_entrada' : 298,
             'P_entrada' : 1.01*10**5}
    
    #teste.add_diffuser()
    #teste.add_compressor(ef_comp= 0.9, rc= 45)
    #teste.add_virtual_component(T = 1500)
    teste.add_fan(eta = 0.9, pi = 1.6)
    teste.add_compressor(ef_comp= 0.9, rc= 35)
    teste.add_combustion_chamber(eta= 0.91, pi = 0.95, **{'T_saida' : 1700})
    teste.add_afterburner(eta_ab= 0.95, pi = 1.2)
    teste.add_turbine(ef_turb= 0.85)
    teste.add_mixer(eta_m= 0.9, pi = 3)
    teste.add_nozzle(eta_n= 0.96, pi_duct= 0)
    # print(teste.components_lines)
    teste.solve(initial_data= args)
    teste.save_simulation()
    #print(teste.components_lines['turbina 0'].get_work_kg())