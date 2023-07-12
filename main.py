import sys
import include.core as core

if __name__ == "__main__": 
    properties = {'aircraft_vel' : 238.4,
                  'm_dot_air' : 13.61,
                  'cp_fuel' : 1244,
                  'LCV' : 43.96*10**6,
                  'bypass' : 0.5,
                  'split_ratio' : 0.5}
    
    teste = core.Motor(height= 1000,
                  properties= properties)
    
    args = {'T_entrada' : 298,
             'P_entrada' : 1.01*10**5}
    
    teste.add_diffuser()
    #teste.add_compressor(ef_comp= 0.9, rc= 45)
    #teste.add_virtual_component(T = 1500)
    #teste.add_fan(eta = 0.9, pi = 1.6)
    teste.add_compressor(ef_comp= 0.88, rc= 6.5)
    teste.add_combustion_chamber(eta= 0.91, pi = 0.95, **{'T_saida' : 1389})
    teste.add_turbine(ef_turb= 0.85)
    #teste.add_mixer(eta_m= 0.9, pi = 3)
    teste.add_nozzle(eta_n= 0.96, pi_duct= 0)
    # print(teste.components_lines)
    teste.solve(initial_data= args)
    teste.save_simulation()
    #print(teste.components_lines['turbina 0'].get_work_kg())

