"""
Polycrude manufacturing via plastic waste hydrocracking

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 03/10/2025
"""

import qsdsan as qs
import numpy as np


qs.CEPCI=567.5  #CEPCI 2017 IN QSDSAN


class Conveyor(qs.SanUnit):
    """ Creates a belt conveyor to transport solids

    These conveyor calculate the inlet stream based on the plant capacity
    parameter, `feedPW`.

    Parameters
    ----------
    ID : str
        ID of the conveyor.
    ins : Iterable(stream)
        Inlet.
    outs: Iterable(stream)
        Outlet.
    T : float, optional
        Temperature, [K] (the default is 298.15 K, which correspond to ambient conditions).
    P : float, optional
        Pressure, [Pa] (the default is 101325 Pa, which correspond to ambient conditions).
    feedPW : float, optional
        Plastic waste feedstock per hour or plant capacity, [kg/h] (the default is 10416.7 kg/h,
        which is the suggested minimum size for economic feasibility in plastic recycling plants [1]).
    purity_factor : float
        Mass fraction of LDPE content in the plastic waste feedstock

    Notes
    -----
    The conveyor design, power consumption and cost equations are based 
    on [2], utilizing typical values for volumetric capacity, belt velocity and belt 
    width as referenced in [2]. Costs are updated to 2023 using CEPCI.

    References
    ----------
    [1] Yadav, G.; Singh, A.; Dutta, A.; Uekert, T.; DesVeaux, J. S.; Nicholson, S. R.; 
        Tan, E. C. D.; Mukarakate, C.; Schaidle, J. A.; Wrasman, C. J.; Carpenter, A. C.; Baldwin, R. M.; 
        Román-Leshkov, Y.; Beckham, G. T. Techno-Economic Analysis and Life Cycle Assessment for Catalytic 
        Fast Pyrolysis of Mixed Plastic Waste. Energy Environ. Sci. 2023, 16 (9), 3638-3653. 
        https://doi.org/10.1039/D3EE00749A.
    
    [2] Seider, W. D.; Lewin, D. R.; Seader, J. D.; Widagdo, S.; Gani, R.;
        Ng, M. K. Purchase Costs of Other Chemical Processing Equipment. In Product
        and Process Design Principles; Wiley, 2017; pp 477-484.

    """

    _units={'Conveyor_feed':'tonnes/h', 'Conveyor_volumetricflow':'ft^3/h',
            'Power':'kw-h', 'Width': 'in', 'Length':'ft'}   

    def __init__(self, ID, ins=None, outs=(), thermo=None, T=298.15, P=101325,  #not LDPE
                feedPW=10416.7,  
                purity_factor=None,  
                 ):
        super().__init__(ID, ins, outs,thermo)  
        self.T=T
        self.P=P  
        self.feedPW=feedPW
        self.purity_factor=purity_factor

    _N_ins=1
    _N_outs=1

    def _run(self):
        #inlet flow rates with composition                                                        
        LDPEin = self.feedPW*(self.purity_factor)
        self.ins[0].imass['LDPE']= LDPEin  #kg/h

        Impuritiesin = self.feedPW*(1-self.purity_factor)
        self.ins[0].imass['impurities']=Impuritiesin   #kg/h

        #No change in flows in the outlets

        self.outs[0].imass['LDPE']=self.ins[0].imass['LDPE'] 
        self.outs[0].imass['impurities']=self.ins[0].imass['impurities']  

    def _design(self):

        mass_flow=self.ins[0].get_total_flow('kg/h')
        density = self.ins[0].rho  #kg/m^3
        volumetric_flow= (mass_flow/density)*35.3147    #ft^3/h

        if density ==0 or density is None:
            raise ValueError('Density is not defined, check input stream properties')
       
        volumetric_capacity = 660    #ft^3/h
        velocity=100  # ft/s, average velocity 
        area= volumetric_flow*60/velocity
        base_width = 14
        
        width=base_width*max(1,volumetric_flow/volumetric_capacity )
        length= area/(width/12)   #ft
        mass_flow_lb_s= self.ins[0].get_total_flow('lb/s')   

        Conveyor_power = 0.00058 * (mass_flow_lb_s**0.82)*length*0.7457*24   #kw-h
        
        self.design_results['Conveyor_feed'] =self.ins[0].get_total_flow('tonnes/h')   
        self.design_results['Conveyor_volumetricflow'] =volumetric_flow       
        self.design_results['Power']=Conveyor_power 
        self.design_results['Width']=width   
        self.design_results['Length']=length                               
        self.add_power_utility(Conveyor_power)    

    def _cost(self):
        W= self.design_results['Width']  
        L=self.design_results['Length']
        V=self.design_results['Conveyor_volumetricflow']
        
        purchase_cost=(24.4*W*L + 1094*V**0.22)*(qs.CEPCI/567)  
        

        self.baseline_purchase_costs['Conveyor'] = purchase_cost   
        self.F_D['Conveyor']=self.F_P['Conveyor']=self.F_M['Conveyor']=1  
        self.F_BM['Conveyor']=1.61
