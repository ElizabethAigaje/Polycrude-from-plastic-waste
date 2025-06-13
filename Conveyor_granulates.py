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

# Inheriting from unit requires initialization, run (what happens in the process), design and costs methods

class Conveyor_granulates(qs.SanUnit):
    """ Creates a belt conveyor to transport solids

    Parameters
    ----------
    ID : str
        ID of the conveyor.
    ins : Iterable(stream)
        Inlet.
    outs: Iterable(stream)
        Outlet.
    T : float, optional
        Temperature (the default is 298.15 K, which correspond to ambient conditions), [K].
    P : float, optional
        Pressure (the default is 101325 Pa, which correspond to ambient conditions), [Pa].
    
    Notes
    -----
    The conveyor design, power consumption and cost equations are based 
    on [1], utilizing typical values for volumetric capacity, belt velocity and belt 
    width as referenced in [1].Costs are updated to 2023 using CEPCI.

    References
    ----------
    [1] Seider, W. D.; Lewin, D. R.; Seader, J. D.; Widagdo, S.; Gani, R.;
        Ng, M. K. Purchase Costs of Other Chemical Processing Equipment. In Product
        and Process Design Principles; Wiley, 2017; pp 477-484.

    """

    #Information in the design results
    _units={'Conveyor_feed':'tonnes/h', 'Conveyor_volumetricflow':'ft^3/h',
            'Power':'kw-h', 'Width': 'in', 'Length':'ft'}      

    def __init__(self, ID, ins=(), outs=(), thermo=None, T=298.15, P=101325):
        super().__init__(ID, ins, outs,thermo)  
        self.T=T
        self.P=P  


    _N_ins=1
    _N_outs=1

    def _run(self):

        self.outs[0].copy_like(self.ins[0])  

    def _design(self):

        mass_flow=self.ins[0].get_total_flow('kg/h')
        density = self.ins[0].rho  #kg/m^3
        volumetric_flow= (mass_flow/density)*35.3147    #ft^3/h

        if density ==0 or density is None:
            raise ValueError('Density is not defined, check input stream properties')
       
        volumetric_capacity = 660    #ft^3/h
        velocity=100  # ft/s  average velocity 
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
        
        purchase_cost=(24.4*W*L + 1094*V**0.22)*(qs.CEPCI/567)  #equation cost for conveyor and hopper from [1]
        

        self.baseline_purchase_costs['Conveyor'] = purchase_cost   
        self.F_D['Conveyor']=self.F_P['Conveyor']=self.F_M['Conveyor']=1  #design, pressure, material factors
        self.F_BM['Conveyor']=1.61
