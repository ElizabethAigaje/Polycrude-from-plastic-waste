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

class Grinder(qs.SanUnit):
    """ Creates a grinder that reduce LDPE into flakes

    Parameters
    ----------
    ID : str
        ID of the grinder.
    ins : Iterable(stream)
        Inlet.
    outs: Iterable(stream)
        Outlet.
    power_requirements: float
        Energy required for size reduction, [kW-h/ton].
    T : float, optional
        Temperature, [K] (the default is 298.15 K, which correspond to ambient conditions).
    P : float, optional
        Pressure, [Pa] (the default is 101325 Pa, which correspond to ambient conditions).

    Notes
    -----
    The grinder cost equation is based on [1] and updated to 2023 cost with CEPCI.

    References
    ----------
    [1] Seider, W. D.; Lewin, D. R.; Seader, J. D.; Widagdo, S.; Gani, R.;
        Ng, M. K. Purchase Costs of Other Chemical Processing Equipment. In Product
        and Process Design Principles; Wiley, 2017; pp 484.

    """

    _N_ins=1
    _N_outs=1
    _units={'Grinder_feed':'tonnes/h', 'Power':'kw'}   #units for the design results

    def __init__(self, ID, ins=(), outs=(), power_requirements=None, thermo=None,  T=298.15, P=101325):
        super().__init__(ID, ins, outs, thermo)    #inherinting the funtionalities of Units
        self.power_requirements=power_requirements   
        
    def _run(self):
        """ Run method to ensure output mirrors the input stream with no change in composition. """
        self.outs[0].copy_like(self.ins[0])  #the output mirrors the input stream. No change in composition, flow or properties.
        
    def _design(self):
        """ Calculate power requirements and store in design results. """
        Grinder_power = self.power_requirements * self.ins[0].get_total_flow('tonnes/h')  
        self.design_results['Grinder_feed'] =self.ins[0].get_total_flow('tonnes/h')        #design results is a built in attribute/dictionary that stores the design parameters. 
        self.design_results['Power']=Grinder_power                               
        self.add_power_utility(Grinder_power)    #add the utilty requirements

    def _cost(self):
        """ Calculate purchase cost based on feed. """
        S=self.design_results['Grinder_feed']
        if S<=200:      #maximum size of a grinder [1]
            purchase_cost=(4310*S**0.78)*(qs.CEPCI/567)   
        else:
            number_grinder= np.ceil(S/200)
            purchase_cost=number_grinder*(4310*200**0.78)*(qs.CEPCI/567) 
        
        self.baseline_purchase_costs['Grinder'] = purchase_cost  
        self.F_D['Grinder']=self.F_P['Grinder']=self.F_M['Grinder']=1  
        self.F_BM['Grinder']=2

