"""
Polycrude manufacturing via plastic waste hydrocracking

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 06/13/2025
"""
import qsdsan as qs
import biosteam as bst
from biosteam.units.design_tools import CEPCI_by_year
from biosteam.units.decorators import cost
from qsdsan.utils import auom

qs.CEPCI=567.5  #CEPCI 2017 IN QSDSAN
class Electrolysis(qs.SanUnit):
    """ Create a electrolysis PEM plant for H2 distributed production 
    using the H2A models reported by NREL [https://www2.nrel.gov/hydrogen/h2a-production-models]
    
    Parameters
    ----------
    ID : str
        ID of the electrolysis system.
    ins : Iterable(stream)
        water_in.
    outs : Iterable(stream)
        hydrogen_out, oxygen_out.
    electricity_usage: float
        Electricity consumption per mass of hydrogen produced, [kWh/kg hydrogen].
    T: float, optional
        Operating temperature of the electrolytic cell , [K] (the default is 353.15 [1]).
    P : float, optional
        Hydrogen outlet pressure, [Pa] (the default is 3102642 [1]).
    stack_oversize : str, optional
        Stack oversize due to degradation (the default is 0.13 [1]).
    

    References
    ---------- 
    [1] James, B.; Colella, W.; Moton, J.; Saur, G.; Ramsden, T. PEM Electrolysis H2A Production
        Case Study Documentation; NREL/TP--5400-61387, 1214980; 2013; p NREL/TP--5400-61387, 1214980. 
        https://doi.org/10.2172/1214980

    """
    _N_ins=1  #Water 'kg/hr'
    _N_outs=2  #H2_makeup to use in the plant

    _units={'H2_produced':'kg/h', 'Electricity_consumption':'kw'}  

    def __init__(self, ID, ins=None, outs=(), electricity_usage=None, thermo=None,  T=80+273.15, P=450*6894.76, stack_oversize= 0.13): #Operating temperature and outlet pressure, electricity usage in kWh/kgH2
        super().__init__(ID, ins, outs, thermo)                                                                                       
        self.T=T
        self.P=P
        self.electricity_usage=electricity_usage
        self.stack_oversize=stack_oversize
      
        
    def _run(self):

        feed, = self.ins           #comma to unpack the a single stream  
        makeup_hydrogen, oxygen = self.outs    

        #makeup_hydrogen.imass['H2']=feed.F_mass/3.5/3.79  
        feed.imass['H2O'] = makeup_hydrogen.F_mass*3.79/3.5     #3.5 kgH2 per gal H2O and 3.79 water density kg/gal
        oxygen.imass['O2'] =7.94*makeup_hydrogen.F_mass       #production of oxigen kgO2/kg H2
        makeup_hydrogen.phase='g'
        oxygen.phase='g'
        feed.phase='l'
        makeup_hydrogen.T=self.T
        makeup_hydrogen.P=self.P
        oxygen.T=self.T
        oxygen.P=101324

    def _design(self):
        Electricity_consumption= self.outs[0].F_mass*self.electricity_usage  #55.8 kWh/kg H2  Fmass return kg/h- total units kWh  
        
        self.add_power_utility(Electricity_consumption)    #add the utilty requirements

        self.design_results['H2_produced'] =self.outs[0].F_mass   
        self.design_results['Electricity_consumption']=Electricity_consumption                            #design resutls is available for reports and integration with cost estimation
        
    def _cost(self):
        stack_input_power_peak=(self.electricity_usage-5.4)*self.outs[0].F_mass*(1+self.stack_oversize)/1000  #stack electrical usage is the total electricity usage-BoP electricity usage (5.4) in kWh/kg/H2, peak is the max production in MW -not affected by degradation rate
        stack_cost=1.3*1000/(2*1.9)      #1.3 $/cm^2, 2 A/cm^2 and 1.9 V/cell, $/kW
        mechanical_BOP=286*self.outs[0].F_mass*(1+self.stack_oversize)*24/stack_input_power_peak/1000    #286 $/kgH2/day
        electronical_BOP=121         #$/KW

        purchase_cost=(stack_cost+mechanical_BOP+electronical_BOP)*stack_input_power_peak*1000*(qs.CEPCI/541.7)  #    Reference year 2016

        self.baseline_purchase_costs['Electrolysis'] = purchase_cost   #store it without accounting material costs
        self.F_BM['Electrolysis']=1.12