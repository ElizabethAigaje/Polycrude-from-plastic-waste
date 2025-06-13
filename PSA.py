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


_kg_to_lb_= auom('kg').conversion_factor('lb')

@cost(basis='PSA_hydrogen', ID='PSA', units='lb/h', cost=345557.72, S=25444, CE=CEPCI_by_year[2017], n=0.6, BM=2.3)  #accounted for adsorption and desorption units
        
class PSA(qs.SanUnit):
    """ Creates a pressure swing adsorption (PSA) unit to separate hdrogen from light hydrocarbons
    
    Parameters
    ----------
    ID : str
        ID of the PSA
    ins : Iterable(stream)
        gases_in.
    outs : Iterable(stream)
        light_gases, hydrogen_recovered.
    T : float, optional
        Operating PSA temperature, [K] (the default is 323.15 [1,2]).
    P : float, optional
        Operating PSA pressure, [Pa] (the default is 1500000 [1,2]).
    efficiency : float, optional
        Hydrogen recovery at PSA (the default is 0.85 [1,2]).
    
    Notes
    -----
    PSA cost is based on reference [3]. This cost account for 2 PSA (accounting for adsorption and desorption)
    We assume all the other light components (except hydrogen) completely adsorbed.

    References
    ---------- 
    [1] Dutta, A.; Sahir, A.; Tan, E. Process Design and Economics for the Conversion of Lignocellulosic Biomass 
        to Hydrocarbon Fuels: Thermochemical Research Pathways with In Situ and Ex Situ Upgrading of Fast Pyrolysis Vapors.
    [2] P. Cavaliere, Water Electrolysis for Hydrogen Production. Cham: Springer International Publishing, 
        2023. doi: 10.1007/978-3-031-37780-8
    [3] Cappello, V.; Sun, P.; Zang, G.; Kumar, S.; Hackler, R.; Delgado, H. E.; Elgowainy, A.; Delferro, M.; Krause, T. 
        Conversion of Plastic Waste into High-Value Lubricants: Techno-Economic Analysis and Life Cycle Assessment. 
        Green Chem. 2022, 24 (16), 6306–6318
    """

    _N_ins=1
    _N_outs=2

    _units={'PSA_hydrogen':'lb/h' }  

    def __init__(self, ID, ins=(), outs=(), Thermo=None, T=50+273.15, P= 15*100000, efficiency=0.85):   
        super().__init__(ID, ins, outs, Thermo)    
        self.T=T
        self.P=P
        self.efficiency=efficiency

         
    def _run(self):
        # Assumption that other components are adsorb and desorb

        feed, = self.ins           #comma to unpack the a single stream  
        effluent, hydrogen = self.outs           #tuple unpacking to check that the #of items in self.outs matches the number of variables on the left side     
        effluent.P=151685                    #Desorption at pressure near to atmospheric pressure, NREL at 22 psig
        effluent.T=self.T
        hydrogen.P=self.P*0.96
        hydrogen.T=self.T

        
        #effluent.copy_like(self.ins)  
        #hydrogen.copy_like(feed)   #1st we pass the properties--As in bst column we are intersted in Hydrogen

        hydrogen.imass['H2']= self.efficiency*feed.imass['H2']  #kg/hr   
        hydrogen.phase='g'
        
        hydrogen_notrecovered = (1-self.efficiency)*feed.imass['H2']
        effluent.imass['H2']= hydrogen_notrecovered  #just update the hydrogen flow

        for component in feed.chemicals.IDs:
            if component != 'H2':
                effluent.imass[component]=feed.imass[component]
        
        #effluent.imass['C12H26'] = feed.imass['C12H26']  #copy the inlet
        #effluent.imass['C3H8'] = feed.imass['C3H8']  #copy the inlet
        #effluent.imass['C4H10'] = feed.imass['C4H10']  #copy the inlet

        effluent.phase='g'

        #effluent.vle(T=effluent.T, P=effluent.P)

    def _design(self):
         
         self.design_results['PSA_hydrogen'] = self.outs[1].F_mass*_kg_to_lb_   #H2 produced is the base for the cost in lb/h
                                   
         



   

    