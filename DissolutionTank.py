"""
Polycrude manufacturing via plastic waste hydrocracking

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 03/17/2025
"""

import qsdsan as qs
import biosteam as bst
from biosteam.units.design_tools import CEPCI_by_year
from biosteam.units.decorators import cost
from qsdsan.utils import auom

class Dissolution_tank(qs.sanunits.MixTank): 
    """ Creates a agitated tank that dissolves LDPE in dodecane for polycrude production.

    Parameters
    ----------
    ID : str
        ID of the dissolution tank.
    ins : Iterable(stream)
        feed, makeup_solvent, recycled_solvent.
    outs: Iterable(stream)
        dissolution, undissolved.
    feedPW : float, optional
        Plastic waste feedstock per hour or plant capacity, [kg/h] (the default is 10416.7 kg/h,
        recommended for economic feasibility in plastic recycling plants [1]).
    LDPE_concentration : float
        Mass fraction of LDPE in solvent.
    max_solubility : float, optional
        Maximum solubility of LDPE in solvent at the dissolution temperature (the default is 0.33,
        which is the maximum solubility of LDPE in dodecane at 120 C [2]).
    tau : float
        Dissolution time, [hr].
    kW_m3 : float
        Power required for mixing (the default is 1.5, which correspond to severe mixing in slurry suspensions [3])

    Notes
    -----
    The total solvent requirement is calculated from the `feedPW` and the `LDPE_concentration`.
    The makeup solvent is calculated as the total amount required minus the amount of recycled solvent.

    References
    ---------
    [1] Yadav, G.; Singh, A.; Dutta, A.; Uekert, T.; DesVeaux, J. S.; Nicholson, S. R.; 
        Tan, E. C. D.; Mukarakate, C.; Schaidle, J. A.; Wrasman, C. J.; Carpenter, A. C.; Baldwin, R. M.; 
        Román-Leshkov, Y.; Beckham, G. T. Techno-Economic Analysis and Life Cycle Assessment for Catalytic 
        Fast Pyrolysis of Mixed Plastic Waste. Energy Environ. Sci. 2023, 16 (9), 3638-3653. 
        https://doi.org/10.1039/D3EE00749A.
    [2] Zhou, P.; Yu, J.; Sánchez-Rivera, K. L.; Huber, G. W.; Van Lehn, R. C. 
        Large-Scale Computational Polymer Solubility Predictions and Applications to Dissolution-Based Plastic Recycling.
        Green Chem. 2023, 25 (11), 4402–4414. https://doi.org/10.1039/D3GC00404J.   
    [3] Zolghadr, A.; Foroozandehfar, A.; Kulas, D. G.; Shonnard, D. Study of the Viscosity and Thermal Characteristics 
        of Polyolefins/Solvent Mixtures: Applications for Plastic Pyrolysis. ACS Omega 2021, 6 (48), 32832–32840.
        https://doi.org/10.1021/acsomega.1c04809.
    """

    _N_ins=3
    _N_outs=2


    def __init__(self, ID, ins=(), outs=(), Thermo=None, 
                 feedPW= 10416.7,  #kg/h also capacity
                 LDPE_concentration=None, max_solubility=0.33, tau=(), kW_per_m3=1.5):   #max solubility at 120 C for LPDE in solvent, power required to mix slurry suspensions 
        super().__init__(ID, ins, outs, Thermo) 
        self.feedPW=feedPW
        self.LDPE_concentration=LDPE_concentration
        self.max_solubility=max_solubility
        self.tau=tau
        self.kW_per_m3=kW_per_m3

    def _run(self):
        """Calculate material flow distribution and phase conditions."""
        feed, makeup_solvent, recycled_solvent = self.ins           
        dissolution, undissolved = self.outs               
        
        # Calculate the amount of solvent and makeup solvent required 

        solvent_needed=self.feedPW/self.LDPE_concentration  #kg/h 
        makeup_solvent.imass['C12H26']=solvent_needed-recycled_solvent.imass['C12H26']    

        # Calculate the amount of LDPE dissolved in the solvent 
        # and the impurites (and LDPE not dissolved) to be removed 
        
        ldpe_mass=feed.imass['LDPE']
        undissolved_mass=feed.imass['impurities']   # It is assumed that the solvent does not dissolve any impurities
        dodecane_mass=makeup_solvent.imass['C12H26']+recycled_solvent.imass['C12H26']

        max_ldpe_dissolved = self.max_solubility*dodecane_mass


        if ldpe_mass <= max_ldpe_dissolved:
            #Fully dissolve LDPE
            dissolution.imass['LDPE']=ldpe_mass
            dissolution.imass['C12H26']=dodecane_mass
            dissolution.phase ='l'
            undissolved.imass['impurities']=undissolved_mass
            undissolved.phase='s'
            dissolution.T= recycled_solvent.T 
            undissolved.T= 25+273.15
        else:
            dissolution.imass['LDPE']=max_ldpe_dissolved
            dissolution.imass['C12H26']=dodecane_mass
            dissolution.phase ='l'
            dissolution.T= recycled_solvent.T 

            undissolved.imass['impurities']=undissolved_mass
            undissolved.imass['LDPE']=ldpe_mass-max_ldpe_dissolved
            undissolved.phase='s'
            undissolved.T=25+273.15


    def _design(self):
        """Inherits design calculation from MixTank."""
        qs.sanunits.MixTank._design(self)

    def _cost(self):
        """Inherits cost calculation from MixTank."""
        qs.sanunits.MixTank._cost(self)
       
         
         



   