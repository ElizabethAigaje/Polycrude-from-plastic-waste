import qsdsan as qs

__all__ = ('create_components', )


def create_components():
#  Setup inlet flows
  
    #heptacontane_chem = qs.Chemical('heptacontane')
    #LDPE = qs.Component.from_chemical(ID='LDPE', chemical=heptacontane_chem, particle_size='Soluble', degradability='Undegradable', organic=True)
    LDPE = qs.Component('LDPE', formula='C2500H5002', measured_as='COD', phase='s', particle_size='Soluble', degradability='Undegradable', organic=True)
    heptacontane_chem = qs.Chemical('heptacontane')
    #LDPE.CAS ==heptacontane_chem.CAS
    LDPE.Tm=heptacontane_chem.Tm
    LDPE.Tb=heptacontane_chem.Tb
    LDPE.Pc=heptacontane_chem.Pc
    LDPE.Tc=heptacontane_chem.Tc
    LDPE.Vc=heptacontane_chem.Vc
    LDPE.Hfus=heptacontane_chem.Hfus
    LDPE.Hf=heptacontane_chem.Hf
    LDPE.copy_models_from(heptacontane_chem)
    # LDPE.Hvap.add_method(heptacontane_chem.Hvap)
    # LDPE.Cp.add_method(heptacontane_chem.Cp)
    # LDPE.Cn.add_model(heptacontane_chem.Cn)
    # LDPE.mu.add_model(0.1)
    # LDPE.V.add_method(0.03804)   #molar volume in solid
    #LDPE.V.l.add_method(0.03804)    #molar volume in liquid
    #LDPE.rho.add_method(heptacontane_chem.rho)
#Polycrude
# Create the surrogate chemical first (ensure it has all needed properties)
    pentatriacontane_chem = qs.Chemical('pentatriacontane')
    polycrude = qs.Component.from_chemical(ID='polycrude', chemical=pentatriacontane_chem, particle_size='Soluble', degradability='Undegradable', organic=True)
    polycrude.S0=379.09
    polycrude.Sfus=0

#Impurities
# Create the surrogate chemical first (ensure it has all needed properties)
    C_chem = qs.Chemical('C')
    impurities = qs.Component.from_chemical(ID='impurities', chemical=C_chem, particle_size='Particulate', degradability='Undegradable', organic=True)


#Light Hydrocracbons
    C1H4_chem=qs.Chemical('C1H4')
    C1H4 = qs.Component.from_chemical(ID='C1H4', chemical=C1H4_chem, particle_size='Soluble', degradability='Undegradable', organic=True)

    C2H6_chem=qs.Chemical('ethane')
    C2H6 = qs.Component.from_chemical(ID='C2H6', chemical=C2H6_chem, particle_size='Soluble', degradability='Undegradable', organic=True)
    # Psat, Tb and Hvap

    C3H8_chem=qs.Chemical('C3H8')
    C3H8 = qs.Component.from_chemical(ID='C3H8', chemical=C3H8_chem, particle_size='Soluble', degradability='Undegradable', organic=True)

    C4H10_chem=qs.Chemical('C4H10')
    C4H10 = qs.Component.from_chemical(ID='C4H10', chemical=C4H10_chem, particle_size='Soluble', degradability='Undegradable', organic=True)

    C5H12_chem=qs.Chemical('C5H12')
    C5H12 = qs.Component.from_chemical(ID='C5H12', chemical=C5H12_chem, particle_size='Soluble', degradability='Undegradable', organic=True)

    C6H14_chem=qs.Chemical('C6H14')
    C6H14 = qs.Component.from_chemical(ID='C6H14', chemical=C6H14_chem, particle_size='Soluble', degradability='Undegradable', organic=True)

#Solvent range hydrocarbons

    C7H16_chem=qs.Chemical('C7H16')
    C7H16 = qs.Component.from_chemical(ID='C7H16', chemical=C7H16_chem, particle_size='Soluble', degradability='Undegradable', organic=True)

    C8H18_chem=qs.Chemical('C8H18')
    C8H18 = qs.Component.from_chemical(ID='C8H18', chemical=C8H18_chem, particle_size='Soluble', degradability='Undegradable', organic=True)

    C9H20_chem=qs.Chemical('C9H20')
    C9H20 = qs.Component.from_chemical(ID='C9H20', chemical=C8H18_chem, particle_size='Soluble', degradability='Undegradable', organic=True)

    C10H22_chem=qs.Chemical('C10H22')
    C10H22 = qs.Component.from_chemical(ID='C10H22', chemical=C10H22_chem, particle_size='Soluble', degradability='Undegradable', organic=True)

    C11H24_chem=qs.Chemical('C11H24')
    C11H24 = qs.Component.from_chemical(ID='C11H24', chemical=C11H24_chem, particle_size='Soluble', degradability='Undegradable', organic=True)

    C12H26_chem=qs.Chemical('C12H26')
    C12H26 = qs.Component.from_chemical(ID='C12H26', chemical=C12H26_chem, particle_size='Soluble', degradability='Undegradable', organic=True)
    C12H26.Psat.method='AMBROSE_WALTON'   # to have access a higher operating pressures
#Light gases
    H2_chem=qs.Chemical('H2')
    H2 = qs.Component.from_chemical(ID='H2', chemical=H2_chem, particle_size='Soluble', degradability='Undegradable', organic=False)
    H2._locked_state='g'
    
#Catalyst
    catalyst_chem=qs.Chemical('Al2O5Si')
    catalyst = qs.Component.from_chemical(ID='catalyst', chemical=catalyst_chem, particle_size='Particulate', degradability='Undegradable', organic=False)  
      #molar volume m^3/mol with density 1249.44 kg/m3 and MW 162 gr/mol
    catalyst.V.s.add_method(0.0001296)
    catalyst._locked_state='s'
    
    
#For H2 production
    H2O_chem=qs.Chemical('H2O')
    H2O = qs.Component.from_chemical(ID='H2O', chemical=H2O_chem, particle_size='Soluble', degradability='Undegradable', organic=False)


    O2_chem=qs.Chemical('O2')
    O2 = qs.Component.from_chemical(ID='O2', chemical=O2_chem, particle_size='Soluble', degradability='Undegradable', organic=False)
   

    cmps = qs.Components([LDPE, polycrude, impurities, catalyst, C1H4, C2H6, C3H8, C4H10, C5H12, C6H14, C7H16, C8H18, C9H20, C10H22, C11H24, C12H26, H2, C4H10, H2O, O2]) # *cmps_default
    qs.set_thermo(cmps)

    return cmps
    

cmps=create_components()
print(cmps)
#H2_component=cmps['catalyst']
#x=H2_component._locked_state
#print(dir(H2_component))
#print(x)