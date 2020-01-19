import cantera as ct
import numpy as np


gas = ct.Solution("gri30.xml")

# internal helper function
def get_number_of_elements(gas):

    nO = np.array([gas.n_atoms(k, 'O') for k in range(gas.n_species)])

    if 'C' in gas.element_names:
        nC = np.array([gas.n_atoms(k, 'C') for k in range(gas.n_species)])
    else:
        nC = np.zeros(gas.n_species)

    if 'H' in gas.element_names:
        nH = np.array([gas.n_atoms(k, 'H') for k in range(gas.n_species)])
    else:
        nH = np.zeros(gas.n_species)

    if 'S' in gas.element_names:
        nS = np.array([gas.n_atoms(k, 'S') for k in range(gas.n_species)])
    else:
        nS = np.zeros(gas.n_species)

    return nO, nC, nH, nS

# Bilger mixture fraction
def get_mixture_fraction(gas, fuel, oxidizer):

    original_state = gas.state
    
    nO, nC, nH, nS = get_number_of_elements(gas)
    Z_C = np.dot(gas.Y / gas.molecular_weights, nC)
    Z_O = np.dot(gas.Y / gas.molecular_weights, nO)
    Z_H = np.dot(gas.Y / gas.molecular_weights, nH)
    Z_S = np.dot(gas.Y / gas.molecular_weights, nS)

    gas.TPX = None, None, fuel  
    Z_Cf = np.dot(gas.Y / gas.molecular_weights, nC)
    Z_Of = np.dot(gas.Y / gas.molecular_weights, nO)
    Z_Hf = np.dot(gas.Y / gas.molecular_weights, nH)
    Z_Sf = np.dot(gas.Y / gas.molecular_weights, nS)

    gas.TPX = None, None, ox 
    Z_Co = np.dot(gas.Y / gas.molecular_weights, nC)
    Z_Oo = np.dot(gas.Y / gas.molecular_weights, nO)
    Z_Ho = np.dot(gas.Y / gas.molecular_weights, nH)
    Z_So = np.dot(gas.Y / gas.molecular_weights, nS)
    
    gas.state = original_state

    denominator = 2.0*(Z_Cf-Z_Co) + 0.5*(Z_Hf-Z_Ho) + 2.0*(Z_Sf-Z_So) - (Z_Of-Z_Oo)
    
    if denominator == 0:
        print("Error: fuel and oxidizer have the same composition!")
        return -1 # better way to handle this? 

    return ( 2.0*(Z_C -Z_Co) + 0.5*(Z_H -Z_Ho) + 2.0*(Z_S -Z_So) - (Z_O -Z_Oo) ) / denominator 


def get_equivalence_ratio_new(gas, fuel=None, oxidizer=None):

    # get non-normalzed element mole fractions
    nO, nC, nH, nS = get_number_of_elements(gas)

    # determine local phi from oxygen avaibale in the mixture and oxygen required for fully burning fuel elements
    if fuel == None and oxidizer == None:

        Xm = gas.X # mole fractions in the mixture

        Cm = nC.dot(Xm)
        Om = nO.dot(Xm)
        Hm = nH.dot(Xm)
        Sm = nS.dot(Xm)

        o2_required =  Hm*0.25 + Sm + Cm # amount of oxygen required to oxidize all H, S and C elements
        o2_present = 0.5*Om              # amount of oxygen actually present in the mixture

        if o2_present == 0.0:
            return float('inf') #  mixture is pure fuel

        locPhi = o2_required/o2_present
        
        return locPhi
    
    # determine phi from prescribed fuel and oxidzer composition
    else:

        original_state = gas.state

        Z = get_mixture_fraction(gas, fuel, oxidizer) # mixture fraction
        gas.set_equivalence_ratio(1.0, fuel, oxidizer)
        Z_st = get_mixture_fraction(gas, fuel, oxidizer) # stoichiometric mixture fraction

        gas.state = original_state

        #handle pure oxidizer and pure fuel conditions
        if Z == 0.0: # pure oxidizer
            return 0.0
        if Z == 1.0: # pure fuel
            return float('inf')

        return Z / (1.0 - Z) * (1.0 - Z_st) / Z_st


# case 1) fuel only contains C,H,S and oxidizer only contains O (and N)
fuel = "CH4:1"
ox = "O2:0.21,N2:0.79"
gas.TP = 300, 1e5

print("case 1a) mixture from pure fuel and pure oxidizer at phi=1")
gas.set_equivalence_ratio(1.0, fuel, ox)
print("original method: {:5.4f}".format(gas.get_equivalence_ratio()))
print("local phi:       {:5.4f}".format(get_equivalence_ratio_new(gas)))
print("phi:             {:5.4f}".format(get_equivalence_ratio_new(gas,fuel=fuel, oxidizer=ox)))
print("mix fraction:    {:5.4f}".format(get_mixture_fraction(gas,fuel=fuel, oxidizer=ox)))
print("")

print("case 1b) mixture from pure fuel and pure oxidizer at phi=1 but burnt state")
gas.set_equivalence_ratio(1.0, fuel, ox)
gas.equilibrate('HP')
print("original method: {:5.4f}".format(gas.get_equivalence_ratio()))
print("local phi:       {:5.4f}".format(get_equivalence_ratio_new(gas)))
print("phi:             {:5.4f}".format(get_equivalence_ratio_new(gas,fuel=fuel, oxidizer=ox)))
print("mix fraction:    {:5.4f}".format(get_mixture_fraction(gas,fuel=fuel, oxidizer=ox)))
print("")

print("case 1c) mixture from pure fuel and pure oxidizer at phi=1.3")
gas.TP = 300, 1e5
gas.set_equivalence_ratio(1.3, fuel, ox)
print("original method: {:5.4f}".format(gas.get_equivalence_ratio()))
print("local phi:       {:5.4f}".format(get_equivalence_ratio_new(gas)))
print("phi:             {:5.4f}".format(get_equivalence_ratio_new(gas,fuel=fuel, oxidizer=ox)))
print("mix fraction:    {:5.4f}".format(get_mixture_fraction(gas,fuel=fuel, oxidizer=ox)))

print("case 1d) mixture from pure fuel and pure oxidizer at phi=1.3 but burnt state")
gas.set_equivalence_ratio(1.3, fuel, ox)
gas.equilibrate('HP')
print("original method: {:5.4f}".format(gas.get_equivalence_ratio()))
print("local phi:       {:5.4f}".format(get_equivalence_ratio_new(gas)))
print("phi:             {:5.4f}".format(get_equivalence_ratio_new(gas,fuel=fuel, oxidizer=ox)))
print("mix fraction:    {:5.4f}".format(get_mixture_fraction(gas,fuel=fuel, oxidizer=ox)))
print("")

#case 2) fuel and oxidizer contain fuel species, O2, burnt products and intermediate products
fuel = "CH4:0.2,O2:0.02,N2:0.1,CO:0.05,CO2:0.02"
ox = "O2:0.21,N2:0.79,CO:0.04,CH4:0.01,CO2:0.03"
gas.TP = 300, 1e5

print("case 2a) general mixture at phi=1")
gas.set_equivalence_ratio(1.0, fuel, ox)
print("original method: {:5.4f}".format(gas.get_equivalence_ratio()))
print("local phi:       {:5.4f}".format(get_equivalence_ratio_new(gas)))
print("phi:             {:5.4f}".format(get_equivalence_ratio_new(gas,fuel=fuel, oxidizer=ox)))
print("mix fraction:    {:5.4f}".format(get_mixture_fraction(gas,fuel=fuel, oxidizer=ox)))
print("")

print("case 2b) general mixture at phi=1 but burnt state")
gas.set_equivalence_ratio(1.0, fuel, ox)
gas.equilibrate('HP')
print("original method: {:5.4f}".format(gas.get_equivalence_ratio()))
print("local phi:       {:5.4f}".format(get_equivalence_ratio_new(gas)))
print("phi:             {:5.4f}".format(get_equivalence_ratio_new(gas,fuel=fuel, oxidizer=ox)))
print("mix fraction:    {:5.4f}".format(get_mixture_fraction(gas,fuel=fuel, oxidizer=ox)))
print("")

print("case 2c) general mixture at phi=1.3")
gas.TP = 300, 1e5
gas.set_equivalence_ratio(1.3, fuel, ox)
print("original method: {:5.4f}".format(gas.get_equivalence_ratio()))
print("local phi:       {:5.4f}".format(get_equivalence_ratio_new(gas)))
print("phi:             {:5.4f}".format(get_equivalence_ratio_new(gas,fuel=fuel, oxidizer=ox)))
print("mix fraction:    {:5.4f}".format(get_mixture_fraction(gas,fuel=fuel, oxidizer=ox)))
print("")

print("case 2d) general mixture at phi=1.3 but burnt state")
gas.set_equivalence_ratio(1.3, fuel, ox)
gas.equilibrate('HP')
print("original method: {:5.4f}".format(gas.get_equivalence_ratio()))
print("local phi:       {:5.4f}".format(get_equivalence_ratio_new(gas)))
print("phi:             {:5.4f}".format(get_equivalence_ratio_new(gas,fuel=fuel, oxidizer=ox)))
print("mix fraction:    {:5.4f}".format(get_mixture_fraction(gas,fuel=fuel, oxidizer=ox)))
print("")

# case 3) pure fuel or puel oxidizer condition
fuel = "CH4:1"
ox = "O2:0.21,N2:0.79"

print("case 3a: pure fuel")
gas.X = fuel
print("original method: {:5.4f}".format(gas.get_equivalence_ratio()))
print("local phi:       {:5.4f}".format(get_equivalence_ratio_new(gas)))
print("phi:             {:5.4f}".format(get_equivalence_ratio_new(gas,fuel=fuel, oxidizer=ox)))
print("mix fraction:    {:5.4f}".format(get_mixture_fraction(gas,fuel=fuel, oxidizer=ox)))
print("")

print("case 3b: pure oxidizer")
gas.X = ox
print("original method: {:5.4f}".format(gas.get_equivalence_ratio()))
print("local phi:       {:5.4f}".format(get_equivalence_ratio_new(gas)))
print("phi:             {:5.4f}".format(get_equivalence_ratio_new(gas,fuel=fuel, oxidizer=ox)))
print("mix fraction:    {:5.4f}".format(get_mixture_fraction(gas,fuel=fuel, oxidizer=ox)))
print("\n")

fuel = "CH4:0.2,O2:0.02,N2:0.1,CO:0.05,CO2:0.02"
ox = "O2:0.21,N2:0.79,CO:0.04,CH4:0.01,CO2:0.03"

print("case 3c: pure fuel general mixture")
gas.X = fuel
print("original method: {:5.4f}".format(gas.get_equivalence_ratio()))
print("local phi:       {:5.4f}".format(get_equivalence_ratio_new(gas)))
print("phi:             {:5.4f}".format(get_equivalence_ratio_new(gas,fuel=fuel, oxidizer=ox)))
print("mix fraction:    {:5.4f}".format(get_mixture_fraction(gas,fuel=fuel, oxidizer=ox)))
print("")

print("case 3d: pure oxidizer general mixture")
gas.X = ox
print("original method: {:5.4f}".format(gas.get_equivalence_ratio()))
print("local phi:       {:5.4f}".format(get_equivalence_ratio_new(gas)))
print("phi:             {:5.4f}".format(get_equivalence_ratio_new(gas,fuel=fuel, oxidizer=ox)))
print("mix fraction:    {:5.4f}".format(get_mixture_fraction(gas,fuel=fuel, oxidizer=ox)))
print("\n")


