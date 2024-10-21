import numpy as np
import pandas as pd
import Thermobar as pt
import julia
# import time
import pyMelt as m
# from scipy.optimize import fsolve
from julia.api import Julia
jl = Julia(compiled_modules=False)
from julia import MAGEMinCalc

def equilibrate_multi(P_bar = None, T_C = None, comp = None):
    Results = pd.DataFrame()

    comp['O'] = comp['Fe3Fet_Liq']*(((159.59/2)/71.844)*comp['FeOt_Liq'] - comp['FeOt_Liq'])

    bulk = comp[['SiO2_Liq', 'Al2O3_Liq', 'CaO_Liq', 'MgO_Liq', 'FeOt_Liq', 'K2O_Liq', 'Na2O_Liq', 'TiO2_Liq', 'O', 'Cr2O3_Liq', 'H2O_Liq']].astype(float).values

    Results = MAGEMinCalc.equilibrate(bulk, P_bar/1000.0, T_C)
    return Results

def equilibrate(P_bar = None, T_C = None, comp = None, fO2_buffer = None, fO2_offset = None):
    Results = pd.DataFrame()

    if comp is None:
        raise Exception("No composition specified")
        
    else:
        if fO2_buffer is not None:
            new_bulk = pd.DataFrame(comp, index = [0])
            new_bulk['Sample_ID_Liq'] = 0
            if fO2_offset is not None:
                new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = T_C+273.15, P_kbar = P_bar/1000, fo2 = fO2_buffer, fo2_offset = fO2_offset, model = "Kress1991", renorm = False)
            else:
                new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = T_C+273.15, P_kbar = P_bar/1000, fo2 = fO2_buffer, model = "Kress1991", renorm = False)

            comp = new[list(comp.keys())].loc[0].to_dict()

        bulk = [comp['SiO2_Liq'], comp['Al2O3_Liq'], comp['CaO_Liq'], comp['MgO_Liq'], comp['FeOt_Liq'], comp['K2O_Liq'], comp['Na2O_Liq'], comp['TiO2_Liq'], comp['Fe3Fet_Liq']*(((159.59/2)/71.844)*comp['FeOt_Liq'] - comp['FeOt_Liq']), comp['Cr2O3_Liq'], comp['H2O_Liq']]

    if type(T_C) == int:
        T_C = float(T_C)
    if type(P_bar) == int:
        P_bar = float(P_bar)

    try:
        Ret = MAGEMinCalc.PT_minimisation(P_bar/1000, T_C, bulk)
    except:
        return Results

    PhaseList = Ret['Phase']
    PhaseList = list(np.unique(PhaseList))

    Results = pd.DataFrame(columns = ['T_C', 'P_bar'] + PhaseList, data = np.zeros((1, len(['T_C', 'P_bar'] + PhaseList))))
    Results['T_C'] = T_C
    Results['P_bar'] = P_bar
    for p in PhaseList:
        Results[p] = 'Y'

    return Results

def findLiq(P_bar = None, T_C_init = None, comp = None): #, fO2_buffer = None, fO2_offset = None):
    '''
    Perform a single find liquidus calculation in MELTS. WARNING! Running this function directly from the command land/jupyter notebook will initiate the MELTS C library in the main python process. Once this has been initiated the MELTS C library cannot be re-loaded and failures during the calculation will likely cause a terminal error to occur.

    Parameters:
    ----------
    P_bar: float
        Specifies the pressure of calculation (bar).

    T_C_init: float
        Initial 'guess' temperature for findLiq calculations (degrees C).

    comp: list or dict
        Input oxide values required for the calculations.

    Returns:
    ---------
    T_Liq_C: np.ndarray
        Array of liquidus temperatures.

    '''

    if P_bar is None:
        raise Exception("Please specify a pressure for calculations")

    if T_C_init is None:
        T_C_init = 1300.0

    if comp is None:
        raise Exception("No composition specified")
    
    if type(comp) == dict:
        bulk = [comp['SiO2_Liq'], comp['Al2O3_Liq'], comp['CaO_Liq'], comp['MgO_Liq'], comp['FeOt_Liq'], comp['K2O_Liq'], comp['Na2O_Liq'], comp['TiO2_Liq'], comp['Fe3Fet_Liq']*(((159.59/2)/71.844)*comp['FeOt_Liq'] - comp['FeOt_Liq']), comp['Cr2O3_Liq'], comp['H2O_Liq']]
        T_Liq = MAGEMinCalc.findliq(bulk, P_bar/1000, T_C_init)
    else:
        T_Liq = MAGEMinCalc.findliq(comp, P_bar/1000, T_C_init)
        
    return T_Liq

    # else:
    #     if fO2_buffer is not None:
    #         new_bulk = pd.DataFrame(comp, index = [0])
    #         new_bulk['Sample_ID_Liq'] = 0
    #         if fO2_offset is not None:
    #             new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = T_C_init+273.15, P_kbar = P_bar/1000, fo2 = fO2_buffer, fo2_offset = fO2_offset, model = "Kress1991", renorm = False)
    #         else:
    #             new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = T_C_init+273.15, P_kbar = P_bar/1000, fo2 = fO2_buffer, model = "Kress1991", renorm = False)

    #         comp = new[list(comp.keys())].loc[0].to_dict()

    #     bulk = [comp['SiO2_Liq'], comp['Al2O3_Liq'], comp['CaO_Liq'], comp['MgO_Liq'], comp['FeOt_Liq'], comp['K2O_Liq'], comp['Na2O_Liq'], comp['TiO2_Liq'], comp['Fe3Fet_Liq']*(((159.59/2)/71.844)*comp['FeOt_Liq'] - comp['FeOt_Liq']), comp['Cr2O3_Liq'], comp['H2O_Liq']]

    # T_Liq = 0

    # T = T_C_init

    # if type(T) == int:
    #     T = float(T)

    # if type(P_bar) == int:
    #     P_bar = float(P_bar)

    # bulk_in = bulk.copy()

    # Liq = ["liq","fl"]

    # start = time.time()

    # Ret = MAGEMinCalc.PT_minimisation(P_bar/1000, T, bulk)
    # PhaseList = Ret['sys']['Phase']

    # i = set.intersection(set(Liq),set(PhaseList))

    # Step = np.array([3,1,0.1])
    # for k in range(len(Step)):
    #     if len(i) == len(PhaseList):
    #         while len(i) == len(PhaseList):
    #             bulk = bulk_in.copy()
    #             if time.time() - start > 60:
    #                 return T_Liq
    #             else:
    #                 T = T - Step[k]
    #                 if fO2_buffer is not None:
    #                     new_bulk = pd.DataFrame(comp, index = [0])
    #                     new_bulk['Sample_ID_Liq'] = 0
    #                     if fO2_offset is not None:
    #                         new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = T+273.15, P_kbar = P_bar/1000, fo2 = fO2_buffer, fo2_offset = fO2_offset, model = "Kress1991", renorm = False)
    #                     else:
    #                         new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = T+273.15, P_kbar = P_bar/1000, fo2 = fO2_buffer, model = "Kress1991", renorm = False)

    #                     comp = new[list(comp.keys())].loc[0].to_dict()

    #                     bulk = [comp['SiO2_Liq'], comp['Al2O3_Liq'], comp['CaO_Liq'], comp['MgO_Liq'], comp['FeOt_Liq'], comp['K2O_Liq'], comp['Na2O_Liq'], comp['TiO2_Liq'], comp['Fe3Fet_Liq']*(((159.59/2)/71.844)*comp['FeOt_Liq'] - comp['FeOt_Liq']), comp['Cr2O3_Liq'], comp['H2O_Liq']]

    #                 Ret = MAGEMinCalc.PT_minimisation(P_bar/1000, T, bulk)
    #                 PhaseList = Ret['sys']['Phase']
    #                 i = set.intersection(set(Liq),set(PhaseList))

    #     if len(i) < len(PhaseList):
    #         while len(i) < len(PhaseList):
    #             bulk = bulk_in.copy()
    #             if time.time() - start > 60:
    #                 return T_Liq
    #             else:
    #                 T = T + Step[k]
    #                 if fO2_buffer is not None:
    #                     new_bulk = pd.DataFrame(comp, index = [0])
    #                     new_bulk['Sample_ID_Liq'] = 0
    #                     if fO2_offset is not None:
    #                         new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = T+273.15, P_kbar = P_bar/1000, fo2 = fO2_buffer, fo2_offset = fO2_offset, model = "Kress1991", renorm = False)
    #                     else:
    #                         new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = T+273.15, P_kbar = P_bar/1000, fo2 = fO2_buffer, model = "Kress1991", renorm = False)

    #                     comp = new[list(comp.keys())].loc[0].to_dict()

    #                     bulk = [comp['SiO2_Liq'], comp['Al2O3_Liq'], comp['CaO_Liq'], comp['MgO_Liq'], comp['FeOt_Liq'], comp['K2O_Liq'], comp['Na2O_Liq'], comp['TiO2_Liq'], comp['Fe3Fet_Liq']*(((159.59/2)/71.844)*comp['FeOt_Liq'] - comp['FeOt_Liq']), comp['Cr2O3_Liq'], comp['H2O_Liq']]

    #                 Ret = MAGEMinCalc.PT_minimisation(P_bar/1000, T, bulk)
    #                 PhaseList = Ret['sys']['Phase']
    #                 i = set.intersection(set(Liq),set(PhaseList))

    # if "liq" in Ret['sys']['Phase']:
    #     if Ret['liq']['Frac'] < 0.8:
    #         return T_Liq
    #     else:
    #         T_Liq = T
    #         return T_Liq
    # else:
        # return T_Liq


def path(comp = None, Frac_solid = None, phases = None, Frac_fluid = None, T_C = None, T_min_C = None, T_path_C = None, T_start_C = None, T_end_C = None, dt_C = None, P_bar = None, P_path_bar = None, P_start_bar = None, P_end_bar = None, dp_bar = None, find_liquidus = None, fO2_buffer = None, fO2_offset = None):
    '''
    Perform a single  calculation in MELTS. WARNING! Running this function directly from the command land/jupyter notebook will initiate the MELTS C library in the main python process. Once this has been initiated the MELTS C library cannot be re-loaded and failures during the calculation will likely cause a terminal error to occur.

    Parameters:
    ----------
    comp: Dict
        Initial compositon for calculations.

    Frac_solid: True/False
        If True, solid phases will be removed from the system at the end of each crystallisation step. Default False.

    Frac_fluid: True/False
        If True, fluid phases will be removed from the system at the end of each crystallisation step. Default False.

    T_C: float
        Calculation temperature - typically used when calculations are performed at a fixed T (e.g.,isothermal degassing).

    T_min_C: float
        Temperature below the liquidus to calculate.

    T_path_C: np.ndarray
        If a specified temperature path is to be used, T_path_C will be used to determine the T at each step of the model. If 2D, this indicates that multiple calculations with different T_path_C arrays are to be performed.

    T_start_C: float
        Initial temperature used for path calculations.

    T_end_C: float
        Final temperature in crystallisation calculations.
melt
    dt_C: float
        Temperature increment during crystallisation calculations.

    P_bar: float
        Calculation pressure - typically used when calculations are performed at a fixed P (e.g.,isobaric crystallisation).

    P_path_bar: np.ndarray
        If a specified pressure path is to be used, P_path_bar will be used to determine the P at each step of the model.

    P_start_bar: float
        Initial pressure used for path calculations.

    P_end_bar: float
        Final pressure in crystallisation calculations.

    dp_bar: float
        Pressure increment during crystallisation calculations.

    find_liquidus: True/False
        If True, the calculations will start with a search for the melt liquidus temperature. Default is False.

    fO2_buffer: string
        If the oxygen fugacity of the system is to be buffered during crystallisation/decompression, then an offset to a known buffer must be specified. Here the user can define the known buffer as either "FMQ" or "NNO".

    fO2_offset: float
        Offset from the buffer spcified in fO2_buffer (log units).

    Returns:
    ----------
    Results: Dict
        Dict containing a series of pandas DataFrames that display the composition and thermodynamic properties of each phase.

    '''
    Results = {}

    if comp is None:
        raise Exception("No composition specified")

    if P_bar is not None and P_path_bar is None:
        P_path_bar = P_bar
    if T_C is not None and T_start_C is None:
        T_start_C = T_C

    if P_path_bar is None and P_start_bar is None:
        raise Exception("Initial P system must be defined")
    if T_path_C is None and T_start_C is None and find_liquidus is None:
        raise Exception("Starting temperature must be specified or the liquidus must be found")

    if comp is None:
        raise Exception("No composition specified")
    # else:
        # if fO2_buffer is not None:
        #     new_bulk = pd.DataFrame(comp, index = [0])
        #     new_bulk['Sample_ID_Liq'] = 0
        #     if fO2_offset is not None:
        #         new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = 1400+273.15, P_kbar = P_bar/1000, fo2 = fO2_buffer, fo2_offset = fO2_offset, model = "Kress1991", renorm = False)
        #     else:
        #         new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = 1400+273.15, P_kbar = P_bar/1000, fo2 = fO2_buffer, model = "Kress1991", renorm = False)

        #     comp = new[list(comp.keys())].loc[0].to_dict()

    bulk = [comp['SiO2_Liq'], comp['Al2O3_Liq'], comp['CaO_Liq'], comp['MgO_Liq'], comp['FeOt_Liq'], comp['K2O_Liq'], comp['Na2O_Liq'], comp['TiO2_Liq'], comp['Fe3Fet_Liq']*(((159.59/2)/71.844)*comp['FeOt_Liq'] - comp['FeOt_Liq']), comp['Cr2O3_Liq'], comp['H2O_Liq']]

    if find_liquidus is not None:
        if P_path_bar is not None:
            try:
                if type(P_path_bar) == np.ndarray:
                    T_Liq = findLiq(P_bar = P_path_bar[0], comp = bulk, T_C_init = 1400.0)
                else:
                    T_Liq = findLiq(P_bar = P_path_bar, comp = bulk, T_C_init = 1400.0)
            except:
                return Results
        elif P_start_bar is not None:
            try:
                T_Liq = findLiq(P_bar = P_start_bar, comp = bulk, T_C_init = 1400.0)
            except:
                return Results

        T_start_C = T_Liq
        if T_min_C is not None:
            T_end_C = T_Liq - T_min_C

        # if fO2_buffer is not None:
        #     new_bulk = pd.DataFrame(comp, index = [0])
        #     new_bulk['Sample_ID_Liq'] = 0
        #     if fO2_offset is not None:
        #         new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = T_Liq+273.15, P_kbar = P_bar/1000, fo2 = fO2_buffer, fo2_offset = fO2_offset, model = "Kress1991", renorm = False)
        #     else:
        #         new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = T_Liq+273.15, P_kbar = P_bar/1000, fo2 = fO2_buffer, model = "Kress1991", renorm = False)

        #     comp = new[list(comp.keys())].loc[0].to_dict()

        bulk = [comp['SiO2_Liq'], comp['Al2O3_Liq'], comp['CaO_Liq'], comp['MgO_Liq'], comp['FeOt_Liq'], comp['K2O_Liq'], comp['Na2O_Liq'], comp['TiO2_Liq'], comp['Fe3Fet_Liq']*(((159.59/2)/71.844)*comp['FeOt_Liq'] - comp['FeOt_Liq']), comp['Cr2O3_Liq'], comp['H2O_Liq']]

    if T_path_C is None:
        if T_end_C is None and dt_C is None:
            T = T_start_C
        elif T_end_C is not None and dt_C is not None:
            T = np.linspace(T_start_C, T_end_C, 1+round((T_start_C-T_end_C)/dt_C))
    elif T_path_C is not None:
        T = T_path_C

    if P_path_bar is None:
        if P_end_bar is None and dp_bar is None:
            P = P_start_bar
        elif P_end_bar is not None and dp_bar is not None:
            P = np.linspace(P_start_bar, P_end_bar, 1+round((P_start_bar-P_end_bar)/dp_bar))
    elif P_path_bar is not None:
        P = P_path_bar

    if type(T) == np.ndarray and P_end_bar is None and dp_bar is not None:
        P = np.linspace(P_start_bar, P_start_bar - dp_bar*(len(T)-1), len(T))
    elif type(P) == np.ndarray and T_end_C is None and dt_C is not None:
        T = np.linspace(T_start_C, T_start_C - dt_C*(len(P)-1), len(P))

    if type(T) == np.ndarray and type(P) == np.ndarray:
        if len(T) != len(P):
            raise Exception("Length of P and T vectors are not the same. Check input parameters")

    if find_liquidus is not None:
        if P_path_bar is not None or T_path_C is not None:
            if type(P) == np.ndarray and type(T) == np.ndarray:
                T_Liq_loc = np.abs(T - T_Liq).argmin()
                if T[T_Liq_loc]>T_Liq:
                    T = T[T_Liq_loc:]
                    P = P[T_Liq_loc:]
                else:
                    T = T[T_Liq_loc-1:]
                    P = P[T_Liq_loc-1:]

    if type(T) == np.ndarray and type(P) != np.ndarray:
        P = np.zeros(len(T)) + P
    if type(P) == np.ndarray and type(T) != np.ndarray:
        T = np.zeros(len(P)) + T

    if phases is None:
        if Frac_solid is not None:
            Results = MAGEMinCalc.path(bulk, T, P/1000, 1, 0)
        else:
            Results = MAGEMinCalc.path(bulk, T, P/1000, 0, 0)
    else:
        if Frac_solid is not None:
            Results = MAGEMinCalc.path(bulk, T, P/1000, 1, phases)
        else:
            Results = MAGEMinCalc.path(bulk, T, P/1000, 0, phases)

    return Results
    # if type(T) == np.ndarray:
    #     length = len(T)
    # else:
    #     length = len(P)

    # Results['Conditions'] = pd.DataFrame(data = np.zeros((length, 3)), columns = ['temperature', 'pressure', 's'])
    # Results['liq'] = pd.DataFrame(data = np.zeros((length, 11)), columns = ['SiO2', 'Al2O3', 'CaO', 'MgO', 'FeO', 'K2O', 'Na2O', 'TiO2', 'O', 'Cr2O3', 'H2O'])
    # Results['liq_prop'] = pd.DataFrame(data = np.zeros((length, 1)), columns = ['mass'])

    # if type(T) != np.ndarray:
    #     T_in = T
    # if type(P) != np.ndarray:
    #     P_in = P/1000

    # bulk_in = bulk.copy()

    # for i in range(length):
    #     if type(T) == np.ndarray:
    #         T_in = T[i]
    #     if type(P) == np.ndarray:
    #         P_in = P[i]/1000

    #     # print(T_in)

    #     if fO2_buffer is not None:
    #         new_comp = pd.DataFrame({'SiO2_Liq': bulk[0], 'Al2O3_Liq': bulk[1], 'CaO_Liq': bulk[2], 'MgO_Liq': bulk[3], 'FeOt_Liq': bulk[4], 'K2O_Liq': bulk[5], 'Na2O_Liq': bulk[6], 'TiO2_Liq': bulk[7], 'Cr2O3_Liq': bulk[9], 'H2O_Liq': bulk[10], 'Fe3Fet_Liq': 0.0}, index = [0])
    #         new_bulk = new_comp.copy()
    #         new_bulk['Sample_ID_Liq'] = 0
    #         if fO2_offset is not None:
    #             new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = T_in+273.15, P_kbar = P_in, fo2 = fO2_buffer, fo2_offset = fO2_offset, model = "Kress1991", renorm = False)
    #         else:
    #             new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = T_in+273.15, P_kbar = P_in, fo2 = fO2_buffer, model = "Kress1991", renorm = False)

    #         new_comp = new[list(new_comp.keys())].loc[0].to_dict()

    #         bulk = [new_comp['SiO2_Liq'], new_comp['Al2O3_Liq'], new_comp['CaO_Liq'], new_comp['MgO_Liq'], new_comp['FeOt_Liq'], new_comp['K2O_Liq'], new_comp['Na2O_Liq'], new_comp['TiO2_Liq'], new_comp['Fe3Fet_Liq']*(((159.59/2)/71.844)*new_comp['FeOt_Liq'] - new_comp['FeOt_Liq']), new_comp['Cr2O3_Liq'], new_comp['H2O_Liq']]

    #     try:
    #        Ret = MAGEMinCalc.PT_minimisation(P_in, T_in, bulk)
    #     except:
    #         return Results

    #     PhaseList = Ret['sys']['Phase']
    #     for phase in PhaseList:
    #         if phase not in list(Results.keys()):
    #             Results[phase] = pd.DataFrame(data = np.zeros((length, 11)), columns = ['SiO2', 'Al2O3', 'CaO', 'MgO', 'FeO', 'K2O', 'Na2O', 'TiO2', 'O', 'Cr2O3', 'H2O'])
    #             Results[phase + '_prop'] = pd.DataFrame(data = np.zeros((length, 1)), columns = ['mass'])

    #     for R in Results['Conditions']:
    #         if R == 'temperature':
    #             Results['Conditions'][R].loc[i] = T_in
    #         elif R == 'pressure':
    #             Results['Conditions'][R].loc[i] = P_in*1000
    #         else:
    #             Results['Conditions'][R].loc[i] = Ret['sys']['Entropy']

    #     for p in PhaseList:
    #         if "Comp" in list(Ret[p].keys()):
    #             # A = dict(zip(['SiO2', 'Al2O3', 'CaO', 'MgO', 'FeO', 'Na2O', 'K2O', 'TiO2', 'O', 'Cr2O3', 'H2O'], Ret[p]['Comp']))
    #             # print(A)
    #             for el in Results[p]:
    #                 Results[p][el].loc[i] = Ret[p]["Comp"][el]

    #         Results[p + '_prop']['mass'].loc[i] = Ret[p]['Frac']

    #     if Frac_solid is not None and Frac_fluid is not None and i > 0:
    #         bulk = np.array([Results['liq']['SiO2'].loc[i], Results['liq']['Al2O3'].loc[i], Results['liq']['CaO'].loc[i], Results['liq']['MgO'].loc[i], Results['liq']['FeO'].loc[i], Results['liq']['K2O'].loc[i], Results['liq']['Na2O'].loc[i], Results['liq']['TiO2'].loc[i], Results['liq']['O'].loc[i], Results['liq']['Cr2O3'].loc[i], Results['liq']['H2O'].loc[i]])
    #         # bulk = Ret['liq']['Comp']
    #         comp = bulk.copy()
    #         bulk = 100*comp/sum(comp)
    #         # print(bulk)
    #     else:
    #         bulk = bulk_in

    # return Results


def AdiabaticDecompressionMelting(comp = None, T_p_C = None, P_start_kbar = None, P_end_kbar = None, dp_kbar = None, Frac = None):
    lz = m.lithologies.matthews.klb1()
    mantle = m.mantle([lz], [1], ['Lz'])
    T_start_C = mantle.adiabat(P_start_kbar/10, T_p_C)
    bulk = [comp['SiO2_Liq'], comp['Al2O3_Liq'], comp['CaO_Liq'], comp['MgO_Liq'], comp['FeOt_Liq'], comp['K2O_Liq'], comp['Na2O_Liq'], comp['TiO2_Liq'], comp['Fe3Fet_Liq']*(((159.59/2)/71.844)*comp['FeOt_Liq'] - comp['FeOt_Liq']), comp['Cr2O3_Liq'], comp['H2O_Liq']]
    
    if Frac is None:
        Frac = 0

    Results = MAGEMinCalc.AdiabaticDecompressionMelting(bulk, T_start_C, P_start_kbar, P_end_kbar, dp_kbar, Frac)

    return Results

    # P_bar = np.linspace(P_start_bar, P_end_bar, round((P_start_bar - P_end_bar)/dp_bar))

    # bulk_in = bulk.copy()

    # Results = {}
    # Results['Conditions'] = pd.DataFrame(data = np.zeros((len(P_bar), 3)), columns = ['temperature', 'pressure', 's'])
    # for i in range(len(P_bar)):
    #     P = P_bar[i]
    #     if i == 0:
    #         try:
    #             Ret = MAGEMinCalc.PT_minimisation(P/1000, float(T_start_C), bulk)
    #             T = T_start_C
    #         except:
    #             return Results
    #     else:
    #         try:
    #             Ret = MAGEMinCalc.PS_minimisation(s, P/1000, float(T), bulk)
    #             T = Ret['sys']['Temperature']
    #         except:
    #             return Results

    #     PhaseList = Ret['sys']['Phase']
    #     for phase in PhaseList:
    #         if phase not in list(Results.keys()):
    #             Results[phase] = pd.DataFrame(data = np.zeros((len(P_bar), 11)), columns = ['SiO2', 'Al2O3', 'CaO', 'MgO', 'FeO', 'K2O', 'Na2O', 'TiO2', 'O', 'Cr2O3', 'H2O'])
    #             Results[phase + '_prop'] = pd.DataFrame(data = np.zeros((len(P_bar), 1)), columns = ['mass'])

    #     for R in Results['Conditions']:
    #         if R == 'temperature':
    #             Results['Conditions'][R].loc[i] = T
    #         elif R == 'pressure':
    #             Results['Conditions'][R].loc[i] = P_bar[i]
    #         else:
    #             Results['Conditions'][R].loc[i] = Ret['sys']['Entropy']

    #     for p in PhaseList:
    #         if "Comp" in list(Ret[p].keys()):
    #             for el in Results[p]:
    #                 Results[p][el].loc[i] = Ret[p]["Comp"][el]

    #         Results[p + '_prop']['mass'].loc[i] = Ret[p]['Frac']

    #     s = Ret['sys']['Entropy']

    #     bulk = bulk_in.copy()
# def phaseSat(comp = None, phases = None, T_initial_C = None, T_step_C = None, dt_C = None, P_bar = None, H2O_Liq = None, fO2_buffer = None, fO2_offset = None):
#     '''
#     Perform a single crystallisation calculation in MELTS. WARNING! Running this function directly from the command land/jupyter notebook will initiate the MELTS C library in the main python process. Once this has been initiated the MELTS C library cannot be re-loaded and failures during the calculation will likely cause a terminal error to occur.
#
#     Parameters:
#     ----------
#     comp: list or dict
#         Input oxide values required for the calculations.
#
#     phases: list
#         phases of interest
#
#     T_initial_C: float
#         Initial temperature used for liquidus calculations.
#
#     T_step_C: float
#         Temperature step at each point of the model.
#
#     dt_C: float
#         Total temperature change allowed in the model.
#
#     P_bar: float
#         Pressure of the calculation.
#
#     Returns:
#     ----------
#     Results: Dict
#         Dict containing a float for each saturation temperature found and the T_Liq and melt H2O values.
#
#     '''
#     Results = {'a_sat': np.nan, 'b_sat': np.nan, 'c_sat': np.nan, 'T_Liq': np.nan}
#     if len(phases) == 2:
#         del Results['c_sat']
#
#     bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], comp['Fe3Fet_Liq']*((159.59/2)/71.844)*comp['FeOt_Liq'], 0.0, (1- comp['Fe3Fet_Liq'])*comp['FeOt_Liq'], comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]
#     bulk = list(100*np.array(bulk)/np.sum(bulk))
#
#     try:
#         Results['T_Liq'] = findLiq(P_bar = P_bar, comp = comp, T_C_init = T_initial_C, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
#     except:
#         return Results
#
#     if type(H2O_Liq) == np.ndarray:
#         if Results['H2O_melt'] < 0.99*bulk[14]:
#             return Results
#
#
#     T = Results['T_Liq']
#     T_final = T - dt_C
#     while T >= T_final:
#         melts = melts.addNodeAfter()
#         melts.engine.setBulkComposition(bulk)
#         melts.engine.pressure = P_bar
#         melts.engine.temperature = T
#
#         try:
#             melts.engine.calcEquilibriumState(1,0)
#         except:
#             return Results
#
#         PhaseList = melts.engine.solidNames
#         print(PhaseList)
#         try:
#             if 'tridymite1' in PhaseList:
#                 PhaseList = ['quartz1'] + PhaseList
#             if 'clinopyroxene2' in PhaseList:
#                 PhaseList = ['orthopyroxene1'] + PhaseList
#
#             if phases[0] in PhaseList and np.isnan(Results['a_sat']):# == 0:
#                 Results['a_sat'] = melts.engine.temperature
#
#             if phases[1] in PhaseList and np.isnan(Results['b_sat']):# == 0:
#                 Results['b_sat'] = melts.engine.temperature
#
#             if len(phases) == 3:
#                 if phases[2] in PhaseList and np.isnan(Results['c_sat']):# == 0:
#                     Results['c_sat'] = melts.engine.temperature
#
#                 if ~np.isnan(Results['a_sat']) and ~np.isnan(Results['b_sat']) and ~np.isnan(Results['c_sat']):# > 0:
#                     break
#
#             if len(phases) == 2:
#                 if ~np.isnan(Results['a_sat']) and ~np.isnan(Results['b_sat']):# > 0:
#                     break
#
#             T = T - T_step_C
#         except:
#             T = T - T_step_C
#
#
#     return Results