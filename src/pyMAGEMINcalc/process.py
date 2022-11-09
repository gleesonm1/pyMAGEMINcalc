import numpy as np
import pandas as pd
import Thermobar as pt
import julia
import time
from julia import MAGEMinCalc

def findLiq(P_bar = None, T_C_init = None, comp = None, fO2_buffer = None, fO2_offset = None):
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
    else:
        if fO2_buffer is not None:
            new_bulk = pd.DataFrame(comp, index = [0])
            new_bulk['Sample_ID_Liq'] = 0
            if fO2_offset is not None:
                new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = T_C_init+273.15, P_kbar = P_bar/1000, fo2 = fO2_buffer, fo2_offset = fO2_offset, model = "Kress1991", renorm = False)
            else:
                new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = T_C_init+273.15, P_kbar = P_bar/1000, fo2 = fO2_buffer, model = "Kress1991", renorm = False)

            comp = new[list(comp.keys())].loc[0].to_dict()

        bulk = [comp['SiO2_Liq'], comp['Al2O3_Liq'], comp['CaO_Liq'], comp['MgO_Liq'], comp['FeOt_Liq'], comp['K2O_Liq'], comp['Na2O_Liq'], comp['TiO2_Liq'], comp['Fe3Fet_Liq']*(((159.59/2)/71.844)*comp['FeOt_Liq'] - comp['FeOt_Liq']), comp['Cr2O3_Liq'], comp['H2O_Liq']]

    T_Liq = 0

    T = T_C_init

    bulk_in = bulk.copy()

    Liq = ["liq","fl"]

    start = time.time()

    PhaseList = MAGEMinCalc.satPhase(P_bar/1000, T, bulk)

    i = set.intersection(set(Liq),set(PhaseList))

    Step = np.array([3,1,0.1])
    for k in range(len(Step)):
        if len(i) == len(PhaseList):
            while len(i) == len(PhaseList):
                bulk = bulk_in.copy()
                if time.time() - start > 60:
                    return T_Liq
                else:
                    T = T - Step[k]
                    if fO2_buffer is not None:
                        new_bulk = pd.DataFrame(comp, index = [0])
                        new_bulk['Sample_ID_Liq'] = 0
                        if fO2_offset is not None:
                            new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = T+273.15, P_kbar = P_bar/1000, fo2 = fO2_buffer, fo2_offset = fO2_offset, model = "Kress1991", renorm = False)
                        else:
                            new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = T+273.15, P_kbar = P_bar/1000, fo2 = fO2_buffer, model = "Kress1991", renorm = False)

                        comp = new[list(comp.keys())].loc[0].to_dict()

                        bulk = [comp['SiO2_Liq'], comp['Al2O3_Liq'], comp['CaO_Liq'], comp['MgO_Liq'], comp['FeOt_Liq'], comp['K2O_Liq'], comp['Na2O_Liq'], comp['TiO2_Liq'], comp['Fe3Fet_Liq']*(((159.59/2)/71.844)*comp['FeOt_Liq'] - comp['FeOt_Liq']), comp['Cr2O3_Liq'], comp['H2O_Liq']]

                    Ret = MAGEMinCalc.satPhase(P_bar/1000, T, bulk)
                    PhaseList = Ret['Phase']
                    i = set.intersection(set(Liq),set(PhaseList))

        if len(i) < len(PhaseList):
            while len(i) < len(PhaseList):
                bulk = bulk_in.copy()
                if time.time() - start > 60:
                    return T_Liq
                else:
                    T = T + Step[k]
                    if fO2_buffer is not None:
                        new_bulk = pd.DataFrame(comp, index = [0])
                        new_bulk['Sample_ID_Liq'] = 0
                        if fO2_offset is not None:
                            new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = T+273.15, P_kbar = P_bar/1000, fo2 = fO2_buffer, fo2_offset = fO2_offset, model = "Kress1991", renorm = False)
                        else:
                            new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = T+273.15, P_kbar = P_bar/1000, fo2 = fO2_buffer, model = "Kress1991", renorm = False)

                        comp = new[list(comp.keys())].loc[0].to_dict()

                        bulk = [comp['SiO2_Liq'], comp['Al2O3_Liq'], comp['CaO_Liq'], comp['MgO_Liq'], comp['FeOt_Liq'], comp['K2O_Liq'], comp['Na2O_Liq'], comp['TiO2_Liq'], comp['Fe3Fet_Liq']*(((159.59/2)/71.844)*comp['FeOt_Liq'] - comp['FeOt_Liq']), comp['Cr2O3_Liq'], comp['H2O_Liq']]

                    Ret = MAGEMinCalc.satPhase(P_bar/1000, T, bulk)
                    PhaseList = Ret['Phase']
                    i = set.intersection(set(Liq),set(PhaseList))

    if "liq" in Ret['Phase']:
        if Ret['Liq_Frac'] < 0.9:
            return T_Liq
        else:
            T_Liq = T
            return T_Liq
    else:
        return T_Liq


def path(comp = None, Frac_solid = None, Frac_fluid = None, T_C = None, T_path_C = None, T_start_C = None, T_end_C = None, dt_C = None, P_bar = None, P_path_bar = None, P_start_bar = None, P_end_bar = None, dp_bar = None, find_liquidus = None, fO2_buffer = None, fO2_offset = None):
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
    else:
        if fO2_buffer is not None:
            new_bulk = pd.DataFrame(comp, index = [0])
            new_bulk['Sample_ID_Liq'] = 0
            if fO2_offset is not None:
                new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = 1400+273.15, P_kbar = P_bar/1000, fo2 = fO2_buffer, fo2_offset = fO2_offset, model = "Kress1991", renorm = False)
            else:
                new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = 1400+273.15, P_kbar = P_bar/1000, fo2 = fO2_buffer, model = "Kress1991", renorm = False)

            comp = new[list(comp.keys())].loc[0].to_dict()

        bulk = [comp['SiO2_Liq'], comp['Al2O3_Liq'], comp['CaO_Liq'], comp['MgO_Liq'], comp['FeOt_Liq'], comp['K2O_Liq'], comp['Na2O_Liq'], comp['TiO2_Liq'], comp['Fe3Fet_Liq']*(((159.59/2)/71.844)*comp['FeOt_Liq'] - comp['FeOt_Liq']), comp['Cr2O3_Liq'], comp['H2O_Liq']]

    if find_liquidus is not None:
        if P_path_bar is not None:
            try:
                if type(P_path_bar) == np.ndarray:
                    T_Liq = findLiq(P_bar = P_path_bar[0], comp = comp, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, T_C_init = 1400.0)
                else:
                    T_Liq = findLiq(P_bar = P_path_bar, comp = comp, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, T_C_init = 1400.0)
            except:
                return Results
        elif P_start_bar is not None:
            try:
                T_Liq = findLiq(P_bar = P_start_bar, comp = comp, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, T_C_init = 1400.0)
            except:
                return Results

        T_start_C = T_Liq

        if fO2_buffer is not None:
            new_bulk = pd.DataFrame(comp, index = [0])
            new_bulk['Sample_ID_Liq'] = 0
            if fO2_offset is not None:
                new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = T_Liq+273.15, P_kbar = P_bar/1000, fo2 = fO2_buffer, fo2_offset = fO2_offset, model = "Kress1991", renorm = False)
            else:
                new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = T_Liq+273.15, P_kbar = P_bar/1000, fo2 = fO2_buffer, model = "Kress1991", renorm = False)

            comp = new[list(comp.keys())].loc[0].to_dict()

        bulk = [comp['SiO2_Liq'], comp['Al2O3_Liq'], comp['CaO_Liq'], comp['MgO_Liq'], comp['FeOt_Liq'], comp['K2O_Liq'], comp['Na2O_Liq'], comp['TiO2_Liq'], comp['Fe3Fet_Liq']*(((159.59/2)/71.844)*comp['FeOt_Liq'] - comp['FeOt_Liq']), comp['Cr2O3_Liq'], comp['H2O_Liq']]

    if T_path_C is None:
        if T_end_C is None and dt is None:
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

    if type(T) == np.ndarray:
        length = len(T)
    else:
        length = len(P)

    Results['Conditions'] = pd.DataFrame(data = np.zeros((length, 2)), columns = ['temperature', 'pressure'])
    Results['liq'] = pd.DataFrame(data = np.zeros((length, 11)), columns = ['SiO2_Liq', 'Al2O3_Liq', 'CaO_Liq', 'MgO_Liq', 'FeOt_Liq', 'K2O_Liq', 'Na2O_Liq', 'TiO2_Liq', 'O_Liq', 'Cr2O3_Liq', 'H2O_Liq'])

    if type(T) != np.ndarray:
        T_in = T
    if type(P) != np.ndarray:
        P_in = P/1000

    for i in range(length):
        if type(T) == np.ndarray:
            T_in = T[i]
        if type(P) == np.ndarray:
            P_in = P[i]/1000

        if fO2_buffer is not None:
            new_comp = pd.DataFrame({'SiO2_Liq': bulk[0], 'Al2O3_Liq': bulk[1], 'CaO_Liq': bulk[2], 'MgO_Liq': bulk[3], 'FeOt_Liq': bulk[4], 'K2O_Liq': bulk[5], 'Na2O_Liq': bulk[6], 'TiO2_Liq': bulk[7], 'Cr2O3_Liq': bulk[9], 'H2O_Liq': bulk[10], 'Fe3Fet_Liq': 0.0}, index = [0])
            new_bulk = new_comp.copy()
            new_bulk['Sample_ID_Liq'] = 0
            if fO2_offset is not None:
                new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = T_in+273.15, P_kbar = P_in, fo2 = fO2_buffer, fo2_offset = fO2_offset, model = "Kress1991", renorm = False)
            else:
                new = pt.convert_fo2_to_fe_partition(liq_comps = new_bulk, T_K = T_in+273.15, P_kbar = P_in, fo2 = fO2_buffer, model = "Kress1991", renorm = False)

            new_comp = new[list(new_comp.keys())].loc[0].to_dict()

            bulk = [new_comp['SiO2_Liq'], new_comp['Al2O3_Liq'], new_comp['CaO_Liq'], new_comp['MgO_Liq'], new_comp['FeOt_Liq'], new_comp['K2O_Liq'], new_comp['Na2O_Liq'], new_comp['TiO2_Liq'], new_comp['Fe3Fet_Liq']*(((159.59/2)/71.844)*new_comp['FeOt_Liq'] - new_comp['FeOt_Liq']), new_comp['Cr2O3_Liq'], new_comp['H2O_Liq']]

        try:
           Ret = MAGEMinCalc.satPhase(P_in, T_in, bulk)
        except:
            return Results

        for R in Results['Conditions']:
            if R == 'temperature':
                Results['Conditions'][R].loc[i] = T_in
            elif R == 'pressure':
                Results['Conditions'][R].loc[i] = P_in*1000

        A = dict(zip(Ret['Oxides'], Ret['Liq_Comp']))
        for el in Results['liq']:
            if el != 'FeOt_Liq':
                Results['liq'][el].loc[i] = A[el[:-4]]
            else:
                Results['liq'][el].loc[i] = A['FeO']

        if Frac_solid is not None and Frac_fluid is not None:
            bulk = np.array(Ret['Liq_Comp'])
        else:
            bulk = bulk_in

    return Results