import numpy as np
from scipy.interpolate import interp1d
import json
from scipy.integrate import simpson, quad
from scipy.special import roots_legendre
from scipy.interpolate import CubicSpline
from .solar_sp import solar_sp

# Constants
PI = np.pi
D2R = PI / 180.0
C = 2.9979e14  # Speed of light (micron/sec)
H = 6.6262e-34  # Planck constant
K = 1.3807e-23  # Boltzmann constant
STEPHAN_BOLTZMANN = 5.6697e-12  # Stephan-Boltzmann constant
COEFF = 15.0 / PI**5 * STEPHAN_BOLTZMANN
RMAX = 5.2  # Jupiter's orbit radius
NUMBER_OF_STEPS = 50


class RZODI:
    """
    A class to represent the RZODI structure.
    """
    _defined = False

    def __init__(self):
        if not RZODI._defined:
            self.PIXEL_NO = np.int64(0)
            self.DAY1990 = 0.0
            self.WAVE_LEN = 0.0
            self.LONGITUDE = 0.0
            self.LATITUDE = 0.0
            self.PHOTOMETRY = 0.0
            RZODI._defined = True


def ZODI_str(Nstruct=1):
    """
    Function to create an array of RZODI structures.
    """
    if not isinstance(Nstruct, int):
        Nstruct = 1
    return [RZODI() for _ in range(Nstruct)]


def mk_zdata(lambda_val, day, lon, lat):
    """
    Generates a zmodel input dataset.
    """
    nwave = 1 if isinstance(lambda_val, (int, float)) else len(lambda_val)
    nday = 1 if isinstance(day, (int, float)) else len(day)
    nlon = 1 if isinstance(lon, (int, float)) else len(lon)
    nlat = 1 if isinstance(lat, (int, float)) else len(lat)

    npts = max(nwave, nday, nlon, nlat)

    if not all(x in [1, npts] for x in [nwave, nday, nlon, nlat]):
        raise ValueError("Inputs must be scalars or arrays of equal dimension.")

    data = ZODI_str(npts)

    if isinstance(lambda_val, (int, float)):
        for i in range(npts):
            data[i].WAVE_LEN = float(lambda_val)
    else:
        for i in range(npts):
            data[i].WAVE_LEN = float(lambda_val[i])

    if isinstance(day, (int, float)):
        for i in range(npts):
            data[i].DAY1990 = float(day)
    else:
        for i in range(npts):
            data[i].DAY1990 = float(day[i])

    if isinstance(lon, (int, float)):
        for i in range(npts):
            data[i].LONGITUDE = float(lon)
    else:
        for i in range(npts):
            data[i].LONGITUDE = float(lon[i])

    if isinstance(lat, (int, float)):
        for i in range(npts):
            data[i].LATITUDE = float(lat)
    else:
        for i in range(npts):
            data[i].LATITUDE = float(lat[i])

    return data


def read_zpars(filename='aend_og_wlabels.json'):
    """
    Reads the ZODI model parameters from a JSON file.
    """
    try:
        with open(filename, 'r') as f:
            data = json.load(f)
        lambda_z = np.array(data['lambda_z'])
        zpar = np.array(data['zpar'])
        return lambda_z, zpar
    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        return None, None
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON from '{filename}': {e}")
        return None, None
    except KeyError as e:
        print(f"Error: Key '{e}' not found in JSON file.")
        return None, None
    except Exception as e:
        print(f"Error reading data from '{filename}': {e}")
        return None, None

def colcorr(det, t, return_derivative=False):
    # Initialize static variables
    if not hasattr(colcorr, 'init'):
        colcorr.init = True

        infile = 'colcorr.tab'
        print('Reading color correction file:', infile)

        # Read data from colcorr.tab
        data = np.loadtxt(infile)
        intemp = data[:, 0]
        incc = data[:, 1:].T  # shape (10, len(intemp))

        # Define interpolation range
        tmax = 1000.
        tmin = 10.
        dt = 5.
        temp = np.arange(tmin, tmax + dt, dt)

        nt = len(temp)

        # Build tables of CC and dCdT for each detector
        colcorr.table = np.zeros((nt, 2, 10))

        for i in range(10):
            spline = CubicSpline(intemp, incc[i])
            cc = spline(temp)
            cc2 = spline(temp + dt)

            colcorr.table[:, 0, i] = cc
            colcorr.table[:, 1, i] = (cc2 - cc) / dt

        colcorr.temp = temp
        colcorr.tmin = tmin
        colcorr.dt = dt
        colcorr.nt = nt

    # Compute CC and dCdT at each temperature in t
    it = np.clip(((t - colcorr.tmin) / colcorr.dt).astype(int), 0, colcorr.nt - 1)
    delta_t = t - colcorr.temp[it]

    dcdt = colcorr.table[it, 1, det]
    cc = colcorr.table[it, 0, det] + dcdt * delta_t

    if return_derivative:
        return cc, dcdt
    else:
        return cc

def zPlanck(T, lam_um, return_derivative=False):
    """
    Planck function in units of W/cm^2/sr, matching the IDL zplanck.
    T       : temperature in Kelvin
    lam_um  : wavelength in microns (μm)
    return_derivative : if True, also returns derivative dB/dT
    """
    PI = np.pi
    c = 2.9979e14  # micron/sec
    h = 6.6262e-34
    k = 1.3807e-23
    stephan_boltzmann = 5.6697e-12
    coeff = 15.0 / PI**5 * stephan_boltzmann

    T = np.asarray(T, dtype=np.float64)
    lam = np.asarray(lam_um, dtype=np.float64)
    chi = h * c / (lam * k * T)
    val = np.zeros_like(chi)
    dBdT = np.zeros_like(chi)

    rj = chi < 0.001
    if np.any(rj):
        val[rj] = coeff * T[rj]**4 * chi[rj]**3
        dBdT[rj] = val[rj] / T[rj]

    exact = (chi >= 0.001) & (chi <= 50)
    if np.any(exact):
        chi_ex = chi[exact]
        val[exact] = coeff * T[exact]**4 * chi_ex**4 / (np.exp(chi_ex) - 1.0)
        dBdT[exact] = val[exact] * chi_ex / T[exact] / (1.0 - np.exp(-chi_ex))

    wien = chi > 50.0
    if np.any(wien):
        chi_wn = chi[wien]
        val[wien] = coeff * T[wien]**4 * chi_wn**4 * np.exp(-chi_wn)
        dBdT[wien] = val[wien] * chi_wn / T[wien] / (1.0 - np.exp(-chi_wn))

    if np.isscalar(chi):
        val = val.item()
        dBdT = dBdT.item()

    if return_derivative:
        return val, dBdT
    else:
        return val

def gaussint(a, b, numpts=24):
    """
    Returns Gauss-Laguerre abscissae and weights for finite interval [a,b].

    CRITICAL FIX: The IDL implementation uses Gauss-Laguerre quadrature with
    specific pre-computed nodes and weights, not Gauss-Legendre. This function
    now exactly matches the IDL gaussint implementation to fix the 22% systematic error.
    """
    import numpy as np

    if numpts == 24:
        # IDL grid24 values (exact from IDL source code)
        grid_raw = np.array([
            0.0567047755, 0.2990108986, 0.7359095554, 1.3691831160,
            2.2013260537, 3.2356758036, 4.4764966151, 5.9290837627,
            7.5998993100, 9.4967492209, 11.629014912, 14.007957977,
            16.647125597, 19.56289801, 22.775241987, 26.308772391,
            30.194291163, 34.471097572, 39.190608804, 44.422349336,
            50.264574993, 56.864967174, 64.466670616, 73.534234792
        ])

        # IDL wts24 values (exact from IDL source code)
        wts_raw = np.array([
            0.1455497377, 0.3393497722, 0.5347365921, 0.7322248725,
            0.9326159014, 1.1367925904, 1.3457293379, 1.560519046,
            1.7824092263, 2.0128491498, 2.2535525026, 2.5065825127,
            2.7744704429, 3.0603848697, 3.3683805666, 3.7037765832,
            4.0737527888, 4.4883345170, 4.9621093140, 5.5174318658,
            6.19095, 7.0486426246, 8.2265276593, 10.0758794961
        ])

        # IDL normalization factor
        scale_factor = 73.534234792

        # IDL scaling procedure:
        # 1. Normalize by dividing by scale_factor
        # 2. Scale to interval [a,b]: grid = grid_norm * (b-a) + a, wts = wts_norm * (b-a)
        grid_norm = grid_raw / scale_factor
        wts_norm = wts_raw / scale_factor

        grid = grid_norm * (b - a) + a
        wts = wts_norm * (b - a)

        return grid, wts

    elif numpts == 50:
        # IDL grid50 values (exact from IDL source code)
        grid_raw = np.array([
            0.0286305183, 0.1508829357, 0.3709487815, 0.6890906999,
            1.1056250235, 1.6209617511, 2.2356103759, 2.9501833666,
            3.7653997744, 4.6820893876, 5.7011975748, 6.8237909098,
            8.0510636694, 9.3843453083, 10.825109031, 12.374981608,
            14.035754599, 15.809397197, 17.698070933, 19.704146535,
            21.830223306, 24.079151444, 26.454057841, 28.958376011,
            31.595880956, 34.370729963, 37.287510610, 40.351297573,
            43.567720270, 46.943043991, 50.484267963, 54.199244880,
            58.096828017, 62.187054175, 66.481373878, 70.992944826,
            75.737011547, 80.731404802, 85.997211136, 91.559690412,
            97.449565614, 103.70489123, 110.37385880, 117.51919820,
            125.22547013, 133.61202792, 142.85832548, 153.26037197,
            165.3856433, 180.69834370
        ])

        # IDL wts50 values (exact from IDL source code)
        wts_raw = np.array([
            0.0734786263, 0.1711133062, 0.2690583985, 0.3672778840,
            0.4658590232, 0.5648992928, 0.6644999757, 0.7647657788,
            0.8658052545, 0.9677314398, 1.0706625889, 1.1747230029,
            1.2800439495, 1.3867646993, 1.4950336890, 1.6050098453,
            1.7168640913, 1.8307810796, 1.9469611911, 2.0656228595,
            2.1870052872, 2.3113716402, 2.4390128272, 2.5702520038,
            2.7054499706, 2.8450116969, 2.9893942558, 3.1391165624,
            3.2947714220, 3.4570405806, 3.6267137126, 3.8047126447,
            3.9921226303, 4.1902332693, 4.4005928365, 4.6250816069,
            4.8660126478, 5.1262732767, 5.4095283377, 5.7205203736,
            6.0655271405, 6.4530854743, 6.8951889633, 7.4093807809,
            8.0226692310, 8.7795282602, 9.7603031266, 11.131390723,
            13.324267692, 18.114146002
        ])

        # IDL normalization factor for 50-point grid
        scale_factor = 180.69834370

        # IDL scaling procedure:
        # 1. Normalize by dividing by scale_factor
        # 2. Scale to interval [a,b]: grid = grid_norm * (b-a) + a, wts = wts_norm * (b-a)
        grid_norm = grid_raw / scale_factor
        wts_norm = wts_raw / scale_factor

        grid = grid_norm * (b - a) + a
        wts = wts_norm * (b - a)

        return grid, wts
    else:
        # Fallback for other sizes
        from scipy.special import roots_legendre
        x, w = roots_legendre(numpts)
        grid = 0.5 * (b - a) * x + 0.5 * (b + a)
        wts = 0.5 * (b - a) * w
        return grid, wts


def hong_phase_func(x, a):
    """
    Hong phase function for SKYSURF model.
    Implementation based on O'Brien+2025 and IDL phase_func_fns.pro

    Parameters:
    -----------
    x : array-like
        Scattering angles
    a : array-like
        Hong parameters [g1, g2, g3, w1, w2, w3]

    Returns:
    --------
    phi : np.ndarray
        Phase function values
    """
    x = np.asarray(x, dtype=np.float64)
    a = np.asarray(a, dtype=np.float64)

    nf = len(a) // 2  # Number of components
    phi = np.zeros_like(x, dtype=np.float64)

    for i in range(nf):
        g = a[i]
        w = a[i + nf]
        # Henyey-Greenstein phase function component
        phi += w * (1 - g**2) / (1 + g**2 - 2*g*np.cos(x))**1.5

    return phi

def phasefunc(x, a):
    """
    Trial Phase Function for ZKERNEL (Kelsall model)
    Returns:
        phi (np.ndarray): phase function
        pder (np.ndarray): partial derivatives w.r.t. each parameter
    """
    x = np.asarray(x, dtype=np.float64)

    C0 = a[0]
    C1 = a[1]
    ka = a[2]

    pisq = PI**2
    picu = PI**3
    q0 = 2.0
    q1 = PI
    qk = (np.exp(ka * PI) + 1.0) / (ka**2 + 1.0)
    v = (2.0 * PI) * (C0 * q0 + C1 * q1 + qk)

    u = C0 + C1 * x + np.exp(ka * x)
    phi = u / v

    # Partial derivatives
    npts = len(x)
    pder = np.zeros((npts, 3), dtype=np.float64)

    # ∂phi/∂C0
    pder[:, 0] = (1.0 - phi * 4.0 * PI) / v

    # ∂phi/∂C1
    pder[:, 1] = (x - 2.0 * PI * phi * q1) / v

    # ∂phi/∂ka
    term1 = x * np.exp(ka * x)
    term2_numer = (ka**2 + 1.0) * PI * np.exp(ka * PI) - 2.0 * ka * (np.exp(ka * PI) + 1.0)
    term2_denom = (ka**2 + 1.0)**2
    term2 = 2.0 * PI * phi * term2_numer / term2_denom
    pder[:, 2] = (term1 - term2) / v

    return phi, pder


def get_parnames():
    npars = 256
    parname = ['Unused'] * npars

    parname[0] = 'FuncIndx'
    parname[1:10] = ['PF1_C0', 'PF1_C1', 'PF1_Ka', 'PF2_C0', 'PF2_C1', 'PF2_Ka', 'PF3_C0', 'PF3_C1', 'PF3_Ka']
    parname[10:14] = ['To1', 'To2', 'K', 'Delta']
    parname[14:33] = ['C_No', 'C_Alpha', 'C_Beta', 'C_Gamma', 'C_Mu', 'C_Mu2', 'C_Mu3', 'C_Mu4', 'C_Mu5', 'C_Mu6', 'C_Mu7', 'C_Mu8', 'C_Mu9', 'C_Mu10', 'C_Omega', 'C_Incl', 'C_Xo', 'C_Yo', 'C_Zo']
    parname[33:45] = ['C_A1', 'C_A2', 'C_A3', 'C_A4', 'C_E3', 'C_E4', 'C_E5', 'C_E6', 'C_E7', 'C_E8', 'C_E9', 'C_E10']
    parname[45:71] = ['B1_No', 'B1_Dz', 'B1_Dr', 'B1_Ro', 'B1_Vi', 'B1_Vr', 'B1_Pi', 'B1_Pr', 'B1_P2r', 'B1_Omega', 'B1_Incl', 'B1_Xo', 'B1_Yo', 'B1_Zo', 'B1_A1', 'B1_A2', 'B1_A3', 'B1_A4', 'B1_E3', 'B1_E4', 'B1_E5', 'B1_E6', 'B1_E7', 'B1_E8', 'B1_E9', 'B1_E10']
    parname[71:97] = ['B2_No', 'B2_Dz', 'B2_Dr', 'B2_Ro', 'B2_Vi', 'B2_Vr', 'B2_Pi', 'B2_Pr', 'B2_P2r', 'B2_Omega', 'B2_Incl', 'B2_Xo', 'B2_Yo', 'B2_Zo', 'B2_A1', 'B2_A2', 'B2_A3', 'B2_A4', 'B2_E3', 'B2_E4', 'B2_E5', 'B2_E6', 'B2_E7', 'B2_E8', 'B2_E9', 'B2_E10']
    parname[97:123] = ['B3_No', 'B3_Dz', 'B3_Dr', 'B3_Ro', 'B3_Vi', 'B3_Vr', 'B3_Pi', 'B3_Pr', 'B3_P2r', 'B3_Omega', 'B3_Incl', 'B3_Xo', 'B3_Yo', 'B3_Zo', 'B3_A1', 'B3_A2', 'B3_A3', 'B3_A4', 'B3_E3', 'B3_E4', 'B3_E5', 'B3_E6', 'B3_E7', 'B3_E8', 'B3_E9', 'B3_E10']
    parname[123:149] = ['B4_No', 'B4_Dz', 'B4_Dr', 'B4_Ro', 'B4_Vi', 'B4_Vr', 'B4_Pi', 'B4_Pr', 'B4_P2r', 'B4_Omega', 'B4_Incl', 'B4_Xo', 'B4_Yo', 'B4_Zo', 'B4_A1', 'B4_A2', 'B4_A3', 'B4_A4', 'B4_E3', 'B4_E4', 'B4_E5', 'B4_E6', 'B4_E7', 'B4_E8', 'B4_E9', 'B4_E10']
    parname[149:183] = ['SR_No', 'SR_R', 'SR_dR', 'SR_dZ', 'LB_No', 'LB_R', 'LB_dR', 'LB_Theta', 'LB_dTheta', 'LB_dZ', 'TB_No', 'TB_R', 'TB_dR', 'TB_Theta', 'TB_dTheta', 'TB_dZ', 'RB_Omega', 'RB_Incl', 'RB_Xo', 'RB_Yo', 'RB_Zo', 'RB_A1', 'RB_A2', 'RB_A3', 'RB_A4', 'RB_E3', 'RB_E4', 'RB_E5', 'RB_E6', 'RB_E7', 'RB_E8', 'RB_E9', 'RB_E10', 'SR_dR2']

    return parname

def get_parindex(parname):
    namelist = get_parnames()
    indx = []
    for name in parname:
        if name in namelist:
            indx.append(namelist.index(name))
        else:
            indx.append(-1)
    return indx

def wiring2tok(wiring):
    ptr1 = 0
    ptr2 = wiring.find('>', ptr1)
    input_token = wiring[ptr1:ptr2].strip()
    output_tokens = wiring[ptr2+1:].split('+')
    output_tokens = [token.strip() for token in output_tokens]
    return input_token, output_tokens

def zpar_wiring(wiring=None, reset=False):
    global wiring_indx
    npar = 256

    if wiring:
        wiring_indx = np.arange(npar)
        for iset in wiring:
            nout = wiring2tok(iset)
            if nout:
                input_token, output_tokens = nout
                indx_in = get_parindex([input_token])[0]
                indx_out = get_parindex(output_tokens)
                for i in range(nout):
                    s = np.where(wiring_indx == indx_out[i])[0]
                    wiring_indx[s] = indx_in
    elif 'wiring_indx' not in globals():
        wiring_indx = np.arange(npar)
    elif reset:
        wiring_indx = np.arange(npar)

    return wiring_indx

def parfix(wl, a, indxpar=None):
    npar = len(a)
    indxa = np.arange(npar)
    if indxpar is None:
        indxpar = indxa

    windx = zpar_wiring()

    pfix = ['FuncIndx', 'B1_Ro', 'B1_Vr', 'B1_P2r', 'B2_Ro', 'B2_Vr', 'B2_P2r', 'B3_Ro', 'B3_Vr', 'B3_P2r', 'B4_Ro', 'B4_Vr', 'B4_P2r', 'LB_Theta', 'TB_Theta', 'RB_Xo', 'RB_Yo', 'RB_Zo']

    if 1.25 not in wl:
        pfix.extend(['C_A1', 'B1_A1', 'B2_A1', 'B3_A1', 'B4_A1', 'RB_A1', 'PF1_C0', 'PF1_C1', 'PF1_Ka'])
    if 2.20 not in wl:
        pfix.extend(['C_A2', 'B1_A2', 'B2_A2', 'B3_A2', 'B4_A2', 'RB_A2', 'PF2_C0', 'PF2_C1', 'PF2_Ka'])
    if 3.50 not in wl:
        pfix.extend(['C_A3', 'B1_A3', 'B2_A3', 'B3_A3', 'B4_A3', 'RB_A3', 'C_E3', 'B1_E3', 'B2_E3', 'B3_E3', 'B4_E3', 'RB_E3', 'PF3_C0', 'PF3_C1', 'PF3_Ka'])
    if 4.90 not in wl:
        pfix.extend(['C_A4', 'B1_A4', 'B2_A4', 'B3_A4', 'B4_A4', 'RB_A4', 'C_E4', 'B1_E4', 'B2_E4', 'B3_E4', 'B4_E4', 'RB_E4'])
    if 12.0 not in wl:
        pfix.extend(['C_E5', 'B1_E5', 'B2_E5', 'B3_E5', 'B4_E5', 'RB_E5'])
    if 25.0 not in wl:
        pfix.extend(['C_E6', 'B1_E6', 'B2_E6', 'B3_E6', 'B4_E6', 'RB_E6'])
    if 60.0 not in wl:
        pfix.extend(['C_E7', 'B1_E7', 'B2_E7', 'B3_E7', 'B4_E7', 'RB_E7'])
    if 100.0 not in wl:
        pfix.extend(['C_E8', 'B1_E8', 'B2_E8', 'B3_E8', 'B4_E8', 'RB_E8'])
    if 140.0 not in wl:
        pfix.extend(['C_E9', 'B1_E9', 'B2_E9', 'B3_E9', 'B4_E9', 'RB_E9'])
    if 240.0 not in wl:
        pfix.extend(['C_E10', 'B1_E10', 'B2_E10', 'B3_E10', 'B4_E10', 'RB_E10'])

    pfixIndx = get_parindex(pfix)

    s = np.where(a == -1)[0]
    if s.size > 0:
        pfixIndx.extend(s)

    s = get_parindex(['K', 'To2'])
    if a[s[0]] == 0.0:
        pfixIndx.extend(s)

    for band in ['B1', 'B2', 'B3', 'B4']:
        s = get_parindex(f'{band}_No')
        if a[s[0]] == 0.0:
            pfixIndx.extend(range(s[0], s[0] + 26))

    for component in ['SR', 'LB', 'TB']:
        s = get_parindex(f'{component}_No')
        if a[s[0]] == 0.0:
            pfixIndx.extend(range(s[0], s[0] + 6))

    s = get_parindex(['SR_No', 'LB_No', 'TB_No'])
    if np.max(a[s]) == 0.0:
        pfixIndx.extend(range(s[0], s[0] + 33))

    pfixIndx = sorted(set(pfixIndx))
    indxpar = sorted(set(indxpar))

    s = np.setdiff1d(indxpar, pfixIndx)
    indxpar = s

    s = np.intersect1d(indxpar, windx)
    if s.size > 0:
        indxpar = s
    else:
        print('Warning: parfix found no free parameters')

    return indxpar

def earthsun(day, lon, lat):
    """
    Returns Solar Elongation and Earth Position from Input Day and Ecliptic Lon/Lat

    Parameters:
    day (float): Day number since Jan 1, 1990 (1 = Jan 1, 1990)
    lon (float): Ecliptic longitude of line-of-sight in degrees
    lat (float): Ecliptic latitude of line-of-sight in degrees

    Returns:
    SolElong (float): Solar elongation in radians
    Earth_Dis (float): Earth-Sun distance in AU
    Earth_Lon (float): Earth's true longitude (radians)
    Earth_Mean_Lon (float): Earth's mean longitude (radians)
    """

    pi2 = 2 * np.pi
    d2r = np.pi / 180.0
    eccen = 0.01671254

    # Solar longitude
    lambda_solar = (-80.598349 + 0.98564736 * day +
                    1.912 * np.cos(pi2 / 365.25 * (day - 94.8))) * d2r

    # Mean anomaly
    mean_anomaly = (356.637087 + 0.98560028 * day) * d2r % pi2

    # Earth-Sun distance (AU)
    Earth_Dis = (1.0 - eccen**2) / (1.0 + eccen * np.cos(mean_anomaly))

    # Earth's true longitude (radians)
    Earth_Lon = (- (np.pi - lambda_solar)) % (2 * np.pi)

    # Solar elongation (radians)
    SolElong = np.arccos(np.cos(lat * d2r) * np.cos(lon * d2r - lambda_solar))

    # Earth's mean longitude (radians)
    Earth_Mean_Lon = (99.403445 + 0.98564736 * day) * d2r % pi2

    return SolElong, Earth_Dis, Earth_Lon, Earth_Mean_Lon

'''
def gaussint(a, b, numpts=50):
    # Placeholder for gaussint function
    # This function should set up the grid and weights for Gauss-Legendre quadrature
    x, w = roots_legendre(numpts)
    x = 0.5 * (b - a) * x + 0.5 * (b + a)
    w = 0.5 * (b - a) * w
    return x, w
'''

def simpint(b, stepsize=0.025):
    """
    Generates Simpson integration abscissae and weights from 0 to b.
    Matches IDL implementation exactly including precision.

    Parameters:
    b (float): Upper limit of integration (lower limit is assumed to be 0)
    stepsize (float): Step size for the integration grid (default: 0.025)

    Returns:
    grid (np.ndarray): Grid of abscissae
    wts (np.ndarray): Weights for Simpson's integration
    """
    h = stepsize
    tsteps = int(round(b/h)) + 1  # total number of grid points

    # CRITICAL: IDL uses SINGLE PRECISION (findgen, fltarr)
    # This must match exactly to get the same numerical results
    grid = np.arange(tsteps, dtype=np.float32) * np.float32(h)
    wts = np.zeros(tsteps, dtype=np.float32)

    # Match IDL's 6-digit truncated constants EXACTLY
    # IDL: .333333, 1.333333, 0.666667
    # Use float32 to match IDL's single precision
    c1 = np.float32(0.333333)
    c2 = np.float32(1.333333)
    c3 = np.float32(0.666667)
    h_f32 = np.float32(h)

    # CRITICAL: IDL sets endpoint AFTER loops to ensure correct value
    wts[0] = c1 * h_f32

    # IDL: for i = 1,tsteps-1,2 do wts(i) = 1.333333 * h
    for i in range(1, tsteps-1, 2):
        wts[i] = c2 * h_f32

    # IDL: for i = 2,tsteps-2,2 do wts(i) = 0.666667 * h
    for i in range(2, tsteps-1, 2):
        wts[i] = c3 * h_f32

    # IDL sets this AFTER the loops (overwrites if tsteps-1 was set in loop)
    wts[tsteps-1] = c1 * h_f32

    return grid, wts


import numpy as np

def scattfunc(phase_type, det, Lambda, LOS, R, Re, SolElong, a, df_out=None, solar_irr=None):
    """
    Compute scattered light contribution of zodiacal dust.

    Parameters:
    - phase_type: str
        'kelsall' or 'skysurf'
    - det: int
        DIRBE detector number (0–9) or -1 if non-DIRBE
    - Lambda: float
        Wavelength in microns
    - LOS: float or np.ndarray
        Line-of-sight distance from Earth (AU)
    - R: float or np.ndarray
        Distance from the Sun to the particle (AU)
    - Re: float
        Distance from Earth to Sun (AU)
    - SolElong: float
        Solar elongation (in radians)
    - a: np.ndarray
        Model parameters (assumed length at least 3 per detector)
    - df_out: np.ndarray or None
        Optional output array to receive partial derivatives [npts, nparams]
    - solar_irr: float or None
        Optional solar irradiance override

    Returns:
    - Scatt: np.ndarray
        Scattered light brightness per LOS element [npts]
    """
    if det <= 2:
        # Solar flux at 1 AU (W/m^2/um) - from IDL code
        SolFlux1AU = np.array([
            2.3405606e+08,  # 1.25 micron
            1.2309874e+08,  # 2.2 micron
            64292872.,      # 3.5 micron
            35733824.,      # 4.9 micron
            5763843.0,      # 12 micron
            1327989.4,      # 25 micron
            230553.73,      # 60 micron
            82999.336,      # 100 micron
            42346.605,      # 140 micron
            14409.608       # 240 micron
        ])

        # Solar flux calculation 
        if solar_irr is not None:
            SolFlux = solar_irr / R**2
        elif det >= 0 and phase_type == 'skysurf':
            SolFlux = solar_sp(Lambda) / R**2
        elif det >= 0 and phase_type == 'kelsall':
            SolFlux = SolFlux1AU[det] / R**2
        else:
            SolFlux = 0.0
        # print(f"SolFlux: {SolFlux}")
        # Phase angle calculation (match IDL logic)
        phase_ang = np.arcsin(np.clip(Re / R * np.sin(SolElong), -1.0, 1.0))
        itest = (LOS >= Re * np.cos(SolElong))
        scat_ang = np.where(itest, np.pi - phase_ang, phase_ang)

        # Phase function parameter selection - match IDL logic
        if len(a) == 9:
            # Original Kelsall: extract 3 parameters based on detector
            aaa = a[3*det:3*det+3]
        elif len(a) % 2 == 0:
            # Hong functions: use all parameters
            aaa = a
        else:
            # Fallback
            aaa = a[0:3]

        # Phase function - use Hong function for SKYSURF, regular phase function for Kelsall
        if phase_type == 'skysurf' and len(a) % 2 == 0:
            # SKYSURF uses Hong phase function
            phase_func = hong_phase_func(scat_ang, aaa)
        else:
            # Use regular Kelsall phase function
            phase_func, _ = phasefunc(scat_ang, aaa)

        phase_func = np.asarray(phase_func)
        # print(f"phase_func: {phase_func}")
# 
        # Scattered light calculation using conditional solar flux
        Scatt = SolFlux * phase_func
        # print(f"Scatt: {Scatt}")
        
        if df_out is not None:
            # Derivatives with respect to phase function parameters
            dphase_da = np.zeros((len(phase), 9))
            for i in range(9):
                da = np.zeros(9)
                da[i] = 1.0
                dphase_da[:, i] = phasefunc(phase, da)
            
            for i in range(9):
                df_out[:, i] = SolFlux1AU[det] * dphase_da[:, i] * r_factor * re_factor
    else:
        Scatt = np.zeros_like(LOS)
        if df_out is not None:
            df_out.fill(0.0)
    
    return Scatt


def thermfunc(det, Lambda, R, a, want_partials=False, no_colcorr=False):
    eps = 1e-20
    To1, To2, Kappa, Delta = a[:4]

    # Conversion factor from W/cm^2/sr to MJy/sr
    cfact = Lambda / 3.0e-10

    # First component
    Temp1 = To1 / R**Delta
    Bnu1, dBdT1 = zPlanck(Temp1, Lambda, return_derivative=True)
    Bnu1 *= cfact
    dBdT1 *= cfact

    if det >= 0 and not no_colcorr:
        CCtherm1 = colcorr(det, Temp1, return_derivative=True)
    else:
        CCtherm1 = 1.0, 0.0

    # Second component
    if Kappa != 0.0:
        Temp2 = To2 / R**Delta
        Bnu2, dBdT2 = zPlanck(Temp2, Lambda, return_derivative=True)
        Bnu2 *= cfact
        dBdT2 *= cfact

        if det >= 0 and not no_colcorr:
            CCtherm2, dCdT2 = colcorr(det, Temp2, return_derivative=True)
        else:
            Temp2 = 0.0
            Bnu2 = eps
            dBdT2 = 0.0
            CCtherm2 = (1.0, 0.0)
            dCdT2 = 0.0
    else:
        Temp2 = 0.0
        Bnu2 = eps
        dBdT2 = 0.0
        CCtherm2 = (1.0, 0.0)
        dCdT2 = 0.0

    # Total thermal brightness
    Therm = Bnu1 * CCtherm1[0] + Kappa * Bnu2 * CCtherm2[0]
    dBdT = dBdT1 + Kappa * dBdT2


    df = None
    if want_partials:
        npts = len(R)
        npar = len(a)
        df = np.zeros((npts, npar), dtype=np.float64)

        dT1 = dBdT1 * CCtherm1[0] + Bnu1 * CCtherm1[1]
        dT2 = (dBdT2 * CCtherm2[0] + Bnu2 * CCtherm2[1]) * Kappa

        df[:, 0] = dT1 / R**Delta                      # To1
        df[:, 1] = dT2 / R**Delta                      # To2
        df[:, 2] = Bnu2 * CCtherm2[0]                     # Kappa
        df[:, 3] = -np.log(R) * (Temp1 * dT1 + Temp2 * dT2)  # Delta

    return Therm, df



def zcloud(x, y, z, R, a, df_out=None, FuncIndx=0):
    """
    Computes the zodiacal cloud number density and optionally its derivatives.

    Parameters:
    - x, y, z: np.ndarray
        Heliocentric Cartesian coordinates
    - R: np.ndarray
        Distance sqrt(x^2 + y^2 + z^2)
    - a: np.ndarray
        Model parameters
    - df_out: np.ndarray or None
        Output for partial derivatives [npts, npar]
    - FuncIndx: int
        Selector for latitudinal density function (0–5)

    Returns:
    - Dens: np.ndarray
        Number density at each point
    - df_out: np.ndarray (if provided)
        Partial derivatives wrt model parameters
    """

    want_partials = df_out is not None
    d2r = np.pi / 180.
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    R = np.asarray(R)
    npts = len(R)

    No = a[0]
    if np.any(No == 0.0):  # Checks if 0.0 is present in the array
        return np.zeros_like(R)


    # Parameters
    Alpha, Beta, Gamma = a[1:4]
    Mu = a[4]
    # print("A:",len(a))  # This will give you the length of `a`
    Mu2, Mu3, Mu4, Mu5, Mu6, Mu7, Mu8, Mu9, Mu10 = a[5:14]
    Omega, Incl = a[14:16] * d2r
    Xo, Yo, Zo = a[16:19]

    # Rotation terms
    sino, coso = np.sin(Omega), np.cos(Omega)
    sini, cosi = np.sin(Incl), np.cos(Incl)

    # Translate to cloud center
    # print(x.shape, y.shape, z.shape)
    # print('x0:', Xo)
    Xp, Yp, Zp = x - Xo, y - Yo, z - Zo
    Rc = np.sqrt(Xp**2+Yp**2+Zp**2)

    # Rotate into cloud coords
    Zc = sino * sini * Xp - coso * sini * Yp + cosi * Zp
    Zeta = np.abs(Zc / Rc)

    # Radial power law (with polynomial mods)
    dR = R - 1.
    AlphaR = Alpha + (Mu6**2 * dR + Mu7**2 * dR**2 + Mu8**2 * dR**3)
    AlphaR = np.minimum(AlphaR, 50.0)  # Prevent overflow

    # Default outputs
    Dens_Cloud_Vert = np.ones_like(R)
    Dens_Cloud_Rad = np.ones_like(R)
    lnR = np.zeros_like(R)

    if FuncIndx == 0:  # John Good
        Dc = np.sqrt(np.clip(Rc**2 - Zc**2, 1e-10, None))
        ZoD = np.abs(Zc / Dc)
        ZoD2G = ZoD**Gamma
        Dens_Cloud_Vert = np.exp(-Beta * ZoD2G)
        Dens_Cloud_Rad = 1. / Dc**AlphaR
        lnR = np.log(Dc)
    elif FuncIndx == 1:  # Modified Fan
        Zeta2G = Zeta**Gamma
        Dens_Cloud_Vert = np.exp(-Beta * Zeta2G)
        Dens_Cloud_Rad = 1. / Rc**AlphaR
        lnR = np.log(Rc)
    elif FuncIndx == 2:  # Widened Modified Fan
        GZR = np.where(Zeta < Mu, 0.5 * Zeta**2 / Mu, Zeta - 0.5 * Mu)
        GZR2G = GZR**Gamma
        Dens_Cloud_Vert = np.exp(-Beta * GZR2G)
        Dens_Cloud_Rad = 1. / Rc**AlphaR
        lnR = np.log(Rc)
    elif FuncIndx == 3:  # Ellipsoid
        Bterm = 1. + (Beta * Zeta)**2
        Dens_Cloud_Vert = 1. / Bterm**Gamma
        Dens_Cloud_Rad = 1. / Rc**AlphaR
        lnR = np.log(Rc)
    elif FuncIndx == 4:  # Sombrero (Cosine)
        cosB = np.sqrt(np.clip(1 - Zeta**2, 1e-10, None))
        cosB2G = cosB**Gamma
        Dens_Cloud_Vert = (1. + Beta * cosB2G) / (1. + Beta)
        Dens_Cloud_Rad = 1. / Rc**AlphaR
        lnR = np.log(Rc)
    elif FuncIndx == 5:  # Polynomial Z-height
        aZ = np.abs(z)
        BetaZ = Beta + Mu3 * aZ + Mu4 * aZ**2 + Mu5 * aZ**3
        Dens_Cloud_Vert = np.exp(-BetaZ * Zeta)
        Dens_Cloud_Rad = 1. / Rc**AlphaR
        lnR = np.log(Rc)

    # Final density
    Dens = No * Dens_Cloud_Vert * Dens_Cloud_Rad
    return Dens #Subject to change

    # If no derivatives requested, return now
    #if not want_partials:
    #    return Dens
'''''
    # Prepare df_out array
    #npar = len(a)
    if df_out.shape != (npts, npar):
        df_out = np.zeros((npts, npar))

    # df wrt No
    df_out[:, 0] = Dens / No

    # Alpha
    df_out[:, 1] = -Dens * lnR

    # Mu6, Mu7, Mu8
    df_out[:, 9] = -Dens * lnR * dR * 2. * Mu6
    df_out[:, 10] = -Dens * lnR * dR**2 * 2. * Mu7
    df_out[:, 11] = -Dens * lnR * dR**3 * 2. * Mu8

    # sZc for directional signs
    sZc = np.sign(Zc)
    dlnf = np.zeros_like(R)

    if FuncIndx == 0:
        df_out[:, 2] = -Dens * ZoD2G
        df_out[:, 3] = -Beta * Dens * ZoD2G * np.log(np.clip(ZoD, 1e-10, None))
        Zsqr = np.clip(1. - Zeta**2, 1e-10, None)
        dlnf = AlphaR * Zeta / Zsqr - Beta * Gamma * Zeta**(Gamma - 1) * Zsqr**(-Gamma / 2 - 1)

    elif FuncIndx == 1:
        df_out[:, 2] = -Dens * Zeta2G
        df_out[:, 3] = -Beta * Dens * Zeta2G * np.log(np.clip(Zeta, 1e-10, None))
        dlnf = -Beta * Gamma * Zeta**(Gamma - 1)

    elif FuncIndx == 2:
        df_out[:, 2] = -Dens * GZR2G
        df_out[:, 3] = -Beta * Dens * GZR2G * np.log(np.clip(GZR, 1e-10, None))
        df_out[:, 4] = Beta * Dens * 0.5 * Gamma * GZR**(Gamma - 1) * np.where(Zeta < Mu, (Zeta / Mu)**2, 1.0)
        dlnf = -Beta * Gamma * GZR**(Gamma - 1) * np.where(Zeta < Mu, Zeta / Mu, 1.0)

    elif FuncIndx == 3:
        Bterm = 1. + (Beta * Zeta)**2
        df_out[:, 2] = -2. * Gamma / Beta * Dens * (1. - 1. / Bterm)
        df_out[:, 3] = -Dens * np.log(Bterm)
        dlnf = -2 * Gamma * Beta**2 * Zeta / Bterm

    elif FuncIndx == 4:
        cosB = np.sqrt(np.clip(1 - Zeta**2, 1e-10, None))
        cosB2G = cosB**Gamma
        df_out[:, 2] = Dens / Dens_Cloud_Vert * (-Dens_Cloud_Vert + cosB2G) / (1 + Beta)
        df_out[:, 3] = Dens / Dens_Cloud_Vert * (Beta * cosB2G / (1 + Beta) * np.log(cosB))
        dlnf = -Gamma * Zeta * Beta / (1 + Beta) / Dens_Cloud_Vert * cosB2G / cosB**2

    elif FuncIndx == 5:
        aZ = np.abs(z)
        BetaZ = Beta + Mu3 * aZ + Mu4 * aZ**2 + Mu5 * aZ**3
        df_out[:, 2] = -Dens * Zeta
        df_out[:, 6] = -Dens * Zeta * aZ
        df_out[:, 7] = -Dens * Zeta * aZ**2
        df_out[:, 8] = -Dens * Zeta * aZ**3
        dlnf = -BetaZ

    # Geometry partials (omega, inclination, Xo, Yo, Zo)
    dZc_dX = sino * sini
    dZc_dY = -coso * sini
    dZc_dZ = cosi

    df_out[:, 14] = Dens * dlnf * sZc / Rc * (coso * sini * Xp + sino * sini * Yp) * d2r
    df_out[:, 15] = Dens * dlnf * sZc / Rc * (sino * cosi * Xp - coso * cosi * Yp - sini * Zp) * d2r
    df_out[:, 16] = Dens * (-dlnf * sZc / Rc * dZc_dX + (dlnf * Zeta + AlphaR) * Xp / Rc**2)
    df_out[:, 17] = Dens * (-dlnf * sZc / Rc * dZc_dY + (dlnf * Zeta + AlphaR) * Yp / Rc**2)
    df_out[:, 18] = Dens * (-dlnf * sZc / Rc * dZc_dZ + (dlnf * Zeta + AlphaR) * Zp / Rc**2)

    return Dens, df_out
'''''


def new_isocloud(x, y, z, R, a):
    """
    Isotropic zodiacal cloud component.
    Uses Widened Modified Fan (FuncIndx=2) with Beta=0, Gamma=0.

    Parameters:
    - x, y, z, R: Heliocentric coordinates and distance
    - a: Cloud parameters (uses Mu6, Mu7, Mu8 for radial correction)

    Returns: Isotropic component density
    """
    # Hardcoded parameters
    No = 1.356e-9  # Normalization
    Alpha = -1.46  # Radial power law
    Beta = 0.0  # No latitude variation
    Gamma = 0.0  # No latitude variation

    # Hardcoded geometry (centered on Sun)
    Omega = 0.0
    Incl = 0.0
    Xo = 0.0
    Yo = 0.0
    Zo = 0.0

    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    R = np.asarray(R)

    # Extract Mu parameters for radial correction
    Mu6 = a[9]
    Mu7 = a[10]
    Mu8 = a[11]

    # Compute radial power-law with polynomial correction
    dR = R - 1.0
    AlphaR = Alpha + (Mu6**2 * dR + Mu7**2 * dR**2 + Mu8**2 * dR**3)
    AlphaR = np.minimum(AlphaR, 50.0)  # Prevent overflow

    # Since Omega=Incl=0 and Xo=Yo=Zo=0, Rc = R
    # With Beta=0, Gamma=0: Dens_Cloud_Vert = 1
    # Final density: pure radial power law
    Dens = No / R**AlphaR

    return Dens


def migband(x, y, z, R, a, want_partials=False):
    """
    Computes number density of a migrating dust band.

    Parameters
    ----------
    x, y, z : ndarray
        Heliocentric Cartesian coordinates.
    R : ndarray
        Radial distance sqrt(x^2 + y^2 + z^2).
    a : ndarray
        Model parameters (length 14 expected).
    want_partials : bool, optional
        If True, computes partial derivatives.

    Returns
    -------
    dens : ndarray
        Number density.
    df : ndarray
        Partial derivatives if want_partials is True, else None.
    """

    No = a[0]
    df = None

    if No != 0.0:
        d2r = np.pi / 180.0

        Dz      = a[1] * d2r
        Dr      = a[2]
        Ro      = a[3]
        Vi      = a[4]
        Vr      = a[5]  # Not used
        Pi      = a[6]
        Pr      = a[7]
        P2r     = a[8]  # Not used
        Omega   = a[9]  * d2r
        Incl    = a[10] * d2r
        Xo      = a[11]
        Yo      = a[12]
        Zo      = a[13]

        # Geometry
        Xp = x - Xo
        Yp = y - Yo
        Zp = z - Zo

        sino = np.sin(Omega)
        coso = np.cos(Omega)
        sini = np.sin(Incl)
        cosi = np.cos(Incl)

        Zc = sino * sini * Xp - coso * sini * Yp + cosi * Zp
        Rc = np.sqrt(Xp**2 + Yp**2 + Zp**2)

        Zeta = np.abs(Zc / Rc)
        ZDz = Zeta / Dz
        RDr = Rc / Dr
        ViTerm = 1.0 + (ZDz ** Pi) / Vi

        # Density
        Dens = np.zeros_like(Rc)
        WtTerm = np.ones_like(RDr)

        arg1 = RDr ** 20
        arg2 = ZDz ** 6

        valid = np.where((arg1 <= 86) & (arg2 <= 86))[0]
        if valid.size > 0:
            WtTerm[valid] = 1.0 - np.exp(-arg1[valid])
            Dens[valid] = No * (Ro / Rc[valid])**Pr * np.exp(-arg2[valid]) * ViTerm[valid] * WtTerm[valid]

        valid = np.where((arg2 <= 86) & (arg1 > 86))[0]
        if valid.size > 0:
            Dens[valid] = No * (Ro / Rc[valid])**Pr * np.exp(-arg2[valid]) * ViTerm[valid]

        if want_partials:
            npts = x.size
            npar = len(a)
            df = np.zeros((npts, npar), dtype=np.float64)

            sZ = np.sign(Zc)
            dlnf = (-6.0 * ZDz**5 + Pi / Vi * ZDz**(Pi - 1) / ViTerm) / Dz

            # df with respect to each parameter
            df[:, 0] = Dens / No  # No
            df[:, 1] = -Dens * ZDz * dlnf * d2r  # Dz

            # Dr
            df[:, 2] = 0.0
            valid = np.where((arg2 <= 86) & (arg1 <= 86) & (Dens != 0.0))[0]
            if valid.size > 0:
                df[valid, 2] = No * (Ro / Rc[valid])**Pr * np.exp(-arg2[valid]) * ViTerm[valid] * \
                               (-20 * RDr[valid]**20 * np.exp(-arg1[valid]))

            # Vi
            df[:, 4] = -Dens * (ZDz**Pi) / (Vi**2 * ViTerm)

            # Pi
            df[:, 6] = Dens * (ZDz**Pi) * np.log(ZDz) / (ViTerm * Vi)

            # Pr
            df[:, 7] = Dens * np.log(Ro / Rc)

            # Omega
            df[:, 9] = Dens * dlnf * sZ / Rc * (coso * sini * Xp + sino * sini * Yp) * d2r

            # Incl
            df[:, 10] = Dens * dlnf * sZ / Rc * (sino * cosi * Xp - coso * cosi * Yp - sini * Zp) * d2r

            # Xo
            df[:, 11] = Dens * (-dlnf * sZ / Rc * sini * sino + (dlnf * Zeta + Pr) * Xp / Rc**2)

            # Yo
            df[:, 12] = Dens * (dlnf * sZ / Rc * coso * sini + (dlnf * Zeta + Pr) * Yp / Rc**2)

            # Zo
            df[:, 13] = Dens * (-dlnf * sZ / Rc * cosi + (dlnf * Zeta + Pr) * Zp / Rc**2)

        return Dens, df

    else:
        return np.zeros_like(x), None


import numpy as np

def solring(x, y, z, R, Theta, a, want_partials=False):
    """
    Computes number density of Earth's resonant dust ring.

    Parameters
    ----------
    x, y, z : ndarray
        Heliocentric Cartesian coordinates.
    R : ndarray
        Radial distance (sqrt(x^2 + y^2 + z^2)).
    Theta : float or ndarray
        Earth's mean longitude in radians.
    a : ndarray
        Model parameters (length 21 expected).
    want_partials : bool
        If True, returns partial derivatives.

    Returns
    -------
    Dens : ndarray
        Number density.
    df : ndarray or None
        Derivatives with respect to parameters (if want_partials).
    """

    d2r = np.pi / 180.0

    SR_No = a[0]
    LB_No = a[4]
    TB_No = a[10]

    Dens_SR = Dens_LB = Dens_TB = np.zeros_like(R)
    df = None

    if SR_No != 0.0 or LB_No != 0.0 or TB_No != 0.0:

        # Ring parameters
        SR_R = a[1]
        SR_dR = a[2]
        SR_dZ = a[3]

        # Leading blob
        LB_R = a[5]
        LB_dR = a[6]
        LB_Theta = a[7] * d2r
        LB_dTheta = a[8] * d2r
        LB_dZ = a[9]

        # Trailing blob
        TB_R = a[11]
        TB_dR = a[12]
        TB_Theta = a[13] * d2r
        TB_dTheta = a[14] * d2r
        TB_dZ = a[15]

        # Symmetry plane
        Omega = a[16] * d2r
        Incl = a[17] * d2r
        Xo, Yo, Zo = a[18], a[19], a[20]

        # Rotation trigs
        sino, coso = np.sin(Omega), np.cos(Omega)
        sini, cosi = np.sin(Incl), np.cos(Incl)

        # Geometry: Z height in rotated plane
        SR_Z = sino * sini * x - coso * sini * y + cosi * z
        # Theta offset from Earth
        dTheta = np.arctan2(y, x) - Theta
        # --- Dust Ring Density ---
        if SR_No != 0.0:
            arg = (R - SR_R)**2 / SR_dR**2 + np.abs(SR_Z) / SR_dZ
            valid = np.where(arg <= 86)
            Dens_SR = np.zeros_like(R)
            if valid[0].size > 0:
                Dens_SR[valid] = SR_No * np.exp(-arg[valid])

        # --- Leading Blob ---
        if LB_No != 0.0:
            LB_Delta = (dTheta - LB_Theta) % (2 * np.pi)
            LB_Delta += 2 * np.pi * (np.where(LB_Delta < -np.pi, 1, 0) - np.where(LB_Delta > np.pi, 1, 0))

            arg = (R - LB_R)**2 / LB_dR**2 + LB_Delta**2 / LB_dTheta**2 + np.abs(SR_Z) / LB_dZ
            valid = np.where(arg <= 86)
            Dens_LB = np.zeros_like(R)
            if valid[0].size > 0:
                Dens_LB[valid] = LB_No * np.exp(-arg[valid])

        # --- Trailing Blob ---
        if TB_No != 0.0:
            TB_Delta = (dTheta - TB_Theta) % (2 * np.pi)
            TB_Delta += 2 * np.pi * (np.where(TB_Delta < -np.pi, 1, 0) - np.where(TB_Delta > np.pi, 1, 0))

            arg = (R - TB_R)**2 / TB_dR**2 + TB_Delta**2 / TB_dTheta**2 + np.abs(SR_Z) / TB_dZ
            valid = np.where(arg <= 86)
            Dens_TB = np.zeros_like(R)
            if valid[0].size > 0:
                Dens_TB[valid] = TB_No * np.exp(-arg[valid])

        # --- Total Density ---
        Dens = Dens_SR + Dens_LB + Dens_TB
        # --- Partials ---
        if want_partials:
            npts = len(x)
            df = np.zeros((npts, len(a)), dtype=np.float64)

            if SR_No != 0.0:
                df[:, 0] = Dens_SR / SR_No  # SR_No
                df[:, 1] = Dens_SR * 2 * (R - SR_R) / SR_dR**2  # SR_R
                df[:, 2] = Dens_SR * 2 * (R - SR_R)**2 / SR_dR**3  # SR_dR
                df[:, 3] = Dens_SR * np.abs(SR_Z) / SR_dZ**2  # SR_dZ

            if LB_No != 0.0:
                df[:, 4] = Dens_LB / LB_No  # LB_No
                df[:, 5] = Dens_LB * 2 * (R - LB_R) / LB_dR**2  # LB_R
                df[:, 6] = Dens_LB * 2 * (R - LB_R)**2 / LB_dR**3  # LB_dR
                df[:, 8] = Dens_LB * 2 * LB_Delta**2 / LB_dTheta**3 * d2r  # LB_dTheta
                df[:, 9] = Dens_LB * np.abs(SR_Z) / LB_dZ**2  # LB_dZ

            if TB_No != 0.0:
                df[:, 10] = Dens_TB / TB_No  # TB_No
                df[:, 11] = Dens_TB * 2 * (R - TB_R) / TB_dR**2  # TB_R
                df[:, 12] = Dens_TB * 2 * (R - TB_R)**2 / TB_dR**3  # TB_dR
                df[:, 14] = Dens_TB * 2 * TB_Delta**2 / TB_dTheta**3 * d2r  # TB_dTheta
                df[:, 15] = Dens_TB * np.abs(SR_Z) / TB_dZ**2  # TB_dZ

            # Common for ring + blobs
            signZ = np.sign(SR_Z)
            dfdz = -signZ * (Dens_SR / SR_dZ + Dens_LB / LB_dZ + Dens_TB / TB_dZ)

            df[:, 16] = dfdz * (coso * sini * x + sino * sini * y) * d2r  # Omega
            df[:, 17] = dfdz * (sino * cosi * x - coso * cosi * y - sini * z) * d2r  # Incl

        return Dens, df

    else:
        return np.zeros_like(R), None


import numpy as np

def zsrcfunc(det, Scatt, Therm, a, df_out=None, phase_type='kelsall'):
    """
    Computes the zodiacal dust source function (brightness per particle).

    Parameters:
    - det: int
        DIRBE detector number (0–9) or -1 if non-DIRBE
    - Scatt: float or np.ndarray
        Scattered light contribution
    - Therm: float or np.ndarray
        Thermal emission contribution
    - a: np.ndarray
        Model parameters [Albedo 0–2, shared Albedo (3), Emissivities 4–11]
    - df_out: np.ndarray or None
        Optional output array for derivatives [npts, nparams]

    Returns:
    - Source: float or np.ndarray
        Brightness per particle (MJy/sr)
    - df_out (if provided): np.ndarray
        Partial derivatives w.r.t. parameters
    - dScatt (if provided): float or np.ndarray
        Derivative of source function w.r.t. scattered light
    - dTherm (if provided): float or np.ndarray
        Derivative w.r.t. thermal emission
    """

    want_partials = df_out is not None
    a = np.asarray(a)
    npar = len(a)

    # Albedo assignment (0–3 per detector, rest share a[3])
    AlbedoDet = np.concatenate([a[0:4], np.full(6, a[3])])

    # Emissivity assignment (detectors 0–1: emiss = 1.0, rest from a[4:11])
    EmissDet = np.concatenate([np.ones(2), a[4:12]])

    # Default if det is invalid
    if det >= 0:
        Albedo = AlbedoDet[det]
        # For skysurf phase type, use emissivity from detector 2 slot (3.5 micron)
        if phase_type == 'skysurf':
            Emiss = EmissDet[2]  # Use slot for 3.5 micron detector
        elif phase_type == 'kelsall':
            Emiss = EmissDet[det]
        else:
            Emiss = EmissDet[det]
    else:
        Albedo = 0.0
        Emiss = 1.0


    # Total brightness
    Source = Albedo * Scatt + Emiss * (1.0 - Albedo) * Therm

    dScatt = Albedo
    dTherm = Emiss * (1.0 - Albedo)

    if want_partials:
        npts = len(np.atleast_1d(Therm))
        if df_out.shape != (npts, npar):
            df_out = np.zeros((npts, npar), dtype=float)
        else:
            df_out.fill(0.0)

        # Derivative w.r.t. albedo
        dA = Scatt - Emiss * Therm
        if det < 3:
            df_out[:, det] = dA
        elif det >= 3:
            df_out[:, 3] = dA

        # Derivative w.r.t. emissivity
        if det >= 2:
            df_out[:, 4 + (det - 2)] = Therm / Emiss

        return Source, df_out, dScatt, dTherm

    return Source


def zkernel(data, a, phase_type='kelsall', indxpar=None, losinfo=False, no_colcorr=False, dbwave=None, solar_irr=None, new_iso_comp=False, iso_comp_only=False):
    want_los_info = losinfo
    # print('a:',a)
    RMAX = 5.2
    NUMBER_OF_STEPS = 50

    npts = len(data)
    d2r = np.pi / 180.0
    eps = 1e-20
    if dbwave is None:
        dbwave = [1.25, 2.2, 3.5, 4.9, 12., 25., 60., 100., 140., 240.]
    nmsg = 50000

    windx = zpar_wiring()

    # Apply wiring to zpars (critical for Kelsall mode)
    # This reindexes the zpars array to wire parameters together
    # For example: C_E3>C_A3 makes emissivity equal to albedo at 3.5 microns
    # Note: When no wiring is specified, this is just identity mapping
    a = a[windx]


    FuncIndx = a[0]

    nScatt = 9
    iScatt = np.arange(1, nScatt + 1)
    aScatt = a[iScatt]

    # Check if Hong parameters are stored at index 183+ (SKYSURF mode)
    # Only use Hong params if we're in skysurf mode AND they're non-zero
    if phase_type == 'skysurf' and len(a) > 183:
        hong_params = a[183:189]  # Get exactly 6 Hong parameters
        if np.any(hong_params != 0):  # Only use if non-zero
            aScatt = hong_params

    nTherm = 4
    iTherm = np.arange(10, 10 + nTherm)
    aTherm = a[iTherm]

    nDens_C = 19
    iDens_C = np.arange(14, 14 + nDens_C)
    aDens_C = a[iDens_C]
    nSrc_C = 12
    iSrc_C = np.arange(33, 33 + nSrc_C)
    aSrc_C = a[iSrc_C]

    nDens_B1 = 14
    iDens_B1 = np.arange(45, 45 + nDens_B1)
    aDens_B1 = a[iDens_B1]
    nSrc_B1 = 12
    iSrc_B1 = np.arange(59, 59 + nSrc_B1)
    aSrc_B1 = a[iSrc_B1]

    nDens_B2 = 14
    iDens_B2 = np.arange(71, 71 + nDens_B2)
    aDens_B2 = a[iDens_B2]
    nSrc_B2 = 12
    iSrc_B2 = np.arange(85, 85 + nSrc_B2)
    aSrc_B2 = a[iSrc_B2]

    nDens_B3 = 14
    iDens_B3 = np.arange(97, 97 + nDens_B3)
    aDens_B3 = a[iDens_B3]
    nSrc_B3 = 12
    iSrc_B3 = np.arange(111, 111 + nSrc_B3)
    aSrc_B3 = a[iSrc_B3]

    nDens_B4 = 14
    iDens_B4 = np.arange(123, 123 + nDens_B4)
    aDens_B4 = a[iDens_B4]
    nSrc_B4 = 12
    iSrc_B4 = np.arange(137, 137 + nSrc_B4)
    aSrc_B4 = a[iSrc_B4]

    nDens_RB = 21
    iDens_RB = np.arange(149, 149 + nDens_RB)
    iDens_RB = np.append(iDens_RB, 182)
    aDens_RB = a[iDens_RB]
    nSrc_RB = 12
    iSrc_RB = np.arange(170, 170 + nSrc_RB)
    aSrc_RB = a[iSrc_RB]

    detnum = np.zeros(npts, dtype=int) - 1
    for i in range(10):
        detnum += (np.array([d['wave_len'] for d in data]) == dbwave[i]) * (i + 1)

    SolElong, Earth_Dis, Earth_Lon, Earth_Mean_Lon = earthsun(data['day1990'], data['longitude'], data['latitude'])

    f = np.zeros(npts)

    if indxpar is not None:
        want_partials = True
        npar = len(a)
        if indxpar is None:
            nterms = npar
            indxpar = np.arange(nterms)
        else:
            nterms = len(indxpar)

        df = np.zeros((npts, nterms))
        parnum = np.arange(npar) - 1
        parnum[indxpar] = np.arange(nterms)
    else:
        want_partials = False

    if want_los_info:
        losdata = []

    # print('aSrc_C:', aSrc_C)
    # print('aDens_C:', aDens_C)
    for ilos in range(npts):
        lat = data[ilos]['latitude'] * d2r
        lon = data[ilos]['longitude'] * d2r

        Re = Earth_Dis[ilos]
        Theta = Earth_Lon[ilos]

        COSTHETA = np.cos(Theta)
        SINTHETA = np.sin(Theta)
        COSLON = np.cos(lon)
        SINLON = np.sin(lon)
        COSLAT = np.cos(lat)
        SINLAT = np.sin(lat)

        X0 = Re * COSTHETA
        Y0 = Re * SINTHETA
        B = 2.0 * (X0 * COSLAT * COSLON + Y0 * COSLAT * SINLON)
        C = Re**2 - RMAX**2
        Q = -0.5 * B * (1.0 + np.sqrt(B**2 - 4.0 * C) / abs(B))
        RANGE = max([Q, C / Q])

        # CRITICAL FIX: IDL code uses simpint for ALL cases (gaussint is commented out)
        # The original Python code had conditional logic, but IDL always uses Simpson integration
        los, gqwts = simpint(RANGE)

        Sxy = los * COSLAT
        X = Re * COSTHETA + Sxy * COSLON
        Y = Re * SINTHETA + Sxy * SINLON
        Z = los * SINLAT
        R = np.sqrt(X**2 + Y**2 + Z**2)

        Lambda = data[ilos]['wave_len']
        # Detector assignment logic - match IDL implementation exactly
        if phase_type == 'skysurf':
            Det = 0
        else:
            Det = detnum[ilos]

        if Det >= 0:
            Scatt = scattfunc(phase_type, Det, Lambda, los, R, Re, SolElong[ilos], aScatt, solar_irr=solar_irr)
            # print(f"Scatt: {Scatt}")
        else:
            Scatt = 0.0
        Therm, _ = thermfunc(Det, Lambda, R, aTherm, no_colcorr=no_colcorr)
        Dens_C = zcloud(X, Y, Z, R, aDens_C, FuncIndx=FuncIndx)
        # Source function
        Albedo = aSrc_C[0] if Det == 0 else (aSrc_C[1] if Det == 1 else (aSrc_C[2] if Det == 2 else 0.0))
        Emiss = aSrc_C[4] if Det == 0 else (aSrc_C[5] if Det == 1 else (aSrc_C[6] if Det == 2 else 1.0))
        Source = Albedo * Scatt + Emiss * (1.0 - Albedo) * Therm
        Src_C= zsrcfunc(Det, Scatt, Therm, aSrc_C, phase_type=phase_type)

        # Isotropic component
        if new_iso_comp:
            Dens_new = new_isocloud(X, Y, Z, R, aDens_C)
        else:
            Dens_new = np.zeros_like(R)

        Dens_B1, dDens_B1 = migband(X, Y, Z, R, aDens_B1)
        Src_B1 = zsrcfunc(Det, Scatt, Therm, aSrc_B1, phase_type=phase_type)

        Dens_B2, dDens_B2 = migband(X, Y, Z, R, aDens_B2)
        Src_B2 = zsrcfunc(Det, Scatt, Therm, aSrc_B2, phase_type=phase_type)

        Dens_B3, dDens_B3 = migband(X, Y, Z, R, aDens_B3)
        Src_B3 = zsrcfunc(Det, Scatt, Therm, aSrc_B3, phase_type=phase_type)

        Dens_B4, dDens_B4 = migband(X, Y, Z, R, aDens_B4)
        Src_B4 = zsrcfunc(Det, Scatt, Therm, aSrc_B4, phase_type=phase_type)

        Dens_RB, dDens_RB = solring(X, Y, Z, R, Earth_Mean_Lon[ilos], aDens_RB)
        Src_RB = zsrcfunc(Det, Scatt, Therm, aSrc_RB, phase_type=phase_type)

        # Isotropic component only mode
        if iso_comp_only:
            Dens_C = np.zeros_like(R)
            Dens_B1 = np.zeros_like(R)
            Dens_B2 = np.zeros_like(R)
            Dens_B3 = np.zeros_like(R)
            Dens_B4 = np.zeros_like(R)
            Dens_RB = np.zeros_like(R)

        Flux = (Src_C * Dens_C + Src_B1 * Dens_B1 + Src_B2 * Dens_B2 + Src_B3 * Dens_B3 + Src_B4 * Dens_B4 + Src_RB * Dens_RB + Src_C * Dens_new)

        # Sum the flux along the line of sight
        f[ilos] = np.sum(gqwts * Flux)

        if want_los_info:
            Dens = Dens_C + Dens_B1 + Dens_B2 + Dens_B3 + Dens_RB + Dens_new
            losdata.append({
                'wave_len': Lambda,
                'pixel_no': data[ilos]['pixel_no'],
                'lon': lon / d2r,
                'lat': lat / d2r,
                'elong': SolElong[ilos] / d2r,
                'day': data[ilos]['day1990'],
                'earth_lon': Theta / d2r,
                'earth_rad': Re,
                'element': {
                    'earth_dist': los,
                    'solar_dist': R,
                    'ecliptic_z': Z,
                    'density': Dens,
                    'flux': Flux,
                    'int_wts': gqwts
                },
                'zodi': f[ilos]
            })

        if want_partials:
            for ii in range(nScatt):
                ipar = parnum[windx[iScatt[ii]]]
                if ipar >= 0:
                    dSrcS = dScatt[:, ii] * (
                        dSrc_dScatt_C * Dens_C +
                        dSrc_dScatt_B1 * Dens_B1 +
                        dSrc_dScatt_B2 * Dens_B2 +
                        dSrc_dScatt_B3 * Dens_B3 +
                        dSrc_dScatt_B4 * Dens_B4 +
                        dSrc_dScatt_RB * Dens_RB
                    )
                    df[ilos, ipar] += np.sum(gqwts * dSrcS)

            for ii in range(nTherm):
                ipar = parnum[windx[iTherm[ii]]]
                if ipar >= 0:
                    dSrcT = dTherm[:, ii] * (
                        dSrc_dTherm_C * Dens_C +
                        dSrc_dTherm_B1 * Dens_B1 +
                        dSrc_dTherm_B2 * Dens_B2 +
                        dSrc_dTherm_B3 * Dens_B3 +
                        dSrc_dTherm_B4 * Dens_B4 +
                        dSrc_dTherm_RB * Dens_RB
                    )
                    df[ilos, ipar] += np.sum(gqwts * dSrcT)

            if aDens_C[0] != 0.0:
                for ii in range(nDens_C):
                    ipar = parnum[windx[iDens_C[ii]]]
                    if ipar >= 0:
                        df[ilos, ipar] += np.sum(gqwts * Src_C * dDens_C[:, ii])
                for ii in range(nSrc_C):
                    ipar = parnum[windx[iSrc_C[ii]]]
                    if ipar >= 0:
                        df[ilos, ipar] += np.sum(gqwts * dSrc_C[:, ii] * Dens_C)

            if aDens_B1[0] != 0.0:
                for ii in range(nDens_B1):
                    ipar = parnum[windx[iDens_B1[ii]]]
                    if ipar >= 0:
                        df[ilos, ipar] += np.sum(gqwts * Src_B1 * dDens_B1[:, ii])
                for ii in range(nSrc_B1):
                    ipar = parnum[windx[iSrc_B1[ii]]]
                    if ipar >= 0:
                        df[ilos, ipar] += np.sum(gqwts * dSrc_B1[:, ii] * Dens_B1)

            if aDens_B2[0] != 0.0:
                for ii in range(nDens_B2):
                    ipar = parnum[windx[iDens_B2[ii]]]
                    if ipar >= 0:
                        df[ilos, ipar] += np.sum(gqwts * Src_B2 * dDens_B2[:, ii])
                for ii in range(nSrc_B2):
                    ipar = parnum[windx[iSrc_B2[ii]]]
                    if ipar >= 0:
                        df[ilos, ipar] += np.sum(gqwts * dSrc_B2[:, ii] * Dens_B2)

            if aDens_B3[0] != 0.0:
                for ii in range(nDens_B3):
                    ipar = parnum[windx[iDens_B3[ii]]]
                    if ipar >= 0:
                        df[ilos, ipar] += np.sum(gqwts * Src_B3 * dDens_B3[:, ii])
                for ii in range(nSrc_B3):
                    ipar = parnum[windx[iDens_B3[ii]]]
                    if ipar >= 0:
                        df[ilos, ipar] += np.sum(gqwts * dSrc_B3[:, ii] * Dens_B3)

            if aDens_B4[0] != 0.0:
                for ii in range(nDens_B4):
                    ipar = parnum[windx[iDens_B4[ii]]]
                    if ipar >= 0:
                        df[ilos, ipar] += np.sum(gqwts * Src_B4 * dDens_B4[:, ii])
                for ii in range(nSrc_B4):
                    ipar = parnum[windx[iSrc_B4[ii]]]
                    if ipar >= 0:
                        df[ilos, ipar] += np.sum(gqwts * dSrc_B4[:, ii] * Dens_B4)

            if aDens_RB[0] != 0.0:
                for ii in range(nDens_RB):
                    ipar = parnum[windx[iDens_RB[ii]]]
                    if ipar >= 0:
                        df[ilos, ipar] += np.sum(gqwts * Src_RB * dDens_RB[:, ii])
                for ii in range(nSrc_RB):
                    ipar = parnum[windx[iSrc_RB[ii]]]
                    if ipar >= 0:
                        df[ilos, ipar] += np.sum(gqwts * dSrc_RB[:, ii] * Dens_RB)

        # if (ilos + 1) % nmsg == 0:
        #     print(f'LOS # {ilos}')

    if want_los_info:
        return losdata

    if want_partials:
        return f, df

    return f


