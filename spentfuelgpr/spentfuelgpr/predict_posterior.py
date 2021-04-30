
import json
import numpy
import os

from kernel import Kernel

def main():
    return

def predict():
    """Calculate the spent fuel composition."""
    # Check if the kernels needed exist.
    # Rationale for the choice of these isotopes:
    # - isotope fraction > 1e-15
    # - stable U isotopes, i.e., 1/2 time > days (omit U230, 231, 237)
    # - stable Pu isotopes
    # - decay into 'interesting' material (relevant for U->Np->Pu)
    # - stable Pu isotopes
    isotopes = ('U232', 'U233', 'U234', 'U235', 'U235m', 'U236', 'U238',
                'U239', 'U240', 'Pu238', 'Pu239', 'Pu240', 'Pu241', 
                'Pu242', 'Pu243', 'Pu244', 'Np239', 'Np240', 'Np240m', 
                'Np241')
    data_dir = os.path.join(os.path.split(__file__)[0], "..", "data")
    kernel_dir = os.path.join(data_dir, "trained_kernels")
    if not os.path.isdir(kernel_dir):
        raise OSError("'trained_kernels' directory not found!")
    for iso in isotopes:
        fname = iso + ".npy"
        if not os.path.isfile(os.path.join(kernel_dir, fname)):
            msg = f"Trained kernel '{fname}' not found!"
            raise FileNotFoundError(msg)

    # Load input parameters, training data and the kernel type.
    kernel_type = "ASQE"
    par = ("enrichment", "temperature", "power_output", "burnup")
    reactor_input_params = get_input_params(par)
    size = 100  # Needed if a subset of the data was used for training.
    training_data = np.load(os.path.join(data_dir, "x_trainingset.npy"),
                            allow_pickle=True)[:size]
    spent_fuel_composition = {}
    for iso in isotopes:
        fname = f"{iso}.npy"
        trained_kernel = np.load(os.path.join(kernel_dir, fname), 
                                 allow_pickle=True).item()
        mass = run_kernel(reactor_input_params, training_data,
                          trained_kernel, kernel_type)
        spent_fuel_composition[iso] = mass

    store_results(spent_fuel_composition)

def get_input_params(pnames):
    """Read in the Gpr input parameters.

    Currently, the input filename is hardcoded. Maybe this will be 
    changed in the future but at the moment I have no idea how this 
    could be done.
    
    Parameters
    ----------
    pnames : tuple
        The names of the parameters to be extracted from the JSON 
        file.

    Returns
    -------
    input_params : dict
        A dictionary with the extracted values and with keys being the 
        names from `pnames`.
    """
    data_dir = os.path.join(os.path.split(__file__)[0], "..", "data")
    fname = "gpr_reactor_input_params.json"
    with open(os.path.join(data_dir, fname), "r") as f:
        data = json.load(f)
        input_params = {}
        
        for param in pnames:
            param = param.lower()
            if param == "enrichment":
                enrich = data["fresh_fuel_composition"]["922350000"]
                input_params["enrichment"] = enrich
            else:
                try:
                    input_params[param] = data[param]
                except KeyError as e:
                    msg = (f"Key '{param}' not found. Enter a value "
                           + "manually to continue or press Enter, in "
                           + "which case an error is raised.\n")
                    manual_entry = input(msg)
                    if manual_entry == "":
                        raise e
                    input_params[param] = manual_entry
    
    return input_params

def store_results(composition):
    """Store the composition in a .json file.

    Parameters
    ----------
    composition : dict
        A dictionary with the keys being the isotopes and the values
        being the corresponding masses.
    """
    data_dir = os.path.join(os.path.split(__file__)[0], "..", "data")
    fname = "gpr_reactor_spent_fuel_composition.json"
    with open(os.path.join(data_dir, fname), "r") as f:
        json.dump(composition, f, indent=2)

def run_kernel(reactor_input_params, x_train, trained_kernel, 
               kernel_type='ASQE'):
    """Calculate the amount of the isotope in question in the spent fuel

    Parameters
    ----------
    reactor_input_params
        The input parameters used in the calculations and as defined
        above in *TODO insert name of the function calling `prediction`*
    x_train
        The input parameters used during training.
    trained_kernel
        The trained kernel, specific to the isotope in question
    kernel_type
        The type of the trained kernel. TODO UPDATE WHICH KERNEL IS USED

    Returns
    -------
    mu_s
        The (mean predicted) mass of the isotope in the spent fuel
        after one irradiation period.
    """
    if (kernel_type != 'ASQE'):
        msg = "Currently, only the 'ASQE' kernel type is supported."
        raise ValueError(msg)
    
    k_s = Kernel(reactor_input_params, x_train, kernel_type, params, 
                 gradient=False)
    mu_s = np.dot(k_s.T, alpha_)
    
    return mu_s[0]

"""
Part of the spent fuel compositions obtained from SERPENT simulations 
of one Savannah River Site reactors. Taken from Max Schalz' master 
thesis, see https://github.com/maxschalz/studious_potato/blob/main/data/SERPENT_outputs_NatU_percentages.npy
Simulations kindly provided by Antonio Figueroa.

BU:        0.5MWd       2MWd
  U230: 2.701e-21  2.701e-21
  U231: 1.422e-22  1.422e-22
  U232: 2.850e-14  2.850e-14
  U233: 3.099e-11  3.099e-11
  U234: 8.513e-05  8.513e-05
  U235: 1.030e-02  1.030e-02
  U236: 9.760e-05  9.760e-05
  U237: 8.535e-07  8.535e-07
  U238: 9.885e-01  9.885e-01
  U239: 5.391e-07  5.391e-07
  U240: 1.744e-10  1.744e-10
  U241: 1.100e-15  1.100e-15
  U242: 3.389e-20  3.389e-20
 Pu238: 9.167e-09  9.167e-09
 Pu239: 4.052e-04  4.052e-04
 Pu240: 9.193e-06  9.193e-06
 Pu241: 5.679e-07  5.679e-07
 Pu242: 5.122e-09  5.122e-09
 Pu243: 1.835e-12  1.835e-12
 Pu244: 4.362e-15  4.362e-15
 Np239: 7.775e-05  7.775e-05
 Np240: 1.160e-09  1.160e-09
 Np241: 3.059e-15  3.059e-15
 Np242: 4.438e-21  4.438e-21
 U235m: 8.643e-10  8.643e-10
Np240m: 2.402e-10  2.402e-10

"""

if __name__=="__main__":
    main()
