
import json
import numpy as np
import os

from . import kernel


#TODO
# - Find a workaround for the ugly calculations that are currently
#   performed in 'get_input_params' for the 'power_output' parameter
#   and that are only relevant (?) for the SRS reactor.

def main():
    return

def predict():
    """Calculate the spent fuel composition."""
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

    # Check if the needed kernels and parameter information exist.
    data_dir = os.path.join(os.path.split(__file__)[0], "..", "data")
    kernel_dir = os.path.join(data_dir, "trained_kernels")
    if not os.path.isdir(kernel_dir):
        raise OSError("'trained_kernels' directory not found!")
    for iso in isotopes:
        fname = f"{iso}.npy"
        if not os.path.isfile(os.path.join(kernel_dir, fname)):
            msg = f"Trained kernel '{fname}' not found!"
            raise FileNotFoundError(msg)

        fname = f"training_params_{iso}.json"
        if not os.path.isfile(os.path.join(kernel_dir, fname)):
            msg = f"Training parameters '{fname}' not found!"
            raise FileNotFoundError(msg)

    # Load input parameters, training data and the kernel type.
    training_data = np.load(os.path.join(data_dir, "x_trainingset.npy"),
                            allow_pickle=True)
    par = ("enrichment", "temperature", "power_output", "burnup")
    reactor_input_params = get_input_params(par)
    reactor_input_params = np.expand_dims(reactor_input_params, axis=0)

    # Calculate the spent fuel composition for all isotopes.
    spent_fuel_composition = {"spent_fuel_composition": {}}
    for iso in isotopes:
        kernel_fname = os.path.join(kernel_dir, f"{iso}.npy")
        params_fname = os.path.join(kernel_dir,
                                    f"training_params_{iso}.json")

        with open(params_fname, "r") as f:
            data = json.load(f)
            kernel_type = data["kernel_type"]
            size = data["size"]
        check_input_params(reactor_input_params, training_data[:size])
        trained_kernel = np.load(kernel_fname,
                                 allow_pickle=True).item()
        mass = run_kernel(reactor_input_params, training_data[:size],
                          trained_kernel, kernel_type)
        spent_fuel_composition["spent_fuel_composition"][iso] = mass

    store_results(spent_fuel_composition)

def check_input_params(params, training_data):
    """Check the validity of the input parameters.

    The function raises a ValueError if the checks fail, else it has no
    return value.
    """
    min_vals = np.min(training_data, axis=0)
    max_vals = np.max(training_data, axis=0)
    is_valid = (np.all(min_vals < params[0])
                and np.all(params[0] < max_vals))

    if not is_valid:
        msg = ("[spentfuelgpr] One or more parameters exceed the bounds.\n"
               + f"Minimum parameter values: {min_vals}\n"
               + f"Actual parameter values:  {params[0]}\n"
               + f"Maximum parameter values: {max_vals}")
        raise ValueError(msg)

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
    input_params : array
        An array containing the extracted values in the same order as
        specified in `pnames`.
    """
    fname = "gpr_reactor_input_params.json"
    with open(fname, "r") as f:
        data = json.load(f)
        input_params = []
        
        for param in pnames:
            param = param.lower()
            if param == "enrichment":
                enrich = data["fresh_fuel_composition"]["922350000"]
                input_params.append(enrich)
            elif param == "power_output":
                # These calculations below have to be performed for the
                # Savannah River Site reactor. Currently, they are
                # hardcoded but this is hopefully subject to change.
                # TODO update this implementation
                n_assemblies_tot = 500;
                n_assemblies_model = 18;
                feet_to_cm = 30.48;
                assembly_length = 12;  # in feet
                power = (data[param] * n_assemblies_model
                         / n_assemblies_tot / assembly_length
                         / feet_to_cm)
                power *= 1e6  # conversion MW to W
                input_params.append(power)
            else:
                input_params.append(data[param])
    
    return np.array(input_params)

def store_results(composition):
    """Store the composition in a .json file.

    Parameters
    ----------
    composition : dict
        A dictionary with the keys being the isotopes and the values
        being the corresponding masses.
    """
    fname = "gpr_reactor_spent_fuel_composition.json"
    with open(fname, "w") as f:
        json.dump(composition, f, indent=2)

def run_kernel(reactor_input_params, x_train, trained_kernel, 
               kernel_params, kernel_type='ASQE'):
    """Calculate the mass of one isotope in the spent fuel.

    Parameters
    ----------
    reactor_input_params
        Input parameters used in the calculations and as defined above
        in 'predict()' and in 'get_input_params'
    x_train
        Input parameters used during training
    trained_kernel
        Trained kernel, specific to the isotope in question
    kernel_type
        The type of the trained kernel

    Returns
    -------
    mu_s
        The (mean predicted) mass of the isotope in the spent fuel
        after one irradiation period.
    """
    if (kernel_type != 'ASQE'):
        msg = "Currently, only the 'ASQE' kernel type is supported."
        raise ValueError(msg)

    kernel_params = trained_kernel["Params"]
    alpha = trained_kernel["alpha_"]
    lambda_ = trained_kernel["LAMBDA"]
    k_s = kernel.Kernel(reactor_input_params, x_train, kernel_type,
                 kernel_params, gradient=False, LAMBDA=lambda_)
    mu_s = np.dot(k_s.T, alpha)
    
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
