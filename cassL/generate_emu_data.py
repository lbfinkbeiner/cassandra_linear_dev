"""
This script is abbreviated in comments and the important statement as ged.

It provides an organized framework within which the user can make repeated calls
to the camb_interface script. The repeated calls are used to build up
collections of CAMB output spectra to use as training data sets for the
emulator.

Docstrings use the NumPy format: https://numpydoc.readthedocs.io/en/latest/format.html
"""

import sys, platform, os
import traceback
import copy as cp

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar

import camb
from camb import model, initialpower, get_matter_power_interpolator

from cassL import camb_interface as ci
from cassL import user_interface as ui
from cassL import utils

h_DEFAULT = ci.default_cosmology()["h"]

# The order of this list must match the order used in priors files, which are
# located in cassandra_linear_dev/priors/
COSMO_PARS_INDICES = [
    'ombh2',
    'omch2',
    'n_s',
    'sigma12',
    'A_s',
    'omnuh2',
    'h',
    'omkh2',
    'w0',
    'wa',
    'z'
]   

def labels_to_mapping(labels):
    """
    Return a mapping array (for use with denormalize_row and build_cosmology)
    based on the order of the parameters in the provided array.

    The point is to provide a set of indices which allow immediate accessing of
    the prior values. All of the prior files must abide by the order shown in
    COSMO_PARS_INDICES. 
    
    The priors indices always correspond to the same parameters, but the
    lhs_row may not contain all such parameters. That's why this mapping is 
    important. The keys (i.e. the indices, since this is a list) are are lhs row 
    indices while the values are prior indices.
    
    Parameters
    ----------
    labels: list
        labels in the format exemplified by COSMO_PARS_INDICES.
        
    Returns
    -------
    mapping: list
        This is what the rest of the script refers to as a "mapping": basically,
        we turn labels into a list of equal size (unless it contains duplicates)
        but where the cosmological variable names are replaced by the index
        associated with that variable name in the COSMO_PARS_INDICES list.
    """
    mapping = []
    for label in labels:
        if label not in COSMO_PARS_INDICES:
            raise ValueError("Unrecognized column: " + label)
        mapping.append(COSMO_PARS_INDICES.index(label))

    if len(mapping) != len(np.unique(mapping)):
        raise ValueError("Duplicate labels detected!")

    return mapping


def denormalize_row(lhs_row, priors, mapping):
    """
    Uses priors to convert the normalized lhs_row values into values that
    directly describe cosmological parameters.
    
    Parameters
    ----------
    lhs_row: list
        A list of unitless values between 0 and 1 taken from a Latin hypercube
        of input values for a data set of power spectra.
    priors: two-dimensional np.ndarray
        set of prior bounds, see ui.prior_file_to_array
    mapping: list.
        See the fn labels_to_mapping for more details concerning this object.
        
    Returns
    -------
    list
        A scaled and offset version of lhs_row which directly describes values
        of cosmological parameters.
    """
    tailored_priors = []
    for i in range(len(lhs_row)):
        tailored_priors.append(priors[mapping[i]])
    
    tailored_priors = np.array(tailored_priors)
    xrange = np.ptp(tailored_priors, axis=1)
    xmin = np.min(tailored_priors, axis=1)
    return lhs_row * xrange + xmin
    

def generate_header(priors, mapping, save_label):
    """
    Writes a header to a file specified by save_label. The header is intended to
    include all information necessary for the reproduction of the data set.
    
    The header shall explain the priors used and which columns map to which
    cosmological parameters. This essentially provides self-documenting data
    sets.
    
    However: although this fn is planned to be incorporated in an essential way
    into the `fill_hypercube*` series of functions, it hasn't been tested or
    fully implemented yet. (TO-DO)
    
    Parameters
    ----------
    priors: np.ndarray (two-dimensional)
        Set of prior bounds, see `ui.prior_file_to_array`.
    mapping: list (one-dimensional)
        See the fn `labels_to_mapping` for more details concerning this object.
    save_label: str
        File name associated with the set of power spectra to be 
        generated, for example using fill_hypercube_with_Pk.    
        
    Returns
    -------
    None
    """
    header = "HEADER FOR " + save_label + "\n------------\n\n"
    
    for i in range(len(mapping)):
        par_label = COSMO_PARS_INDICES[mapping[i]]
        header += "Column #" + str(i) + "of the LHS describes the parameter: " \
            + par_label + ".\n"
        
        these_priors = priors[mapping[i]]
        header += "its priors are [" +  these_priors[0] + ", " + \
            these_priors[1] + "].\n\n"
    
    with open(save_label + ".csv", "w") as file:
        file.write(header)
    
    # I am currently in the middle of developing the header system as a
    # replacement for the awful mapping system. Therefore, we print out the
    # header here to facilitate debugging.
    print("The requested header is as follows:\n\n", header)
    

def build_cosmology(lhs_row, mapping):
    """
    Generation of the training/testing data sets for the emulator works by
    iterating through rows of a Latin hypercube matrix and converting these
    values into values for cosmological parameters.
    
    This function creates a cosmology dictionary from the values in a row
    (lhs_row) taken from the Latin hypercube of input parameters. We do this
    because the cosmology dictionary attaches raw floating-point values to
    specific cosmological variables. The camb_interface script is written to
    call CAMB to obtain power spectra using these cosmology dictionaries.
    
    For a detailed explanation of the cosmology dictionary, see the text file
    "standards.txt" in the same directory as this script.
    
    Parameters
    ----------
    lhs_row: list (one-dimensional)
        A list of unitless values between 0 and 1 taken from a Latin hypercube
        of input values for a data set of power spectra.
    mapping: list (one-dimensional)
        See the fn `labels_to_mapping` for more details concerning this object.    
    
    Returns
    -------
    dict
        A cosmology dictionary where values of cosmological parameters are
        determined by the normalized values in lhs_row and mapped to specific
        cosmological parameters using mapping.
    """
    
    # Use Aletheia model 0 as a base
    cosmology = ci.default_cosmology(z_comparisons=False)
    
    for i in range(len(lhs_row)):
        # Two parameters require special handling: curvature and neutrinos
        par_label = COSMO_PARS_INDICES[mapping[i]]
        if par_label == 'omnuh2': # special because we have a helper fn here
            cosmology = ci.specify_neutrino_mass(cosmology, lhs_row[i], 1)
        elif par_label == 'omkh2': # special because camb doesn't like omkh2
            # h ought to have already been specified by now.
            cosmology['OmK'] = lhs_row[i] / cosmology['h'] ** 2
        else:
            cosmology[COSMO_PARS_INDICES[mapping[i]]] = lhs_row[i]
    
    return cosmology

def broadcast_unsolvable(input_cosmology, list_sigma12=None):
    """
    Helper fn for the fill_hypercube fns, called when a valid power spectrum
    cannot be obtained for a given cosmology. Warns the user and returns
    formatted Nones to be written to the output data set.
    
    The formatting of the Nones is important to ensuring that the fill_hypercube
    fns do not crash. Furthermore, by catching CAMB failures and noting them in
    the generated data sets, I've been able to automate removal of failed data 
    points later in the pipeline with the fn
    `train_emu.eliminate_unusable_entries`.
    
    Parameters
    ----------
    input_cosmology: dict
        This is only used to provide debugging output. This fn does not actually
        test whether input_cosmology indeed cannot be used to generate a power
        spectrum.
    list_sigma12:
        TO-DO: how do I articulate what this is?
    """
    print("\nThis cell cannot be solved with a nonnegative redshift.")
    print("This is the failed cosmology:\n")
    ui.print_cosmology(input_cosmology)
    
    # Something about this calculation seems bogus to me since the values given
    # are frequently absurd. Don't take this value so seriously. TO-DO
    if list_sigma12 is not None:
        print("\nThe closest feasible sigma12 value would yield an error of:",
            utils.percent_error(input_cosmology["sigma12"],
                list_sigma12[len(list_sigma12) - 1]), "%\n")

    return None, np.array([np.nan, np.nan, np.nan])

def direct_eval_cell(input_cosmology, standard_k_axis):
    """
    Helper fn for fill_hypercube_with_Pk: generates a single power spectrum
    from a single cosmology dictionary. Modifies h and z values so that the
    power spectrum exhibits a desired sigma12 value as specified in the
    cosmology dictionary.
    
    This script offers two approaches to evaluating individual cosmologies;
    they're intended to give the same results. Why did I write both? To aid in
    debugging; the two approaches differ in what parts of the CAMB code are 
    ultimately accessed.
    
    This fn first generates a power spectrum and then builds an interpolator
    from it.
        
    Parameters
    ----------
    input_cosmology: dict
        Defines a particular cosmology. See the docstring for
        `camb_interface.input_cosmology` for important details.
    standard_k_axis: list or one-dimensional np.ndarray
        The k values at which the power spectrum for `input_cosmology` should be
        evaluated. I don't believe there's a way to feed these directly into
        CAMB; one has to use interpolation.
    
    Returns
    -------
    p: np.ndarray
        Power spectrum in Mpc units computed by CAMB for the cosmology defined
        by input_cosmology. The length of this array is equal to that of
        standard_k_axis.
        
    np.ndarray
        A triplet of values useful for debugging:
        actual_sigma12: int
            The value of h and z are modified in order to get a power spectrum
            with a desired value of sigma12. If the desired value diverges
            significantly from `actual_sigma12`, something has gone wrong.
        float
            The value of h at which the cosmology given by input_cosmology has
            been evaluated to arrive at `p`.
        float
            The value of z at which the cosmology given by input_cosmology has
            been evaluated to arrive at `p`.
            
        The user is advised to write these triplets to a file in order to
        investigate any issues with the accuracy of the power spectra generated
        via this fn.
    """
    num_k_points = len(standard_k_axis)
    
    # redshifts at which to test the cosmology, to best match desired sigma12
    _redshifts=np.flip(np.geomspace(1, 11, 150) - 1)

    z_best = None
    p = None
    list_sigma12 = [None]
    
    while True:
        MEMNeC = ci.balance_neutrinos_with_CDM(input_cosmology, 0)

        # Attempt to solve for the power spectrum at all test redshifts. This
        # call to CAMB shouldn't throw an error on the first iteration of this
        # loop, but as you'll see a few lines later, we may begin to decrease
        # h in order to match the desired sigma12 value. This can sometimes
        # cause CAMB to break even if we satisfy CAMB's nominal requirement
        # that h >= 0.01.
        try:
            _, _, _, list_sigma12 = ci.evaluate_cosmology(MEMNeC, _redshifts,
                fancy_neutrinos=False, k_points=num_k_points,
                hubble_units=False)
        except Exception:
            return broadcast_unsolvable(input_cosmology, list_sigma12)

        interpolator = interp1d(np.flip(list_sigma12), np.flip(_redshifts),
            kind='cubic')

        try:
            # In the vast majority of cases, the following line will generate
            # the ValueError caught in the except clause. However, on rare
            # occasions a ValueError is generated when obtaining the
            # interpolated P(k). Regardless of which line raises the Exception,
            # the solution appears to be the same: just decrease h.
            z_best = interpolator(input_cosmology["sigma12"])
            
            k, _, p, actual_sigma12 = ci.evaluate_cosmology(input_cosmology,
                redshifts=np.array([z_best]), fancy_neutrinos=False,
                k_points=num_k_points) 
                
            if input_cosmology['h'] != h_DEFAULT: # we've touched h,
                # we need to interpolate
                interpolator = interp1d(k, p, kind="cubic")
                p = interpolator(standard_k_axis)
            
            break
            
        except ValueError:
            # The try block almost certainly failed because there is no
            # solution. That is to say, varying redshift alone, and avoiding
            # all negative values of redshift, the desired value of sigma12
            # cannot be reached. Therefore, we need to start playing with h.
            if input_cosmology['h'] > 0.1:
                input_cosmology['h'] -= 0.1
            elif input_cosmology['h'] >= 0.02:
                # Finer-grained decreases will save the most extreme
                # cosmologies in our priors
                input_cosmology['h'] -= 0.01
            else: # We can't decrease h any further.
                return broadcast_unsolvable(input_cosmology, list_sigma12)
    
    if len(p) == 1: # de-nest the power spectrum
        p = p[0]
    
    # If neutrinos are not massless, we have to run again in order to
    # get the correct sigma12 value...
    if input_cosmology['omnuh2'] != 0:
        _, _, _, actual_sigma12 = ci.evaluate_cosmology(MEMNeC,
            redshifts=np.array([z_best]), fancy_neutrinos=False,
            k_points=num_k_points)
    # De-nest
    actual_sigma12 = actual_sigma12[0]

    if input_cosmology['h'] != h_DEFAULT: # announce that we've touched h
        print("We had to move h to", np.around(input_cosmology['h'], 3))

    # We don't need to return k because we take for granted that all
    # runs will have the same k axis.

    return p, np.array((actual_sigma12, input_cosmology['h'], float(z_best)))


def interpolate_cell(input_cosmology, standard_k_axis):
    """
    Helper fn for `fill_hypercube_with_Pk`: generates a single power spectrum
    from a single cosmology dictionary. Modifies h and z values so that the
    power spectrum exhibits a desired sigma12 value as specified in the
    cosmology dictionary.
    
    This script offers two approaches to evaluating individual cosmologies;
    they're intended to give the same results. Why did I write both? To aid in
    debugging; the two approaches differ in what parts of the CAMB code are 
    ultimately accessed.
    
    This fn directly uses a CAMB interpolator. 
    
    TO-DO: this fn is similar to `direct_eval_cell` in many regards, and should
    be combined with it; a fn parameter should allow the user to choose between
    the two functionalities currently afforded by their separation.
                
    Parameters
    ----------
    input_cosmology: dict
        Defines a particular cosmology. See the docstring for
        `camb_interface.input_cosmology` for important details.
    standard_k_axis: list or one-dimensional np.ndarray
        The k values at which the power spectrum should be evaluated. I don't
        believe there's a way to feed these directly into CAMB; one has to use
        interpolation.
    
    Returns
    -------
    p: np.ndarray
        Power spectrum computed by CAMB for the cosmology defined by
        input_cosmology. The length of this array is equal to that of
        standard_k_axis.
        
    np.ndarray
        A triplet of values useful for debugging:
        actual_sigma12: int
            The value of h and z are modified in order to get a power spectrum
            with a desired value of sigma12. If the desired value diverges
            significantly from `actual_sigma12`, something has gone wrong.
        float
            The value of h at which the cosmology given by input_cosmology has
            been evaluated to arrive at `p`.
        float
            The value of z at which the cosmology given by input_cosmology has
            been evaluated to arrive at `p`.
            
        The user is advised to write these triplets to a file in order to
        investigate any issues with the accuracy of the power spectra generated
        via this fn.
    """
    # This allows us to roughly find the z corresponding to the sigma12 that we
    # want.
    _redshifts=np.flip(np.geomspace(1, 11, 150) - 1)
    z_best = None
    p = None
    list_sigma12 = [None]
    #! Hard code
    k_max = 1.01 * max(standard_k_axis)

    while True:
        MEMNeC = ci.balance_neutrinos_with_CDM(input_cosmology, 0)
       
        try:
            MEMNeC_p_interpolator = ci.cosmology_to_Pk_interpolator(MEMNeC,
                    redshifts=_redshifts, kmax=k_max, hubble_units=False)

            s12intrp = ci.sigma12_from_interpolator
            sigma12 = lambda z: s12intrp(MEMNeC_p_interpolator, z)
            list_sigma12 = np.array([sigma12(z) for z in _redshifts])
        except Exception:
            return broadcast_unsolvable(input_cosmology, list_sigma12)

        interpolator = interp1d(np.flip(list_sigma12), np.flip(_redshifts),
            kind='cubic')

        try:
            z_best = interpolator(input_cosmology["sigma12"])
            interpolation_redshifts = np.flip(np.linspace(max(0, z_best - 1),
                                                          z_best + 1, 150))

            get_intrp = ci.cosmology_to_Pk_interpolator
            p_intrp = get_intrp(input_cosmology,
                                redshifts=interpolation_redshifts,
                                kmax=k_max, hubble_units=False)
            p = p_intrp.P(z_best, standard_k_axis)

            break

        except ValueError:
            # we need to start playing with h.
            if input_cosmology['h'] > 0.1:
                input_cosmology['h'] -= 0.1
            elif input_cosmology['h'] >= 0.02:
                # Finer-grained decreases might save a couple of weird cosmologies
                input_cosmology['h'] -= 0.01
            else: # We can't decrease any further.
                return broadcast_unsolvable(input_cosmology, list_sigma12)

    actual_sigma12 = ci.sigma12_from_interpolator(MEMNeC_p_interpolator,
                                                  z_best)

    if input_cosmology['h'] != h_DEFAULT: # announce that we've touched h
        print("We had to move h to", np.around(input_cosmology['h'], 3))

    # We don't need to return k because we take for granted that all
    # runs will have the same k axis.

    return p, np.array((actual_sigma12, input_cosmology['h'], float(z_best)))


def fill_hypercube_with_sigma12(lhs, mapping, priors, samples=None,
                                write_period=None, cell_range=None,
                                save_label="sigma12",
                                crash_when_unsolvable=False):
    """
    Given a Latin hypercube of cosmologies, return an array of sigma12 values.
    This fn is used to build the sigma12 emulator.
    
    This fn is probably obsolete in light of the newer, generalized fn 
    `fill_hypercube_with_sigmaR`.
        
    Parameters
    ----------
    lhs: np.ndarray (of uniform two-dimensional shape)
        Normalized values corresponding to cosmological parameters used to
        define cosmologies to be evaluated. The width of this matrix corresponds
        to the number of cosmological parameters used to define each cosmology.
        The height of this matrix corresponds to the number of different
        cosmologies defined.
    mapping: list (one-dimensional)
        Used to build cosmologies from `lhs`. See the fn `labels_to_mapping` for
        more details concerning `mapping` objects.
    priors: np.ndarray (two-dimensional)
        Collection of bounds on cosmological priors. Used to build cosmologies
        from `lhs`. For more details, see `ui.prior_file_to_array`.
    save_label: str
        Prefix used to name files written to disk by this fn as part of saving
        its progress through the batch of cosmologies defined by `lhs`,
        `mapping`, and `priors`. In other words, this is the "name" of the run
        as a whole. 
    
        TO-DO: according to the NumPy docstring example, this parameter should
        go in the section "Other Parameters" since it comes after at least one
        parameter (in this case, for example, `samples`) explained in the
        section "Other Parameters". I'm explaining it here since it's an
        important parameter, but I should fix the order and then correct any
        references to this fn.
        
    Returns
    -------
    samples: np.ndarray (one-dimensional)
        Contains the calculated sigma12 values for each cosmology jointly
        specified by `lhs`, `mapping`, and `priors`. The length of `samples`
        therefore corresponds to the height of `lhs`.
        
    Other Parameters
    ----------------
    samples: np.ndarray, optional
        TO-DO: explain what this controls
        Collection of already-computed results. Useful
        When specifying this parameter, `cell_range` should also be set.
        The default value is None, which indicates to this fn that the hypercube
        should be generated from scratch.
        
        NB: This parameter makes sense for fill_hypercube_with_Pk, which is a
        much more computationally intensive fn. It is of comparatively little
        use here, since a complete batch of results can be generated in but few
        minutes.
    write_period: int, optional
        TO-DO: explain what this controls
    cell_range:
        TO-DO: explain what this controls
    
        See the NB 
    crash_when_unsolvable:
        TO-DO: explain what this controls
    """    
    def eval_func(cosmology):
        # De-nesting
        return ci.evaluate_sigma12(cosmology, [12.], [1.])[0][0]

    if cell_range is None:
        cell_range = range(len(lhs))
    if samples is None:            
        samples = np.zeros(len(lhs))

    unwritten_cells = 0
    for i in cell_range:
        this_denormalized_row = denormalize_row(lhs[i], priors, mapping)
        this_cosmology = build_cosmology(this_denormalized_row, mapping)

        try:
            samples[i] = eval_func(this_cosmology)
        except camb.CAMBError:
            print("This cell is unsolvable. Since this function requires " + \
                  "no rescaling, your priors are probably extreme.")
            if crash_when_unsolvable:
                raise ValueError("Cell unsolvable.")

        print(i, "complete")
        unwritten_cells += 1

        # Here we use much simpler logic compared to the analogous section of
        # the other two fill_hypercube functions, because the user should very
        # rarely be using cell_range, and writing the whole set of results
        # should not be burdensome.
        if write_period is not None and unwritten_cells >= write_period:
            file_suffix = "_backup_i{}"
            np.save(save_label + file_suffix.format(i), samples)
            unwritten_cells = 0

    return samples
    
    
def fill_hypercube_with_sigmaR(lhs, R_axis, mapping, priors, cell_range=None,
                               write_period=None, save_label="sigmaR",
                               crash_when_unsolvable=False):
    """
    @lhs: this is a list of tuples with which @eval_func is to be evaluated.

    @cell_range: a range object specifying the indices of lhs which still need
        to be evaluated. By default, it is None, which means that the entire
        lhs will be evaluated. This parameter can be used to pick up from where
        previous runs left off, and to run this method in saveable chunks.
        
    Parameters
    ----------
    
    Returns
    -------
    """
    if cell_range is None:
        cell_range = range(len(lhs))

    samples = np.zeros((len(lhs), len(R_axis)))

    unwritten_cells = 0
    for i in cell_range:
        this_denormalized_row = denormalize_row(lhs[i], priors, mapping)
        this_cosmology = build_cosmology(this_denormalized_row, mapping)

        try:
            samples[i] = ci.evaluate_sigmaR(this_cosmology, R_axis, [1.])[0]

            if samples[i] is None and crash_when_unsolvable:
                raise ValueError("Cell unsolvable.")
        except camb.CAMBError:
            print("This cell is unsolvable. Since this function requires " + \
                  "no rescaling, your priors are probably extreme.")
            if crash_when_unsolvable:
                raise ValueError("Cell unsolvable.")

        print(i, "complete")
        unwritten_cells += 1
        if write_period is not None and unwritten_cells >= write_period:
            # We add one because the current cell is also unwritten
            save_start = i - unwritten_cells + 1
            save_end = i + 1

            file_suffix = "_backup_i{}_through_{}_{}.npy"
            file_suffix = file_suffix.format(save_start, i, save_label)

            np.save("sigmaR" + file_suffix, samples[save_start:save_end])

            unwritten_cells = 0

    return samples    


def fill_hypercube_with_Pk(lhs, standard_k_axis, mapping, priors,
                           eval_func=direct_eval_cell, cell_range=None,
                           write_period=None, save_label="Pk",
                           crash_when_unsolvable=False):
    """
    @lhs: this is a list of tuples with which @eval_func is to be evaluated.

    @cell_range: a range object specifying the indices of lhs which still need
        to be evaluated. By default, it is None, which means that the entire
        lhs will be evaluated. This parameter can be used to pick up from where
        previous runs left off, and to run this method in saveable chunks.
        
    Parameters
    ----------
    
    Returns
    -------
    """
    if cell_range is None:
        cell_range = range(len(lhs))

    samples = np.zeros((len(lhs), len(standard_k_axis)))

    # The rescaling parameters are true_sigma12, h and z. Only the sigma12
    # value is used to train the emu but the rest provide debugging information
    # and sure that the output spectra are easily reproducible.
    rescalers_arr = np.zeros((len(lhs), 3))

    unwritten_cells = 0
    for i in cell_range:
        this_p = None
        this_denormalized_row = denormalize_row(lhs[i], priors, mapping)
        this_cosmology = build_cosmology(this_denormalized_row, mapping)

        try:
            samples[i], rescalers_arr[i] = \
                eval_func(this_cosmology, standard_k_axis)

            if samples[i] is None and crash_when_unsolvable:
                raise ValueError("Cell unsolvable.")
        except camb.CAMBError:
            print("This cell is unsolvable. However, in this case, we " + \
                  "observed a CAMBError rather than a negative redshift. " + \
                  "This suggests that there is a problem with the input " + \
                  "hypercube.")
            if crash_when_unsolvable:
                raise ValueError("Cell unsolvable.")

        actual_sigma12 = rescalers_arr[i][0]
        if not np.isnan(actual_sigma12): # we need to re-normalize
            prior = priors[3]
            normalized = (actual_sigma12 - prior[0]) / (prior[1] - prior[0])
            rescalers_arr[i][0] = normalized

        print(i, "complete")
        unwritten_cells += 1
        if write_period is not None and unwritten_cells >= write_period:
            # We add one because the current cell is also unwritten
            save_start = i - unwritten_cells + 1
            save_end = i + 1

            file_suffix = "_backup_i{}_through_{}_{}.npy"
            file_suffix = file_suffix.format(save_start, i, save_label)

            np.save("Pk" + file_suffix, samples[save_start:save_end])
            np.save("rescalers" + file_suffix,
                    rescalers_arr[save_start:save_end])

            unwritten_cells = 0

    return samples, rescalers_arr


def interpolate_nosigma12(input_cosmology, standard_k_axis):
    """
    Returns the power spectrum in Mpc units and the actual sigma12_tilde value
    to which it corresponds. The `input_cosmology` dictionary must specify a
    redshift value, of course.

    TO-DO: as of 03.04.2025, I don't remember what I was doing with this fn
    (it's not called by any other part of at least this script). My suspicion is
    that I meant to use this fn to offer a version of
    `fill_hypercube_with_sigma12` using CAMB interpolators rather than the
    direct-evaluation approach, since there are (negligible) differences in
    the answers given by the two approaches.
    
    In case I eventually have time to investigate and check my hypothesis, I
    will leave this remnant from the previous iteration of this docstring:
    "This is a demo function until we figure out how to apply the interpolation
    approach to the generation of emu data. Once we have that, we can figure
    out how to re-combine this function with the previous one." Before I moved
    this fn, "the previous one" referred to `interpolate_cell`, so perhaps my
    hypothesis doesn't really make sense...
        
    Parameters
    ----------
    input_cosmology: dict
        Defines a particular cosmology. See the docstring for
        `camb_interface.input_cosmology` for important details.
    standard_k_axis: list or one-dimensional np.ndarray
        The k values at which the power spectrum for `input_cosmology` should be
        evaluated.
    
    Returns
    -------
    tuple:
        the first return is an np.ndarray representing the power spectrum
        evaluated at each point given in `standard_k_axis`
        
        the second return is an np.ndarray where the first index corresponds to
        the sigma12 value when evaluating `input_cosmology` (which should
        contain a prescribed redshift value) and the remaining indices are
        filled in with np.nan. The purpose of this arrangement is to conform
        the output of this fn to conventional formats used by other fns in this
        script.
    """
    # This allows us to roughly find the z corresponding to the sigma12 that we
    # want.
    #! Hard code
    k_max = 1.01 * np.max(standard_k_axis)
    z = input_cosmology["z"]

    get_intrp = ci.cosmology_to_Pk_interpolator
    
    interpolation_redshifts = np.flip(np.linspace(max(0, z - 1), z + 1, 150))
    p_intrp = get_intrp(input_cosmology, redshifts=interpolation_redshifts,
                        kmax=k_max, hubble_units=False)
    p = p_intrp.P(z, standard_k_axis)

    actual_sigma12 = ci.sigma12_from_interpolator(p_intrp, z)

    # We don't need to return k because we take for granted that all
    # runs will have the same k axis.
    return p, np.array([actual_sigma12, np.nan, np.nan])
