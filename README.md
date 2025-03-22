# cassandra_linear_dev
 Development tools for the emulation suite Cassandra Linear. To get started, install this repo in the typical Pythonic fashion (e.g. navigate to directory and run `pip install -e .`).

See also the repository https://github.com/lbfinkbeiner/Cassandra.linear
That repository is used for the user-facing code, which is supposed to be all that anyone will ultimately need in order to work with emulator files (.cle).

This repo is used to generate the .cle files in the first place. It is a significantly more complicated code.

Currently, the repo is not in a good state for three reasons:
1. I never finished documenting it / removing legacy code, so some parts of the code are basically unreadable and not all of the code is still relevant to the current functioning.

2. (This is the important one) there is almost certainly at least one serious bug in the code. Although the individual stages in the creation of an emulator were tested and approved, the full pipeline is not working correctly. At this point, even the individual stages cannot necessarily be trusted because I cannot remember how much I have changed in my attempts at fixing the problem. So whatever happens, we need to run extensive tests and investigate until the full pipeline can be trusted.

3. As an example of a big change mentioned in point 2: when I had to stop working on the code and switch to patent studies, I was still in the middle of switching the way that the emulator builder reads in parameter values in order to generate training power spectra with CAMB. The point of this exercise is to greatly simplify debugging and future emulator creation, because a developer will be able to look directly at headers in the training data to see what cosmology they describe.

As of 22.03.2025, I think I may yet have enough time to address points 1 and 3, but 2 could be a significant undertaking depending on the nature of the problem.

## Using this code: (TO-DO: I really should make new Jupyter notebooks to make it easier to follow these instructions...)

(The code is of course capable of much more than the following, but the following examples should get one started if one has no idea what to do with the code. After that, one can simply look around in the code to see what else there is to do.)

# To generate a single power spectrum

1. Create a cosmology dictionary. A great way to start is to simply call camb_interface.default_cosmology(), which gives the Aletheia model 0 cosmology, and the user can modify it according to the desired parameter values.

2. Call camb_interface.evaluate_cosmology, passing in the cosmology dictionary. The user gets the power spectrum as well as various other information described in the evaluate_cosmology docstring. Alternatively, the user can call camb_interface.get_CAMB_interpolator, which is a spiritually analogous function that Andrea Pezzotta probably would have preferred. However, in the latter case the output is an interpolator rather than a power spectrum. This is a regular CAMB interpolator and should be accessed as one, i.e. interpolator.P(z, k).

# To generate a data set

1. Create a Latin hypercube sample with lhc.multithread_unit_LHC_builder. It's multithreaded because I wrote the fn back when maximizing the minimum separation between Latin hypercube samples appeared important to the performance of the emulator, so I sought to harness as much compute as possible to generate a high-quality hypercube (the approach is to just randomly generate Latin hypercube samples and keep the one with the largest minimum separation). The results in my thesis paper seemed to indicate that the impact is negligible. Anyway, the user should generate a Latin hypercube sample with the desired sample size (# rows of the sample) and with the desired number of cosmological parameters (in massless-neutrino cosmologies, 4, otherwise 6).

Since the fn is designed to run indefinitely, one cannot use it to directly load a Latin hypercube sample into main memory. The fn will automatically update a file on disk with the results of the best run yet. I recommend that the user simply kill the shell after letting it run for a minute at most, since as I've said, the minimum separation seems to have little effect. Then, load the file into memory.
-> TO-DO: improve the documentation on this fn to at least explain all of the fn parameters.

2. Choose a set of priors. The sets are shown in the dirrectory cassandra_linear_dev/priors/. "FP" is short for "full pipeline" and at this point applies to all of the prior files despite some not including this string. It just means that it also includes evolution parameters so that one can test the full Cassandra-linear pipeline. "2PERCENT" refers to a debugging idea in which we shrink each parameter range by 2% in order to see if emulator stability can be improved by throwing out extreme cosmologies. I don't remember if this shrinking has helped in any way with the current crisis in the code.

3. The priors are to be read into main memory with the fn user_interface.prior_file_to_array.

4. Call generate_emu_data.fill_hypercube_with_Pk with the Latin hypercube sample array and the priors. The user will also need to specify a list of k values which shall be common to all output spectra; by default, CAMB will output different k-axis arrays for different power spectra (I think it has to do with the input cosmology), so to ensure comparability of the different spectra in the training set, we need to apply an interpolation layer. The interpolation layer ensures that for each column in the matrix of output spectra, a single common value of k is implied.

This is the slowest step and the various auto-saving parameters of the fn should be configured so that no progress is lost if the fn is interrupted while generating all of the output spectra.

# To make an emulator
-> TO-DO: Lukas has to fill in this section still...
