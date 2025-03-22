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

## Using this code: 

(The code is of course capable of much more than the following, but the following examples should get one started if one has no idea what to do with the code. After that, one can simply look around in the code to see what else there is to do.)

# To generate a single power spectrum

1. create a cosmology dictionary. A great way to start is to simply call camb_interface.default_cosmology(), which gives the Aletheia model 0 cosmology, and the user can modify it according to the desired parameter values.

2. Call camb_interface.evaluate_cosmology, passing in the cosmology dictionary. The user gets the power spectrum as well as various other information described in the evaluate_cosmology docstring. Alternatively, the user can call camb_interface.get_CAMB_interpolator, which is a spiritually analogous function that Andrea Pezzotta probably would have preferred. However, in the latter case the output is an interpolator rather than a power spectrum. This is a regular CAMB interpolator and should be accessed as one, i.e. interpolator.P(z, k).

