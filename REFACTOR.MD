

Some big TODOs for the refactor:

- switch to hdf5 (faster, cleaner, more compact)
- rather than directly suppling rates, the user should supply plasma variables, with rates then calculated by the code
(though directly supplying rates should still be possible)
- Don't always tally everything, allow the user to select which variables they want estimated.
- Change name to ORCA-MC to avoid name conflict.
