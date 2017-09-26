These are the C-generated files from neurodamus branch `coreneuronsetup`, commit `16859cd29b78a5d251a0ab2c23dbff43d423a773` from 10 Aug 2016.
Few minor changes
- Added initialization of shadow rhs and shadow didv on `nrn_init` methods
- added `nrn_current_parallel` methods
- StochKv only: added `#define usingR123 0` to disable Random123 generation, to fix crash
- For Point Processes only: added `net_receive2` methods as a wrapper of `net_receive`
- Changed `ProbAMPANMDA_EMS.mod` and `ProbGABAAB_EMS.mod` to allow for CVODE-based interpolation of `states`
