# NueXSecModules
Larsoft modules for use in the Nue Cross Section Analysis

## Running and Using the Module

The module is a simple method to create a structured data products in an intiutive manner and links important information like the reconstructed and true quantities, such that the end user does not need to run any reco-true matching code prior to analysis.



The total package is currently composed of three modules:

1. Reco-true matching module - this module constructed the "MC Ghosts".
2. Orders the recob::Track/Shower objects into TPC Objects and creates associations.
3. Analysis module to save TPC Objects and recob::Particles into an output ROOT file.

The initial idea for this is based on Marco's code, which constructs similar TPC Objects, which contain primarily associations to recob::Tracks/Showers. My module expands to save the dataproducts of the recob::Tracks/Showers as well as their associated MC Truth variables.

The output root file is composed of `std::vector<xsecAna::ParticleContainer>` which compose a `xsecAna::TPCObjectContainer` object. Based on the Pandora particle hierarchy, the event will have any number of `xsecAna::TPCObjectContainer`, which are saved in an `std::vector<xsecAna::TPCObjectContainer>`.

The `xsecAna` directory contains the modules which are used in LArsoft. In this directory you can find the fcl file among the class definitions. You can run the module just like any other larsoft module: `lar -c run_Xsec.fcl -s source_file_path`.

If anyone intends to apply these two container classes outside of my module's framework, then simply be aware of the following: as I have constructed two custom data products, we need to make sure that the proper library information is generated. This means looking at the Linkdef.h, classes_def.h, and classes.h files.

The scripts directory is where I have some test and analysis scripts. `out_inspect.cxx` is a simple ROOT script, which will print out the information contained in your generated file. Here you can check to make sure the modules are doing what you wish.

I will also be adding an analysis script which performs some cuts and generates some plots.

**Important**: as we are using a set of custom data products in the ROOT output file, we need to tell ROOT how to read these. If you try to run the scripts before doing this step, you will likely be told to generate a dictionary file or two. To avoid this problem first you'll want to make sure you `make clean` and then `make` in the `xsecAna` directory. This make file will compile the `TPCObjectContainer` and `ParticleContainer` classes and construct two libraries: `libxsecAna.so` and `libxsecAna_dict.so`. These libraries are what will tell ROOT how to understand the output file. If your `$LD_LIBRARY_PATH` is set to the directory where you run your scripts, you can either copy the libraries over, or append to your library path by doing:

```bash
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/library/
export LD_LIBRARY_PATH
./script.exe
```

This means, however, whenever you change the `TPCObjectContainer` or `ParticleContainer` classes you need to make sure you recompile and ensure the new libraries are visible to your script.



## Open Tickets

This is a place where I will be tracking some issues:

- Adding all MC Truth information in a new tree - this is useful to get true event numbers for efficiencies.
- Add optical information (`simpleFlashBeam`) - this is needed on a per-event basis to perform the flash-matching as well as an in-time cut.
- A good amount of true variables are not being filled properly - some look like they are pointing to addresses, others are simply zero. Does this come from the MC Ghost module? I notice not many `kBeamNeutrino` events have matched MC Ghosts, but more `kCosmicRay` events do.
- Some MC Truth variables are simply not set - `xsecAna::TPCObjectContainer.MCHits()` 

## Acknowledgements

Big thanks to everyone who helped with my getting this working and for ideas:

Gean-Luca, Marco, Corey, Elena.


