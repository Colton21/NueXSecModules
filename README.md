# NueXSecModules
Larsoft modules for use in the Nue Cross Section Analysis

## Running and Using the Module

The module is a simple method to create a structured data products in an intiutive manner and links important information like the reconstructed and true quantities, such that the end user does not need to run any reco-true matching code prior to analysis.

Be sure that you add NueXSecModules to your CMakeLists.txt before building with larsoft and then build with `mrb i`.

The total package is currently composed of three modules:

1. Reco-true matching module - this module constructed the "MC Ghosts".
2. Orders the recob::Track/Shower objects into TPC Objects and creates associations.
3. Analysis module to save TPC Objects and recob::Particles into an output ROOT file.

The initial idea for this is based on Marco's code, which constructs similar TPC Objects, which contain primarily associations to recob::Tracks/Showers. My module expands to save the dataproducts of the recob::Tracks/Showers as well as their associated MC Truth variables.

The output root file is composed of `std::vector<xsecAna::ParticleContainer>` which compose a `xsecAna::TPCObjectContainer` object. Based on the Pandora particle hierarchy, the event will have any number of `xsecAna::TPCObjectContainer`, which are saved in an `std::vector<xsecAna::TPCObjectContainer>`.

The `xsecAna` directory contains the modules which are used in LArsoft. In this directory you can find the fcl file among the class definitions. You can run the module just like any other larsoft module: `lar -c run_Xsec.fcl -s source_file_path`.

If anyone intends to apply these two container classes outside of my module's framework, then simply be aware of the following: as I have constructed two custom data products, we need to make sure that the proper library information is generated. This means looking at the Linkdef.h, classes_def.h, and classes.h files.

The scripts directory is where I have some test and analysis scripts. `out_inspect.cxx` is a simple ROOT script, which will print out the information contained in your generated file. Here you can check to make sure the modules are doing what you wish.

I will also be adding an analysis script which performs some cuts and generates some plots. If you find that you'd like to include all of the MC Truth information, then be sure to set the fcl parameter. There are a very large number of MCParticles and this inflates the size of the output file, as such they're only saved when requested. This should keep the size of the output file more managable and improve speed of the analysis scritps.

**Important**: as we are using a set of custom data products in the ROOT output file, we need to tell ROOT how to read these. If you try to run the scripts before doing this step, you will likely be told to generate a dictionary file or two. To avoid this problem first you'll want to make sure you `make clean` and then `make` in the `xsecAna` directory. This make file will compile the `TPCObjectContainer` and `ParticleContainer` classes and construct two libraries: `libxsecAna.so` and `libxsecAna_dict.so`. These libraries are what will tell ROOT how to understand the output file. If your `$LD_LIBRARY_PATH` is set to the directory where you run your scripts, you can either copy the libraries over, or append to your library path by doing:

```bash
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/library/
export LD_LIBRARY_PATH
./script.exe
```

This means, however, whenever you change the `TPCObjectContainer` or `ParticleContainer` classes you need to make sure you recompile and ensure the new libraries are visible to your script.

###Calculating the Cross Section

The cross section is calculated at the end of the scripts/selection.cxx code. Using whatever the final results of the selection are, the values are plugged into the standard formula to caculate the cross section. This is mostly striaght-forward, however calculating the flux can be tricky.

To properly calculate the flux we use three separate scripts. The first is found in potFinder, where you can run `potfinder.fcl` over the art ROOT files to produce `pot_extraction.root`, which contains the POT histogram for your files. This will take several minutes given the speed of access to your files. Then we can use  `scripts/potSum.cxx` to extract the total POT from the histogram.

Lastly, we need to consulte GENIE and the NuMI flux histograms to find the proper scaling value. This particular scaling factor will change based on your desired signal and beam. If you'd like to use it, I'll be happy to provide you with the scaling factor I use or the module and source files. Those on the gpvms can access `/uboone/app/users/chill2/flux_scripts`. The script is simple enough so I simply execute using root: `root -l flux_calc.cxx+`.


## Open Tickets

This is a place where I will be tracking some issues:

- A good amount of true variables are not being filled properly - some look like they are pointing to addresses, others are simply zero.
- Some MC Truth variables are simply not set - `xsecAna::TPCObjectContainer.MCHits()` 
- Certain truth variables being filled have odd results - why can MC Direction be greater than 1?
- Reco Momentum for recob::Track objects using StartMomentum() is always 1.
- Cosmic information using MCParticles is not great (vtx where Corsika generated, etc), need to use MCTrack/MCShower for these, but no associations in the MC files.

## To-Be Verified
- All NuMI Nue + Cosmic sample.
- All NuMI In-Time Cosmic sample - looking at the timing distribution, I think it might be wrong.
- All NuMI + Cosmic sample.


## Acknowledgements

Big thanks to everyone who helped with my getting this working and for ideas:

Gianluca, Marco, Corey, Elena.


