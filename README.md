# NueXSecModules
Larsoft modules for use in the Nue Cross Section Analysis

## Running and Using the Module

The module is a simple method to create a structured data products in an intiutive manner and links important information like the reconstructed and true quantities, such that the end user does not need to run any reco-true matching code prior to analysis.

Be sure that you add NueXSecModules to your CMakeLists.txt before building with larsoft and then build with `mrb i`.

*For First Time Use assuming  MRB is available*

You can simply do: `source ana_setup_script.sh`. This cleans and builds both the larsoft repository (where the libraries live) and the analysis scripts. This **does not** build your release of larsoft! - You will still have to do an `mrb i` or `make install` as normal in your larsoft build dir if you want to work with any larsoft products/modules.

If you plan to do development of the analysis scripts (stand-alone C++ code), then simply build either of the `xsecAna` and `scripts` repositories by hand or run `source update_libraries.sh`. This ensures that the libraries are up-to-date, copied to all analysis directories, and all the c++ code is freshly built.

If you ever change any of the classes in `/xsecAna/` be sure to `source update_libraries.sh` again to build and update the libraries in `/scripts` and `/python`.
With the last update to my Makefile, if you are on the GPVM or on a device with SIP (system integrity protection) disabled, then ROOT should be able to find the libraries located in `xsecAna`. If this is not the case for you, simply copy the relevant libraries and ROOTMAP to the `scripts` repository. You'll need to do this any time you make a change to the library files.

##Larsoft Moduels
The total package is currently composed of three modules:

1. Reco-true matching module - this module constructed the "MC Ghosts".
2. Orders the recob::Track/Shower objects into TPC Objects and creates associations.
3. Analysis module to save TPC Objects and recob::Particles into an output ROOT file.

The initial idea for this is based on Marco's code, which constructs similar TPC Objects, which contain primarily associations to recob::Tracks/Showers. My module expands to save the dataproducts of the recob::Tracks/Showers as well as their associated MC Truth variables.

The output root file is composed of `std::vector<xsecAna::ParticleContainer>` which compose a `xsecAna::TPCObjectContainer` object. Based on the Pandora particle hierarchy, the event will have any number of `xsecAna::TPCObjectContainer`, which are saved in an `std::vector<xsecAna::TPCObjectContainer>`.

The `xsecAna` directory contains the modules which are used in LArsoft. In this directory you can find the fcl file among the class definitions. You can run the module just like any other larsoft module: `lar -c run_Xsec.fcl -s source_file_path`.

If anyone intends to apply these two container classes outside of my module's framework, then simply be aware of the following: as I have constructed two custom data products, we need to make sure that the proper library information is generated. This means looking at the Linkdef.h, classes_def.h, and classes.h files.

##Analysis and Test Scripts

Be sure to read the comments in the relevant files you're using!

The scripts directory is where I have some test and analysis scripts. `out_inspect.cxx` is a simple ROOT script, which will print out the information contained in your generated file. Here you can check to make sure the modules are doing what you wish.
As for performing the $\nu_{e}$ selection, you can run one of the `selection` scripts. The `selection_slim` script is simplest and a good place to understand the individual cuts implemented. Other versions of the `selection` scripts involve a good deal of plotting code. If you want to create ROOT plots based on the cut variables you should run one of these.
There are more details included in `main.h` which handles all of the various classes.
While I've tried to control for input handling `main.exe` expects any of the data/mc files to be listed in a particular order: MC, EXT, Data.
The selection is able to run without any of the listed files, however will *not* work as intended with (for example) only EXT and Data or MC and Data. But MC and EXT only works as intended.

```
./main.exe path_to_mc_file path_to_ext_file path_to_on_beam_data_file
```

Additional functionality includes providing the path to a custom configuration text file using `-c config.txt` (see utility.h for the formatting of the `configure` function used to set the cut values). If this is not set, the code will use a set of default parameters in `main.h`.
Also, the default running condition is to use the full selection and produce many many plots - use `--slim` to run more quickly and simply see the cut's performances.

##Python

The python directory is a new addition. These scripts can be used to more easily interface with the c++ code, and avoiding needing to recompile when changing cut values.
Look into the code to see which versions are being run (i.e. slim or not), but by default the `run_selection.py` acts as the simplest way of running the selection code. It creates a formatted config file based on the parameters in the python code and the only required inputs are as you would run the selection from the `.exe` but without worrying about the options.

```python
python run_selection.py mc_file ext_file data_file
```

Notice how the path of the `config.txt` is not included here - the script handles this for you.

If you're looking to vary cut parameters quickly, `iterate_cuts.py` is a skeleton (though working) version of various scripts which simply run a for loop over the c++ selection code. This can be somewhat time-consuming, so beware. But it also contains some of the matplotlib functionality for making efficiency and purity plots following the results of the selection code. (slim mode is default here to increase speed).

##Further Notes

If you find that you'd like to include all of the MC Truth information, then be sure to set the fcl parameter. There are a very large number of MCParticles and this inflates the size of the output file, as such they're only saved when requested. This should keep the size of the output file more managable and improve speed of the analysis scritps.

**Important**: If your analysis scripts are unable to find the library files read on: as we are using a set of custom data products in the ROOT output file, we need to tell ROOT how to read these. If you try to run the scripts before doing this step, you will likely be told to generate a dictionary file or two. To avoid this problem first you'll want to make sure you `make clean` and then `make` in the `xsecAna` directory. This make file will compile the `TPCObjectContainer` and `ParticleContainer` classes and construct a number of libraries. Simply make sure these are copied over to the `scripts` directory if you are unable to direct ROOT to the `xsecAna` directory (due to SIP for example).

Setting your $LD_LIBRARY_PATH may also solve this problem.
```bash
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/library/
export LD_LIBRARY_PATH
```

This means, however, whenever you change the `TPCObjectContainer` or `ParticleContainer` classes you need to make sure you recompile and ensure the new libraries are visible to your script.
This includes the python scripts, as they directly call the c++ modules.

###Calculating the Cross Section

The cross section is calculated at the end of the scripts/selection.cxx code. Using whatever the final results of the selection are, the values are plugged into the standard formula to caculate the cross section. This is mostly striaght-forward, however calculating the flux can be tricky.

To properly calculate the flux we use three separate scripts. The first is found in potFinder, where you can run `potfinder.fcl` over the art ROOT files to produce `pot_extraction.root`, which contains the POT histogram for your files. This will take several minutes given the speed of access to your files. Then we can use  `scripts/potSum.cxx` to extract the total POT from the histogram.
I also have a custom_potSum.cxx which is an interation of finding the POT of MC files, as I now save the POT and subrun information.

Lastly, we need to consulte GENIE and the NuMI flux histograms to find the proper scaling value. This particular scaling factor will change based on your desired signal and beam. If you'd like to use it, I'll be happy to provide you with the scaling factor I use or the module and source files. Those on the gpvms can access `/uboone/app/users/chill2/flux_scripts`. The script is simple enough so I simply execute using root: `root -l flux_calc.cxx+`.


## Open Tickets

This is a place where I will be tracking some issues:

- Some MC Variables I'm not using are sometimes confused when the pfparticle comes from a cosmic.
- Some MC Truth variables are simply not set - `xsecAna::TPCObjectContainer.MCHits()` -- this is unused.
- Reco Momentum for recob::Track objects using StartMomentum() is always 1.
- Cosmic information using MCParticles is not great (vtx where Corsika generated, etc), need to use MCTrack/MCShower for these, but no associations in the MC files.
