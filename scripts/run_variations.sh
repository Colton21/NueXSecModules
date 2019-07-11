# bash script to run all the variations in one go
# To run do: source run_variations.sh
run(){
    COLOR="\033[1;33m"
    COLOR="\033[1;33m"
    DEFAULT="\033[0m"
    echo -e "${COLOR}-> ${1}${DEFAULT}";
    eval ${1};
}

run './main.exe --var_mode files/filter_CV.root "";'
run './main.exe --var_mode files/filter_withDIC.root "";'
run './main.exe --var_mode files/filter_stretchResp.root "";'
run './main.exe --var_mode files/filter_deadSaturatedChannels.root "";'
run './main.exe --var_mode files/filter_altDeadChannels.root "";'
run './main.exe --var_mode files/filter_EnhancedTPCVis.root "";'
run './main.exe --var_mode files/filter_squeezeResp.root "";'
run './main.exe --var_mode files/filter_upPEnoise.root "";'
run './main.exe --var_mode files/filter_noiseAmpDown.root "";'
run './main.exe --var_mode files/filter_downPEnoise.root "";'
run './main.exe --var_mode files/filter_noiseAmpUp.root "";'
run './main.exe --var_mode files/filter_DTdown.root "";'
run './main.exe --var_mode files/filter_DTup.root "";'
run './main.exe --var_mode files/filter_DLdown.root "";'
run './main.exe --var_mode files/filter_DLup.root "";'
run './main.exe --var_mode files/filter_dataSCE.root "";'
run './main.exe --var_mode files/filter_LArG4BugFix.root "";'
run './main.exe --var_mode files/filter_BirksRecomb.root "";'

# run this again for saving the plot
run './main.exe --var_mode files/filter_EnhancedTPCVis.root "same";'
