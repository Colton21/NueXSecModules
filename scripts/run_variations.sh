# bash script to run all the variations in one go
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

# run this again for saving the plot
run './main.exe --var_mode files/filter_EnhancedTPCVis.root "same";'
