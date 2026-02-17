#pragma once
#include <string>

struct PlotCommonConfig
{
	std::map<std::string, std::string> var_title_map = {
        {"pair_pt", "p_T^{pair} [GeV]"},
        {"minv_zoomin", "m_{#mu#mu} [GeV]"},
        {"dr_zoomin", "#Delta R"}
    };

    std::map<std::string, bool> log_map = {
        {"pair_pt", true},
        {"minv_zoomin", false},
        {"dr_zoomin", false},
    };
};