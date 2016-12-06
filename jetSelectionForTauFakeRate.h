#ifndef JETSEL_H
#define JETSEL_H

/*
 * couple useful functions for jet-tau fakerates
 */

#include <iostream>
#include <map>
#include <string>

#include "TH3F.h"

using namespace std;




// good bins 1, 2
// Float_t bins_pt[11] = { 0, 29, 33, 37, 40, 43, 45, 48, 56, 63, 500 }; // 10 bins, 11 edges
//Float_t bins_pt[11] = { 0, 30, 33, 37, 40, 43, 45, 48, 56, 63, 500 }; // 10 bins, 11 edges
Float_t bins_pt[12] = { 0, 20, 25, 30, 35, 40, 45, 50, 60, 80, 150, 500 }; // 11 bins, 12 edges
int n_bins_pt = 11;
// bins 3 (exactly from AN489)
//Float_t bins_pt[12] = { 0, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100, 150, 500 }; // 13 bins, 14 edges
//int n_bins_pt = 13;

// Float_t bins_eta[6] = { -3, -1.5, -0.45, 0.45, 1.5, 3 }; // 5 bins, 6 edges
//Float_t bins_eta[8] = { -3, -2.5, -1.5, -0.45, 0.45, 1.5, 2.5, 3 }; // 7 bins, 8 edges
//int n_bins_eta = 7;
Float_t bins_eta[8] = { -3, -2.4, -1.5, -0.45, 0.45, 1.5, 2.4, 3 }; // 7 bins, 8 edges
int n_bins_eta = 7;

//Float_t bins_rad[16] = { 0, 0.06, 0.07, 0.08, 0.087, 0.093, 0.1, 0.107, 0.113, 0.12,
	//0.127, 0.133, 0.14, 0.15, 0.16, 2 }; // 15 bins, 16 edges
//int n_bins_rad = 15;
// exactly from AN489:
Float_t bins_rad[10] = { 0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 2 }; // 9 bins, 10 edges
int n_bins_rad = 9;



extern std::map<string, std::map<string, TH3D>> th3d_distr_maps_control;
extern std::map<string, TH3D> th3d_distr_maps_control_headers;



int fill_jet_distr(string control_point_name, Double_t weight, Double_t pt, Double_t eta, Double_t radius)
	{
	// for tau (and other) fake-rates
	// check if the key (mc_decay, control point) has been initialized
	// std::pair <string,string> key (mc_decay, control_point_name);

	// create channel map in th3d_distr_maps_control
	if (th3d_distr_maps_control.find(mc_decay) == th3d_distr_maps_control.end() )
		{
		//
		th3d_distr_maps_control.insert( std::make_pair(mc_decay, std::map<string, TH3D>()));
		}

	if (th3d_distr_maps_control[mc_decay].find(control_point_name) == th3d_distr_maps_control[mc_decay].end() )
		{
		// the control point distr has not been created/initialized
		// THE HISTOGRAM NAMES HAVE TO BE DIFFERENT
		// BECAUSE O M G ROOT USES IT AS A REAL POINTER TO THE OBJECT
		//         O M G
		th3d_distr_maps_control[mc_decay][control_point_name] = TH3D((control_point_name + mc_decay).c_str(), ";;", n_bins_pt, bins_pt, n_bins_eta, bins_eta, n_bins_rad, bins_rad);
		// th3f_distr_control.insert( std::make_pair(key, TH3F(control_point_name.c_str(),      ";;", 10, bins_pt, n_bins_eta, bins_eta, 15, bins_rad)));
		// later on, when writing to the file,
		// I'll have to rename histograms on each write
		// and probably delete them along the way, so that they don't collide...
		// ROOT SUCKS
		//cout << "creating " << control_point_name << endl;
		}


	// fill the distribution:
	th3d_distr_maps_control[mc_decay][control_point_name].Fill(pt, eta, radius, weight);

	if (th3d_distr_maps_control_headers.find(control_point_name) == th3d_distr_maps_control_headers.end() )
		{
		// th3d_distr_maps_control_headers[control_point_name] = TH3D("Header of Pt/E distributions", ";;Pt/E(GeV)", 400, 0., 400.);
		// th3f_distr_control_headers.insert( std::make_pair(string("j_distr"), TH3F("Header of jets distribution",      ";;", 10, bins_pt, 5, bins_eta, 15, bins_rad)));
		th3d_distr_maps_control_headers.insert( std::make_pair(control_point_name, TH3D((string("Header of ") + control_point_name).c_str(), ";;", n_bins_pt, bins_pt, n_bins_eta, bins_eta, n_bins_rad, bins_rad)));
		}

	// return success:
	return 0;
	}



#endif /* JETSEL_H */

