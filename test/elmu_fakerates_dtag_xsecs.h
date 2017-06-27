#include "TString.h"

// nick and colour for dtags
std::pair<TString, Color_t> dtag_nick_colour_elmufakerates(TString dtag)
	{
	if (dtag.Contains("Data")) return std::make_pair("data", kWhite);
	else if(dtag.Contains("DYJets")) return std::make_pair("DY-jets", kGray);
	else if(dtag.Contains("W0Jets") ||dtag.Contains("W4Jets") ||dtag.Contains("W1Jets") ||dtag.Contains("W2Jets") ||dtag.Contains("W3Jets") ||dtag.Contains("WJets") ) return std::make_pair("W-jets", kRed+1);
	else if(dtag.Contains("WW") ||dtag.Contains("WZ") ||dtag.Contains("ZZ")) return std::make_pair("dibosons", kCyan);
	else if(dtag.Contains("Single") || dtag.Contains("schannel") ||dtag.Contains("tchannel")) return std::make_pair("single top", kAzure);
	else if(dtag.Contains("TT"))
		{
		if (dtag.Contains("qq")) return std::make_pair("tt \\rightarrow jj", kGreen+4);
		else if (dtag.Contains("elq") || dtag.Contains("mqu") || dtag.Contains("qtaubar")) return std::make_pair("tt \\rightarrow \\ell j", kGreen+3);
		else if (dtag.Contains("elmu") || dtag.Contains("elmubar") || dtag.Contains("muelbar")) return std::make_pair("tt \\rightarrow e\\mu", kOrange);
		else if (dtag.Contains("eell")) return std::make_pair("tt \\rightarrow ee", kAzure-9);
		else if (dtag.Contains("mmuu")) return std::make_pair("tt \\rightarrow \\mu\\mu", kGreen-9);
		else if (dtag.Contains("aehltu") || dtag.Contains("ahmtuu") || dtag.Contains("aaehlttuu") || dtag.Contains("aahmttuuu")) // hadronic tau ("true" tau in the reconstruction)
			return std::make_pair("tt \\rightarrow \\ell\\tau_{h}}", kOrange+2);
		else if (dtag.Contains("aelmtuu") || dtag.Contains("aaelmttuuu")) // leptonic taus (one or two)
			return std::make_pair("tt \\rightarrow \\tau_{\\ell} (\\ell)", kOrange+3);
		else if (dtag.Contains("atqu")) return std::make_pair("tt \\rightarrow \\tau j", kCyan-5);
		else if (dtag.Contains("aeltqu") || dtag.Contains("amtquu")) return std::make_pair("tt \\rightarrow \\tau_{\\ell} j", kCyan-5);
		else return std::make_pair("tt_{other}", kCyan-5);
		}
	else if(dtag.Contains("QCD")) return std::make_pair("qcd", kViolet);
	else return std::make_pair("other", kBlack);
	}

// nick and colour for dtags
std::pair<TString, Color_t> dtag_nick_colour_elmufakerates_reduced(TString dtag)
	{
	if (dtag.Contains("Data")) return std::make_pair("data", kWhite);
	else if(dtag.Contains("DYJets")) return std::make_pair("DY-jets", kGray);
	else if(dtag.Contains("W0Jets") ||dtag.Contains("W4Jets") ||dtag.Contains("W1Jets") ||dtag.Contains("W2Jets") ||dtag.Contains("W3Jets") ||dtag.Contains("WJets") ) return std::make_pair("W-jets", kRed+1);
	else if(dtag.Contains("WW") ||dtag.Contains("WZ") ||dtag.Contains("ZZ")) return std::make_pair("dibosons", kCyan);
	else if(dtag.Contains("Single") || dtag.Contains("schannel") ||dtag.Contains("tchannel")) return std::make_pair("single top", kAzure);
	else if(dtag.Contains("TT"))
		{
		if (dtag.Contains("elmu") || dtag.Contains("elmubar") || dtag.Contains("muelbar")) return std::make_pair("tt \\rightarrow e\\mu", kOrange);
		else if (dtag.Contains("aelmtuu") || dtag.Contains("aaelmttuuu")) // leptonic taus (one or two)
			return std::make_pair("tt \\rightarrow \\tau_{\\ell} (\\ell)", kOrange+3);
		else return std::make_pair("tt_{other}", kCyan-5);
		}
	else if(dtag.Contains("QCD")) return std::make_pair("qcd", kViolet);
	else return std::make_pair("other", kBlack);
	}

/*
		if (dtag.Contains("qq")) return std::make_pair("tt \\rightarrow jj", kGreen+4);
		else if (dtag.Contains("elq") || dtag.Contains("mqu") || dtag.Contains("qtaubar")) return std::make_pair("tt \\rightarrow \\ell j", kGreen+3);
		else if (dtag.Contains("eell")) return std::make_pair("tt \\rightarrow ee", kAzure-9);
		else if (dtag.Contains("mmuu")) return std::make_pair("tt \\rightarrow \\mu\\mu", kGreen-9);
		else if (dtag.Contains("aehltu") || dtag.Contains("ahmtuu") || dtag.Contains("aaehlttuu") || dtag.Contains("aahmttuuu")) // hadronic tau ("true" tau in the reconstruction)
			return std::make_pair("tt \\rightarrow \\ell\\tau_{h}}", kOrange+2);
		else if (dtag.Contains("atqu")) return std::make_pair("tt \\rightarrow \\tau j", kCyan-5);
		else if (dtag.Contains("aeltqu") || dtag.Contains("amtquu")) return std::make_pair("tt \\rightarrow \\tau_{\\ell} j", kCyan-5);
*/

