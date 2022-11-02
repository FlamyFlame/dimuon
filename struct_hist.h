#ifndef struct_hist_h
#define struct_hist_h

struct Hist2D{
	std::string name;
	bool projx;
	bool projy;
	std::string hxName;
	std::string hyName;
	bool logx;
	bool logy;
	bool logz;
	std::string name_specifier="";
};

struct Hist1D{
	std::string name;
	bool logx;
	bool logy;
	std::string name_specifier="";
};

#endif