#include "HistFillingBaseClass.h"
#include <iostream>
#include <sstream>

/* ===========================================================
 *                      helpers
 * =========================================================*/
TH1D* HistFillingBaseClass::MakeHist1D(const var1D& v,
                                       const std::string& suffix,
                                       int hist_sign)
{
    if (!v.isValid())
        throw std::runtime_error("Invalid var1D definition for "+v.name);

    std::string sfx = suffix;
    if (hist_sign!=-1) sfx += std::to_string(hist_sign+1);

    std::string hname = "h_"+v.name+sfx;
    std::string htitle = ";"+v.title+";d#sigma/d{"+v.title+"}";

    TH1D* h{nullptr};
    if (v.bins.empty())
        h = new TH1D(hname.c_str(),htitle.c_str(),v.nbins,v.vmin,v.vmax);
    else
        h = new TH1D(hname.c_str(),htitle.c_str(),v.nbins,&v.bins[0]);
    h->Sumw2();
    return h;
}

TH2D* HistFillingBaseClass::MakeHist2D(const var1D& vx,const var1D& vy,
                                       const std::string& suffix,int hist_sign)
{
    if (!vx.isValid()||!vy.isValid())
        throw std::runtime_error("Invalid var1D in MakeHist2D");

    std::string sfx = suffix;
    if (hist_sign!=-1) sfx += std::to_string(hist_sign+1);

    std::string hname = "h_"+vy.name+"_vs_"+vx.name+sfx;
    std::string htitle = ";"+vx.title+";"+vy.title;

    TH2D* h{nullptr};
    if (vx.bins.empty() && vy.bins.empty())
        h = new TH2D(hname.c_str(),htitle.c_str(),
                     vx.nbins,vx.vmin,vx.vmax,
                     vy.nbins,vy.vmin,vy.vmax);
    else if (vx.bins.empty())
        h = new TH2D(hname.c_str(),htitle.c_str(),
                     vx.nbins,vx.vmin,vx.vmax,
                     vy.nbins,&vy.bins[0]);
    else if (vy.bins.empty())
        h = new TH2D(hname.c_str(),htitle.c_str(),
                     vx.nbins,&vx.bins[0],
                     vy.nbins,vy.vmin,vy.vmax);
    else
        h = new TH2D(hname.c_str(),htitle.c_str(),
                     vx.nbins,&vx.bins[0],
                     vy.nbins,&vy.bins[0]);
    h->Sumw2();
    return h;
}

/* A generic THn creator (edges only if all vars have custom bins) */
THn* HistFillingBaseClass::MakeHistnD(const std::vector<var1D>& vars,
                                      const std::string& suffix,int hist_sign)
{
    int ndim = static_cast<int>(vars.size());
    std::vector<int>    nbins(ndim);
    std::vector<double> xmin(ndim),xmax(ndim);
    std::vector<const double*> edges(ndim);

    std::ostringstream oss;
    oss<<"h_nd";
    for (auto& v:vars) oss<<"_"<<v.name;
    std::string sfx = suffix;
    if (hist_sign!=-1) sfx+=std::to_string(hist_sign+1);
    std::string hname = oss.str()+sfx;
    std::string htitle=";";

    for (int i=0;i<ndim;++i){
        auto& v=vars[i];
        if (!v.isValid()) throw std::runtime_error("Invalid var in MakeHistnD");
        nbins[i]=v.nbins;
        if (v.bins.empty()){
            xmin[i]=v.vmin; xmax[i]=v.vmax; edges[i]=nullptr;
        } else {
            xmin[i]=v.bins.front(); xmax[i]=v.bins.back(); edges[i]=v.bins.data();
        }
        htitle+=v.title+";";
    }

    THn* h=nullptr;
    if (std::all_of(edges.begin(),edges.end(),[](auto p){return p==nullptr;})){
        h = new THnD(hname.c_str(),htitle.c_str(),ndim,
                     nbins.data(),xmin.data(),xmax.data());
    } else {
        h = new THnD(hname.c_str(),htitle.c_str(),ndim,
                     nbins.data(),edges.data(),edges.data()+ndim);
    }
    h->Sumw2();
    return h;
}

/* ===========================================================
 *            default implementations (override if needed)
 * =========================================================*/
void HistFillingBaseClass::InitInput()
{
    std::string fullpath = in_out_dir_+infile_name_;
    inFile_.reset( TFile::Open(fullpath.c_str(),"READ") );
    if (!inFile_ || inFile_->IsZombie())
        throw std::runtime_error("Cannot open "+fullpath);

    for (int s=0;s<ParamsSet::nSigns;++s){
        std::string tname = "muon_pair_tree_sign"+std::to_string(s+1);
        inTree_[s]=dynamic_cast<TTree*>( inFile_->Get(tname.c_str()) );
        if (!inTree_[s])
            throw std::runtime_error("Missing tree "+tname);

        reader_[s]=std::make_unique<TTreeReader>(inTree_[s]);

        /* Example branch hookup */
        /* auto* br = reader_[s]->GetTree()->GetBranch("weight"); */
        /* ... push address into histvar2branch_ as needed       */
    }
}

void HistFillingBaseClass::InitOutput()
{
    std::string fullpath = in_out_dir_+outfile_name_;
    outFile_.reset( TFile::Open(fullpath.c_str(),"RECREATE") );
    if (!outFile_ || outFile_->IsZombie())
        throw std::runtime_error("Cannot create "+fullpath);
}

void HistFillingBaseClass::InitHists()
{
    /* ---- build suffix map first (child may have provided partial) ---- */
    for (auto& f : filter_list_unsigned_)
        filter2suffix_[f] = DefaultSuffixMapper(f);
    for (auto& f : filter_list_signed_)
        filter2suffix_[f] = DefaultSuffixMapper(f);

    /* ---------- UNSIGNED filters ---------- */
    for (auto& f : filter_list_unsigned_){
        auto sfx = filter2suffix_[f];
        /* 1D */
        for (auto& varname : filter2var1Ds_[f]){
            auto& v=var1D_list_.at(varname);
            auto* h = MakeHist1D(v,sfx);
            u_filt_h1_[f].push_back({varname,h});
        }
        /* 2D */
        for (auto& vv : filter2var2Ds_[f]){
            auto* h = MakeHist2D(var1D_list_.at(vv.first),
                                 var1D_list_.at(vv.second),sfx);
            u_filt_h2_[f].push_back({vv,h});
        }
        /* nD */
        for (auto& vect : filter2varnDs_[f]){
            std::vector<var1D> vlist;
            for (auto& n: vect) vlist.push_back(var1D_list_.at(n));
            auto* h = MakeHistnD(vlist,sfx);
            u_filt_hn_[f].push_back({vect,h});
        }
    }

    /* ---------- SIGNED filters ---------- */
    for (auto& f : filter_list_signed_){
        auto sfx = filter2suffix_[f];
        /* 1D */
        for (auto& varname : filter2var1Ds_[f]){
            std::pair<TH1D*,TH1D*> hh{nullptr,nullptr};
            for (int s=0;s<ParamsSet::nSigns;++s)
                hh.first  = (s==0)? MakeHist1D(var1D_list_.at(varname),sfx,s)
                                  : hh.first,
                hh.second = (s==1)? MakeHist1D(var1D_list_.at(varname),sfx,s)
                                  : hh.second;
            s_filt_h1_[f].push_back({varname,hh});
        }
        /* 2D */
        for (auto& vv : filter2var2Ds_[f]){
            std::pair<TH2D*,TH2D*> hh{nullptr,nullptr};
            for (int s=0;s<ParamsSet::nSigns;++s)
                hh.first  = (s==0)? MakeHist2D(var1D_list_.at(vv.first),
                                               var1D_list_.at(vv.second),sfx,s)
                                  : hh.first,
                hh.second = (s==1)? MakeHist2D(var1D_list_.at(vv.first),
                                               var1D_list_.at(vv.second),sfx,s)
                                  : hh.second;
            s_filt_h2_[f].push_back({vv,hh});
        }
        /* nD */
        for (auto& vect : filter2varnDs_[f]){
            std::vector<var1D> vlist;
            for (auto& n: vect) vlist.push_back(var1D_list_.at(n));
            std::pair<THn*,THn*> hh{nullptr,nullptr};
            for (int s=0;s<ParamsSet::nSigns;++s)
                hh.first  = (s==0)? MakeHistnD(vlist,sfx,s) : hh.first,
                hh.second = (s==1)? MakeHistnD(vlist,sfx,s) : hh.second;
            s_filt_hn_[f].push_back({vect,hh});
        }
    }
}

void HistFillingBaseClass::ProcessData()
{
    for (int s=0;s<ParamsSet::nSigns;++s){
        auto& rdr = reader_[s];
        while (rdr->Next()){
            FillHistograms(s);
        }
    }
}

bool HistFillingBaseClass::PassSingleMuonGapCut(float,float,int){ return true; }

void HistFillingBaseClass::FillHistograms(int nsign)
{
    FillHistogramsPerFilter("DEFAULT",nsign,false);  // unsigned DEFAULT
}

void HistFillingBaseClass::FillHistogramsPerFilter(const std::string& filter,
                                                   int nsign,bool filter_signed)
{
    if (!filter_signed){   /* ------- unsigned ------- */
        auto it=u_filt_h1_.find(filter);
        if (it==u_filt_h1_.end()) { std::cerr<<"Missing filter "<<filter<<"\n"; return; }
        for (auto& [vname,h] : it->second){
            float* val = histvar2branch_[vname];
            if (!val){ std::cerr<<"Missing branch "<<vname<<"\n"; continue; }
            h->Fill(*val/*,weight*/);   // add weight array if needed
        }
        /* idem for 2D/ND … */
    } else {               /* ------- signed ------- */
        auto it=s_filt_h1_.find(filter);
        if (it==s_filt_h1_.end()){ std::cerr<<"Missing signed filter "<<filter<<"\n"; return; }
        for (auto& [vname,hh] : it->second){
            float* val = histvar2branch_[vname];
            if (!val){ std::cerr<<"Missing branch "<<vname<<"\n"; continue; }
            TH1D* h = (nsign==0)? hh.first : hh.second;
            h->Fill(*val);
        }
        /* idem for 2D/ND … */
    }
}

void HistFillingBaseClass::WriteOutput()
{
    outFile_->cd();
    /* walk maps & Write() every histogram */
    auto writeVec=[&](auto& mp){
        for (auto& kv : mp)
            for (auto& pr : kv.second)
                pr.second->Write();
    };
    writeVec(u_filt_h1_); writeVec(u_filt_h2_); writeVec(u_filt_hn_);
    for (auto& kv : s_filt_h1_)
        for (auto& pr : kv.second){ pr.second.first->Write(); pr.second.second->Write(); }
    /* idem for h2/hN */
    outFile_->Write();
    outFile_->Close();
}

void HistFillingBaseClass::Run()
{
    InitInput();
    FillVar1DList();
    InitFilterMaps();
    InitOutput();
    InitHists();
    ProcessData();
    WriteOutput();
}

void HistFillingBaseClass::PrintInstructions()
{
    std::cout<<"Child class must override:\n"
             <<"  • FillVar1DList() – enumerate all var1D objects\n"
             <<"  • InitFilterMaps() – list filters & variable combos\n"
             <<"  • (optionally) PassSingleMuonGapCut, FillHistograms, etc.\n";
}
