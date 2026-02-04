#pragma once

#include <TChain.h>
#include <TBranch.h>
#include <stdexcept>
#include <string>
#include <sstream>

inline void bind_branch(TChain* ch,
                        const char* name,
                        void* addr,
                        bool use_TBranch = false, // this seems to fail when th stored data type in branch isn't expected when doing make class
                        bool write_debug = false)
{
    if (!ch) {
        throw std::runtime_error("bind_branch: TChain* is null");
    }

    TBranch* b = nullptr;
    int rc = -100;
    if (use_TBranch){
        rc = ch->SetBranchAddress(name, addr, &b);
        if(write_debug) std::cout << "rc: " << rc << ", name: " << name << ", b: " << std::endl;
    } else{
        rc = ch->SetBranchAddress(name, addr);
        if(write_debug) std::cout << "rc: " << rc << ", name: " << name << ", b: " << b << std::endl;
    }

    bool invalid = (rc != 0 && rc != 3);
    if (use_TBranch) invalid |= !b;
    if (invalid) {
        std::ostringstream oss;
        oss << "bind_branch: SetBranchAddress failed for '" << name
            << "' rc=" << rc
            << " branch_ptr=" << b
            << " (chain=" << ch << ", treename=" << ch->GetName() << ")";
        throw std::runtime_error(oss.str());
    }
}

inline void enable_and_bind(TChain* ch, const char* name, void* addr)
{
    if (!ch) throw std::runtime_error("enable_and_bind: TChain* is null");

    ch->SetBranchStatus(name, 1);
    bind_branch(ch, name, addr);
}
