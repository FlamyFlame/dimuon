#include "HFtrigValidation/MyUtils.h"
#include "TMath.h"


namespace MyUtils{

  int TrackQuality(const xAOD::TrackParticle* track,float z_vtx){
    if(!track) return -1;


//-------------------------------------------------------------------------------------------------
      float pt      = track->pt();
      float eta     = track->eta();
      //float phi     = track->phi();

      int   n_Ipix_hits     =track->auxdata<unsigned char>("numberOfInnermostPixelLayerHits"      );
      int   n_Ipix_expected =track->auxdata<unsigned char>("expectInnermostPixelLayerHit"         );
      int   n_NIpix_hits    =track->auxdata<unsigned char>("numberOfNextToInnermostPixelLayerHits");
      int   n_NIpix_expected=track->auxdata<unsigned char>("expectNextToInnermostPixelLayerHit"   );
      int   n_sct_hits      =track->auxdata<unsigned char>("numberOfSCTHits");
      int   n_pix_hits      =track->auxdata<unsigned char>("numberOfPixelHits");
      int   n_sct_holes     =track->auxdata<unsigned char>("numberOfSCTHoles");
      //int   n_pix_holes     =track->auxdata<unsigned char>("numberOfPixelHoles");
      int   n_sct_dead      =track->auxdata<unsigned char>("numberOfSCTDeadSensors");
      int   n_pix_dead      =track->auxdata<unsigned char>("numberOfPixelDeadSensors");

      float chi2=track->chiSquared();
      float ndof=track->numberDoF();
      //float chi2=track->auxdata<float>("chiSquared");
      //float ndof=track->auxdata<float>("numberDoF");

      float d0      = track->d0();
      float z0_wrtPV= track->z0()+track->vz()-z_vtx;
      float theta   = track->theta();

      //ATH_MSG_INFO("pt=%f eta=%f phi=%f BLh=%d BLe=%d pixh=%d pixho=%d scth=%d sctho=%d chi2=%f dof=%d d0=%d
      //z0=%f",pt,eta,phi,n_Ipix_hits,n_Ipix_expected,n_pix_hits,n_pix_holes,n_sct_hits,n_sct_holes,chi2,ndof,); 

      //float e_z0_wrtPV    = track->
      //float e_theta = track->
      //float e_d0    = track->
      //float cov     = track->
      //float e_z0_wrtPV_sin= sqrt( fabs( pow( (double)(e_z0_wrtPV*sin(theta)),2.0 ) +
      // pow( (double) (z0_wrtPV*cos(theta)*e_theta),2.0) + 2*sin(theta)*z0_wrtPV*cos(theta)*cov ));
//-------------------------------------------------------------------------------------------------




//-------------------------------------------------------------------------------------------------
    if(fabs(eta)>2.5) return 0;

    //---------------------------------------------------------------
    bool pass_min_bias=true;
    {
      if(n_Ipix_expected>0){
        if (n_Ipix_hits==0) pass_min_bias=false;
      }
      else{
        if(n_NIpix_expected>0 && n_NIpix_hits==0) pass_min_bias=false;
      }

      int n_sct=n_sct_hits+n_sct_dead;
      if     (pt<=300) {if (n_sct <2)  pass_min_bias=false;}
      else if(pt<=400) {if (n_sct <4)  pass_min_bias=false;}
      else if(pt> 400) {if (n_sct <6)  pass_min_bias=false;}

      int n_pix=n_pix_hits+n_pix_dead;
      if(n_pix<=0) pass_min_bias=false;

      if(fabs(d0)>1.5) pass_min_bias=false;
      if(fabs(z0_wrtPV*sin(theta))>1.5) pass_min_bias=false;

      if(pt>10000 && TMath::Prob(chi2,ndof)<=0.01) pass_min_bias=false;
      //if(n_sct_holes>1 || n_pix_holes>0) continue;
      //if(n_pix_hits<3 || n_sct_hits<8) continue;
    }
    //---------------------------------------------------------------


    //---------------------------------------------------------------
    bool pass_hi_loose=true;
    {
      if(n_Ipix_expected>0){
        if (n_Ipix_hits==0) pass_hi_loose=false;
      }
      else{
        if(n_NIpix_expected>0 && n_NIpix_hits==0) pass_hi_loose=false;
      }

      if(n_pix_hits==0) pass_hi_loose=false;
      if(n_sct_hits< 6) pass_hi_loose=false;
      if(pt>10000 && TMath::Prob(chi2,ndof)<=0.01) pass_hi_loose=false;
      if(fabs(d0) >1.5) pass_hi_loose=false;
      if(fabs(z0_wrtPV*sin(theta))>1.5) pass_hi_loose=false;
    }
    //---------------------------------------------------------------


    //---------------------------------------------------------------
    bool pass_hi_loose_additional_SCT_hit=true;
    if(!pass_hi_loose) pass_hi_loose_additional_SCT_hit=false;
    else{
      if(n_sct_hits<7) pass_hi_loose_additional_SCT_hit=false;
    }
    //---------------------------------------------------------------


    //---------------------------------------------------------------
    bool pass_hi_tight_loose_d0_z0=true;
    if(!pass_hi_loose) pass_hi_tight_loose_d0_z0=false;
    else{
      if(n_pix_hits <2  ) pass_hi_tight_loose_d0_z0=false;
      if(n_sct_hits <8  ) pass_hi_tight_loose_d0_z0=false;
      if(n_sct_holes>1  ) pass_hi_tight_loose_d0_z0=false;
      if(ndof==0) pass_hi_tight_loose_d0_z0=false;
      else if(chi2/ndof>6) pass_hi_tight_loose_d0_z0=false;
    }
    //---------------------------------------------------------------



    //---------------------------------------------------------------
    bool pass_hi_tight=true;
    if(!pass_hi_loose) pass_hi_tight=false;
    else{
      if(n_pix_hits <2  ) pass_hi_tight=false;
      if(n_sct_hits <8  ) pass_hi_tight=false;
      if(n_sct_holes>1  ) pass_hi_tight=false;
      if(fabs(d0)   >1.0) pass_hi_tight=false;
      if(fabs(z0_wrtPV*sin(theta))>1.0) pass_hi_tight=false;
      if(ndof==0) pass_hi_tight=false;
      else if(chi2/ndof>6) pass_hi_tight=false;
    }
    //---------------------------------------------------------------


    //---------------------------------------------------------------
    bool pass_hi_tight_tighter_d0_z0=true;
    if(!pass_hi_tight) pass_hi_tight_tighter_d0_z0=false;
    else{
      if(fabs(d0)>0.5 || fabs(z0_wrtPV*sin(theta))>0.5) pass_hi_tight_tighter_d0_z0=false;
    }
    //---------------------------------------------------------------
//-------------------------------------------------------------------------------------------------

    int quality=0;
    if(pass_min_bias)                    quality+=PP_MIN_BIAS;
    if(pass_hi_loose)                    quality+=HI_LOOSE;
    if(pass_hi_tight)                    quality+=HI_TIGHT;
    if(pass_hi_tight_tighter_d0_z0)      quality+=HI_TIGHT_TIGHTER_D0_Z0;
    if(pass_hi_loose_additional_SCT_hit) quality+=HI_LOOSE_7SCT_HITS;
    if(pass_hi_tight_loose_d0_z0)        quality+=HI_TIGHT_LOOSE_D0_Z0;

    return quality;
  }
}


