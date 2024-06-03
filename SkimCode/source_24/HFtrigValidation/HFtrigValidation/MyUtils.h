#ifndef __MYUTILS_H__
#define __MYUTILS_H__

#include "xAODTracking/TrackParticle.h"

namespace MyUtils{
   enum{
      PP_MIN_BIAS=2,
      HI_LOOSE=4,
      HI_TIGHT=8,
      HI_TIGHT_TIGHTER_D0_Z0=16,
      HI_LOOSE_7SCT_HITS    =32,
      HI_TIGHT_LOOSE_D0_Z0  =64,
   };

   int TrackQuality(const xAOD::TrackParticle* track,float z_vtx);
}
#endif
