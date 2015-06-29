#!/bin/bash

ls /direct/phenix+spin/phnxsp01/rseidl/Run13pp510Muon/ana103/3*.root \
   > filelist_run13pp510_data_reproduction.txt

ls /direct/phenix+spin2/rseidl/Run12pp510Muon/1655/data/3* \
   > filelist_run12pp510_data.txt

#ls /direct/phenix+spin2/rseidl/Run13pp510Muon/2103/data/3* \
ls /direct/phenix+spin2/rseidl/taxi/Run13pp510Muon/3437/data/3* \
   > filelist_run13pp510_data.txt

#ls /direct/phenix+spin/phnxsp01/rseidl/taxi/Run13pp510OT/2765/data/3* \
ls /direct/phenix+spin2/rseidl/taxi/Run13pp510OT/3438/data/3* \
   > filelist_run13pp510_data_OT.txt

ls /direct/phenix+spin2/rseidl/Wsims/run12sim/muonsims/simdata/pytune100_367593/old2/w_sum/pdst*.root \
   > filelist_run12_sim367593_w.txt

ls /direct/phenix+spin2/rseidl/Wsims/run12sim/muonsims/simdata/pytune100_368630/old2/w_sum/pdst*.root \
   > filelist_run12_sim368630_w.txt

ls /direct/phenix+spin2/rseidl/Wsims/run12sim/muonsims/simdata/pytune100_367466/old2/w_sum/pdst*.root \
   > filelist_run12_sim367466_w.txt

ls /direct/phenix+spin2/rseidl/Wsims/run12sim/muonsims/simdata/pytune100_367593/old2/dy_sum/pdst*.root \
        > filelist_run12_sim367593_dy.txt

ls /direct/phenix+spin2/rseidl/Wsims/run12sim/muonsims/simdata/pytune100_367593/old2/light_sum/pdst*.root \
        > filelist_run12_sim367593_light.txt

ls /direct/phenix+spin2/rseidl/Wsims/run12sim/muonsims/simdata/pytune100_367593/old2/onium_sum/pdst*.root \
        > filelist_run12_sim367593_onium.txt

ls /direct/phenix+spin2/rseidl/Wsims/run12sim/muonsims/simdata/pytune100_367593/old2/onlyz_sum/pdst*.root \
        > filelist_run12_sim367593_onlyz.txt

ls /direct/phenix+spin2/rseidl/Wsims/run12sim/muonsims/simdata/pytune100_367593/old2/openbottom_sum/pdst*.root \
        > filelist_run12_sim367593_openbottom.txt

ls /direct/phenix+spin2/rseidl/Wsims/run12sim/muonsims/simdata/pytune100_367593/old2/opencharm_sum/pdst*.root \
        > filelist_run12_sim367593_opencharm.txt

ls /direct/phenix+spin2/rseidl/Wsims/run12sim/muonsims/simdata/pytune100_367593/old2/whad_sum/pdst*.root \
        > filelist_run12_sim367593_whad.txt

ls /direct/phenix+spin2/rseidl/Wsims/run12sim/muonsims/simdata/pytune100_367593/old2/wjet_sum/pdst*.root \
        > filelist_run12_sim367593_wjet.txt

ls /direct/phenix+spin2/rseidl/Wsims/run12sim/muonsims/simdata/pytune100_367593/old2/wtau_sum/pdst*.root \
        > filelist_run12_sim367593_wtau.txt

ls /direct/phenix+spin2/rseidl/Wsims/run12sim/muonsims/simdata/pytune100_367593/old2/z_sum/pdst*.root \
        > filelist_run12_sim367593_z.txt

ls /direct/phenix+spin2/rseidl/Wsims/run12sim/muonsims/simdata/pytune100_367593/old2/zjet_sum/pdst*.root \
        > filelist_run12_sim367593_zjet.txt



ls /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune100_393888/old2/w_sum/pdst*.root \
   > filelist_run13_py100_sim393888_w.txt

ls /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune100_393888/old2/dy_sum/pdst*.root \
        > filelist_run13_py100_sim393888_dy.txt

find /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune100_393888/old2/light_sum/ -name 'pdst*.root' \
        > filelist_run13_py100_sim393888_light.txt

ls /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune100_393888/old2/onium_sum/pdst*.root \
        > filelist_run13_py100_sim393888_onium.txt

ls /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune100_393888/old2/onlyz_sum/pdst*.root \
        > filelist_run13_py100_sim393888_onlyz.txt

ls /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune100_393888/old2/openbottom_sum/pdst*.root \
        > filelist_run13_py100_sim393888_openbottom.txt

find /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune100_393888/old2/opencharm_sum -name  'pdst*.root' \
        > filelist_run13_py100_sim393888_opencharm.txt

ls /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune100_393888/old2/whad_sum/pdst*.root \
        > filelist_run13_py100_sim393888_whad.txt

ls /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune100_393888/old2/wjet_sum/pdst*.root \
        > filelist_run13_py100_sim393888_wjet.txt

ls /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune100_393888/old2/wtau_sum/pdst*.root \
        > filelist_run13_py100_sim393888_wtau.txt

ls /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune100_393888/old2/z_sum/pdst*.root \
        > filelist_run13_py100_sim393888_z.txt

ls /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune100_393888/old2/zjet_sum/pdst*.root \
        > filelist_run13_py100_sim393888_zjet.txt


ls /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune101_393888/old2/w_sum/pdst*.root \
   > filelist_run13_py101_sim393888_w.txt

ls /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune101_393888/old2/dy_sum/pdst*.root \
        > filelist_run13_py101_sim393888_dy.txt

ls /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune101_393888/old2/light_sum/pdst*.root \
        > filelist_run13_py101_sim393888_light.txt

ls /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune101_393888/old2/onium_sum/pdst*.root \
        > filelist_run13_py101_sim393888_onium.txt

ls /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune101_393888/old2/onlyz_sum/pdst*.root \
        > filelist_run13_py101_sim393888_onlyz.txt

ls /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune101_393888/old2/openbottom_sum/pdst*.root \
        > filelist_run13_py101_sim393888_openbottom.txt

ls /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune101_393888/old2/opencharm_sum/pdst*.root \
        > filelist_run13_py101_sim393888_opencharm.txt

ls /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune101_393888/old2/whad_sum/pdst*.root \
        > filelist_run13_py101_sim393888_whad.txt

ls /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune101_393888/old2/wjet_sum/pdst*.root \
        > filelist_run13_py101_sim393888_wjet.txt

ls /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune101_393888/old2/wtau_sum/pdst*.root \
        > filelist_run13_py101_sim393888_wtau.txt

ls /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune101_393888/old2/z_sum/pdst*.root \
        > filelist_run13_py101_sim393888_z.txt

ls /direct/phenix+spin2/rseidl/Wsims/run13sim/muonsims/simdata/pytune101_393888/old2/zjet_sum/pdst*.root \
        > filelist_run13_py101_sim393888_zjet.txt

