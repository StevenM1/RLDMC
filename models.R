## The following models are available:

setupModel <- function(modelName) {
  # Loads model, makes model.dmc object
  # alba-RL -----------------------------------------------------------------
  if(modelName == 'alba-RL-mag') {
    load_model ("LBA", "alba-RL-mag.R")
    model <- model.dmc(p.map=list(A="1",t0="1",st0="1",
                                  sd_v="1",
                                  B0="1",
                                  SR="1", aV="1",
                                  V0="1", wV="1", wS='1'),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(st0=0, sd_v=1,
                                   SR=-10),
                       factors=list(S=c("s1")), 
                       responses=c("r1","r2"),
                       type="norm")
  # alba-RL-mag -----------------------------------------------------------------
  } else if(modelName == 'alba-RL') {
    load_model ("LBA", "alba-RL.R")
    model <- model.dmc(p.map=list(A="1",t0="1",st0="1",
                                  sd_v="1",
                                  B0="1",
                                  SR="1", aV="1",
                                  V0="1", wV="1"),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(st0=0, sd_v=1,
                                   SR=-10),
                       factors=list(S=c("s1")), 
                       responses=c("r1","r2"),
                       type="norm")
  # arw-RL -----------------------------------------------------------------
  } else if(modelName == 'arw-RL') {
    load_model ("RW", "arw-RL.R")
    model <- model.dmc(p.map=list(A="1",     # start punt distribution, not implemented!
                                  t0="1",
                                  st0="1",    # not implemented!
                                  s="1",
                                  B0="1",
                                  SR="1", aV="1",
                                  V0="1", wV="1"),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(st0=0, s=1,
                                   SR=-10, A=0),  # <- keep these fixed!
                       factors=list(S=c("s1")), 
                       responses=c("r1","r2"),
                       type="norm")
  # arw-RL-mag -----------------------------------------------------------------
  } else if(modelName == 'arw-RL-mag') {
    load_model("RW", "arw-RL-mag.R")
    model <- model.dmc(p.map=list(A="1",t0="1",st0="1",s="1",
                                  B0="1", wS='1',
                                  SR="1", aV="1",
                                  V0="1", wV="1"),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(st0=0, s=1,
                                   SR=-10, A=0),
                       factors=list(S=c("s1")), 
                       responses=c("r1","r2"),
                       type="norm")
  # arw-RL-mag-SAT-B2 -----------------------------------------------------------------
  } else if(modelName == 'arw-RL-mag-SAT-B2') {
    load_model ("RW", "arw-RL-mag-mods.R")
    model <- model.dmc(p.map=list(A="1",st0="1",s="1",
                                  t0="1",
                                  B0="1",
                                  SR="1", aV="1",
                                  V0="1", wV="1",
                                  driftMod="1",
                                  V0Mod="1",
                                  B0Mod="cue",
                                  aVMod="1",
                                  t0Mod="1",
                                  wS="1"),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(st0=0, s=1, 
                                   driftMod=0,
                                   B0Mod.ACC=0,
                                   V0Mod=0,
                                   t0Mod=0,
                                   aVMod=0,
                                   SR=-10, A=0),
                       factors=list(S=c("s1"), cue=c('ACC', 'SPD')), 
                       responses=c("r1","r2"),
                       type="norm")
  # arw-RL-mag-SAT-BV02 -----------------------------------------------------------------
  } else if(modelName == 'arw-RL-mag-SAT-BV02') {
    load_model ("RW", "arw-RL-mag-mods.R")
    model <- model.dmc(p.map=list(A="1",st0="1",s="1",
                                  t0="1",
                                  B0="1",
                                  SR="1", aV="1",
                                  V0="1", wV="1",
                                  driftMod="1",
                                  V0Mod="cue",
                                  B0Mod="cue",
                                  aVMod="1",
                                  t0Mod="1",
                                  wS="1"),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(st0=0, s=1, 
                                   driftMod=0,
                                   B0Mod.ACC=0,
                                   V0Mod.ACC=0,
                                   t0Mod=0,
                                   aVMod=0,
                                   SR=-10, A=0),
                       factors=list(S=c("s1"), cue=c('ACC', 'SPD')), 
                       responses=c("r1","r2"),
                       type="norm")
  # arw-RL-mag-SAT-BV0W2 -----------------------------------------------------------------
  } else if(modelName == 'arw-RL-mag-SAT-BV0W2') {
    load_model ("RW", "arw-RL-mag-mods.R")
    model <- model.dmc(p.map=list(A="1",st0="1",s="1",
                                  t0="1",
                                  B0="1",
                                  SR="1", aV="1",
                                  V0="1", wV="1",
                                  driftMod="cue",
                                  V0Mod="cue",
                                  B0Mod="cue",
                                  aVMod="1",
                                  t0Mod="1",
                                  wS="1"),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(st0=0, s=1, 
                                   driftMod.ACC=0,
                                   B0Mod.ACC=0,
                                   V0Mod.ACC=0,
                                   t0Mod=0,
                                   aVMod=0,
                                   SR=-10, A=0),
                       factors=list(S=c("s1"), cue=c('ACC', 'SPD')), 
                       responses=c("r1","r2"),
                       type="norm")
  # arw-RL-mag-SAT-BW2 -----------------------------------------------------------------
  } else if(modelName == 'arw-RL-mag-SAT-BW2') {
    load_model ("RW", "arw-RL-mag-mods.R")
    model <- model.dmc(p.map=list(A="1",st0="1",s="1",
                                  t0="1",
                                  B0="1",
                                  SR="1", aV="1",
                                  V0="1", wV="1",
                                  driftMod="cue",
                                  V0Mod="1",
                                  B0Mod="cue",
                                  aVMod="1",
                                  t0Mod="1",
                                  wS="1"),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(st0=0, s=1, 
                                   driftMod.ACC=0,
                                   B0Mod.ACC=0,
                                   V0Mod=0,
                                   t0Mod=0,
                                   aVMod=0,
                                   SR=-10, A=0),
                       factors=list(S=c("s1"), cue=c('ACC', 'SPD')), 
                       responses=c("r1","r2"),
                       type="norm")
  # arw-RL-mag-SAT-V02 -----------------------------------------------------------------
  } else if(modelName == 'arw-RL-mag-SAT-V02') {
    load_model ("RW", "arw-RL-mag-mods.R")
    model <- model.dmc(p.map=list(A="1",st0="1",s="1",
                                  t0="1",
                                  B0="1",
                                  SR="1", aV="1",
                                  V0="1", wV="1",
                                  driftMod="1",
                                  V0Mod="cue",
                                  B0Mod="1",
                                  aVMod="1",
                                  t0Mod="1",
                                  wS="1"),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(st0=0, s=1, 
                                   driftMod=0,
                                   B0Mod=0,
                                   V0Mod.ACC=0,
                                   t0Mod=0,
                                   aVMod=0,
                                   SR=-10, A=0),
                       factors=list(S=c("s1"), cue=c('ACC', 'SPD')), 
                       responses=c("r1","r2"),
                       type="norm")
  # arw-RL-mag-SAT-V0W2 -----------------------------------------------------------------
  } else if(modelName == 'arw-RL-mag-SAT-V0W2') {
    load_model ("RW", "arw-RL-mag-mods.R")
    model <- model.dmc(p.map=list(A="1",st0="1",s="1",
                                  t0="1",
                                  B0="1",
                                  SR="1", aV="1",
                                  V0="1", wV="1",
                                  driftMod="cue",
                                  V0Mod="cue",
                                  B0Mod="1",
                                  aVMod="1",
                                  t0Mod="1",
                                  wS="1"),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(st0=0, s=1, 
                                   driftMod.ACC=0,
                                   B0Mod=0,
                                   V0Mod.ACC=0,
                                   t0Mod=0,
                                   aVMod=0,
                                   SR=-10, A=0),
                       factors=list(S=c("s1"), cue=c('ACC', 'SPD')), 
                       responses=c("r1","r2"),
                       type="norm")
  # arw-RL-mag-SAT-W2 -----------------------------------------------------------------
  } else if(modelName == 'arw-RL-mag-SAT-W2') {
    load_model ("RW", "arw-RL-mag-mods.R") 
    model <- model.dmc(p.map=list(A="1",st0="1",s="1",
                                  t0="1",
                                  B0="1",
                                  SR="1", aV="1",
                                  V0="1", wV="1",
                                  driftMod="cue",
                                  V0Mod="1",
                                  B0Mod="1",
                                  aVMod="1",
                                  t0Mod="1",
                                  wS="1"),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(st0=0, s=1, 
                                   driftMod.ACC=0,
                                   B0Mod=0,
                                   V0Mod=0,
                                   t0Mod=0,
                                   aVMod=0,
                                   SR=-10, A=0),
                       factors=list(S=c("s1"), cue=c('ACC', 'SPD')), 
                       responses=c("r1","r2"),
                       type="norm")
  }
  # arw-RL-risk-mag-Niek ----------------------------------------------------
  else if(modelName == "arw-RL-risk-mag-Niek"){
    # get model specification
    load_model("RW", paste0(modelName, ".R"))
    model <- model.dmc(p.map=list(A="1",t0="1",st0="1",s="1",
                                  B0="1",
                                  SR="1", aV="1",
                                  RR="1", aR="1",
                                  V0="1", wV="1",
                                  wS="1"),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(st0=0, s=1, RR=-10,
                                   SR=-10, A=0),
                       factors=list(S=c("s1")), 
                       responses=c("r1","r2"),
                       type="norm")
  
  # ddm-RL -----------------------------------------------------------------
  } else if(modelName == 'ddm-RL') {
    load_model ("ddm", "ddm-RL.R")
    model <- model.dmc(p.map=list(aV="1", SR='1', m='1',
                                  t0='1', st0='1', d='1',
                                  z='1', sz='1', a='1',
                                  sv='1'),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(SR=-10, z=0.5, d=0, sv=0, sz=0, st0=0),
                       factors=list(S=c("s1")), 
                       responses=c("r1","r2"),
                       type="norm")
  # ddm-RL-st0 & ddm-RL-st02 -----------------------------------------------------------------
  } else if(modelName == 'ddm-RL-st0' | modelName == 'ddm-RL-st02') {
    load_model ("ddm", "ddm-RL.R")
    model <- model.dmc(p.map=list(aV="1", SR='1', m='1',
                                  t0='1', st0='1', d='1',
                                  z='1', sz='1', a='1',
                                  sv='1'),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(SR=-10, z=0.5, d=0),
                       factors=list(S=c("s1")), 
                       responses=c("r1","r2"),
                       type="norm")
  # ddm-RL-svsz -----------------------------------------------------------------
  } else if(modelName == 'ddm-RL-svsz') {
    load_model ("ddm", "ddm-RL.R")
    model <- model.dmc(p.map=list(aV="1", SR='1', m='1',
                                  t0='1', st0='1', d='1',
                                  z='1', sz='1', a='1',
                                  sv='1'),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(SR=-10, z=0.5, d=0, st0=0),
                       factors=list(S=c("s1")), 
                       responses=c("r1","r2"),
                       type="norm")
  # ddm-RL-nonlinear-svszst0 -----------------------------------------------------------------
  } else if(modelName == 'ddm-RL-nonlinear-svszst0') {
    load_model ("ddm", "ddm-RL-nonlinear.R")
    model <- model.dmc(p.map=list(aV="1", SR='1', m='1',
                                  t0='1', st0='1', d='1',
                                  z='1', sz='1', a='1',
                                  vmax='1',
                                  sv='1'),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(SR=-10, z=0.5, d=0),
                       factors=list(S=c("s1")), 
                       responses=c("r1","r2"),
                       type="norm")
  # ddm-RL-nonlinear -----------------------------------------------------------------
  } else if(modelName == 'ddm-RL-nonlinear') {
    load_model ("ddm", "ddm-RL-nonlinear.R")
    model <- model.dmc(p.map=list(aV="1", SR='1', m='1',
                                  t0='1', st0='1', d='1',
                                  z='1', sz='1', a='1',
                                  vmax='1',
                                  sv='1'),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(SR=-10, z=0.5, d=0, sv=0, sz=0, st0=0),
                       factors=list(S=c("s1")), 
                       responses=c("r1","r2"),
                       type="norm")
  # ddm-RL-SAT-a -----------------------------------------------------------------
  } else if(modelName == 'ddm-RL-SAT-a') {
    load_model ("ddm", "ddm-RL-mod.R")
    model <- model.dmc(p.map=list(aV="1", SR='1', m='1',
                                  t0='1', st0='1', d='1',
                                  z='1', sz='1', a='1',
                                  aMod='cue', driftMod='1',
                                  sv='1'),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(SR=-10, z=0.5, d=0, st0=0, sv=0, sz=0, 
                                   aMod.ACC=0, driftMod=0),
                       factors=list(S=c("s1"), cue=c('ACC', 'SPD')), 
                       responses=c("r1","r2"),
                       type="norm")
  # ddm-RL-SAT-a-st0 -----------------------------------------------------------------
  } else if(modelName == 'ddm-RL-SAT-a-st0') {
    load_model ("ddm", "ddm-RL-mod.R")
    model <- model.dmc(p.map=list(aV="1", SR='1', m='1',
                                  t0='1', st0='1', d='1',
                                  z='1', sz='1', a='1',
                                  aMod='cue', driftMod='1',
                                  sv='1'),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(SR=-10, z=0.5, d=0, 
                                   aMod.ACC=0, driftMod=0),
                       factors=list(S=c("s1"), cue=c('ACC', 'SPD')), 
                       responses=c("r1","r2"),
                       type="norm")
  # ddm-RL-SAT-m -----------------------------------------------------------------
  } else if(modelName == 'ddm-RL-SAT-m') {
    load_model ("ddm", "ddm-RL-mod.R")
    model <- model.dmc(p.map=list(aV="1", SR='1', m='1',
                                  t0='1', st0='1', d='1',
                                  z='1', sz='1', a='1',
                                  aMod='1', driftMod='cue',
                                  sv='1'),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(SR=-10, z=0.5, d=0, st0=0, sv=0, sz=0, 
                                   aMod=0, driftMod.ACC=0),
                       factors=list(S=c("s1"), cue=c('ACC', 'SPD')), 
                       responses=c("r1","r2"),
                       type="norm")
  # ddm-RL-SAT-m-st0 -----------------------------------------------------------------
  } else if(modelName == 'ddm-RL-SAT-m-st0') {
    load_model ("ddm", "ddm-RL-mod.R")
    model <- model.dmc(p.map=list(aV="1", SR='1', m='1',
                                  t0='1', st0='1', d='1',
                                  z='1', sz='1', a='1',
                                  aMod='1', driftMod='cue',
                                  sv='1'),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(SR=-10, z=0.5, d=0,
                                   aMod=0, driftMod.ACC=0),
                       factors=list(S=c("s1"), cue=c('ACC', 'SPD')), 
                       responses=c("r1","r2"),
                       type="norm")
  # ddm-RL-SAT-am -----------------------------------------------------------------
  } else if(modelName == 'ddm-RL-SAT-am') {
    load_model ("ddm", "ddm-RL-mod.R")
    model <- model.dmc(p.map=list(aV="1", SR='1', m='1',
                                  t0='1', st0='1', d='1',
                                  z='1', sz='1', a='1',
                                  aMod='cue', driftMod='cue',
                                  sv='1'),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(SR=-10, z=0.5, d=0, st0=0, sv=0, sz=0, 
                                   aMod.ACC=0, driftMod.ACC=0),
                       factors=list(S=c("s1"), cue=c('ACC', 'SPD')), 
                       responses=c("r1","r2"),
                       type="norm")
  # ddm-RL-SAT-am-st0 -----------------------------------------------------------------
  } else if(modelName == 'ddm-RL-SAT-am-st0') {
    load_model ("ddm", "ddm-RL-mod.R")
    model <- model.dmc(p.map=list(aV="1", SR='1', m='1',
                                  t0='1', st0='1', d='1',
                                  z='1', sz='1', a='1',
                                  aMod='cue', driftMod='cue',
                                  sv='1'),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(SR=-10, z=0.5, d=0, #st0=0, sv=0, sz=0, 
                                   aMod.ACC=0, driftMod.ACC=0),
                       factors=list(S=c("s1"), cue=c('ACC', 'SPD')), 
                       responses=c("r1","r2"),
                       type="norm")
  # rw-RL -----------------------------------------------------------------
  } else if(modelName == 'rw-RL') {
    load_model ("RW", "rw-RL.R")
    model <- model.dmc(p.map=list(A="1",t0="1",st0="1",s="1",
                                  B0="1",
                                  SR="1", aV="1",
                                  V0="1", wV="1"),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(st0=0, s=1, A=0, SR=-10),
                       factors=list(S=c("s1")), 
                       responses=c("r1","r2"),
                       type="norm")
  # softmax-RL2 -----------------------------------------------------------------
  } else if(modelName == 'softmax-RL2') {
    load_model ("Softmax", "softmax_RL2.R")
    model <- model.dmc(p.map=list(Beta="1", aV="1", SR='1'),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(SR=0),
                       factors=list(S=c("s1")), 
                       responses=c("r1","r2"),
                       type="norm")
  # softmax-RL2-SAT-beta -----------------------------------------------------------------
  } else if(modelName == 'softmax-RL2-SAT-beta') {
    load_model ("Softmax", "softmax_RL2.R")
    model <- model.dmc(p.map=list(Beta="cue", aV="1", SR='1'),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(SR=0),
                       factors=list(S=c("s1"), cue=c('ACC', 'SPD')), 
                       responses=c("r1","r2"),
                       type="norm")
    # softmax-RL2-SAT-none -----------------------------------------------------------------
  } else if(modelName == 'softmax-RL2-SAT-none') {
    load_model ("Softmax", "softmax_RL2.R")
    model <- model.dmc(p.map=list(Beta="1", aV="1", SR='1'),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(SR=0),
                       factors=list(S=c("s1")), 
                       responses=c("r1","r2"),
                       type="norm")
  }
    else if(modelName == 'arw-RL-timing') {
    load_model("RW-Timing", "arw-RL-timing.R")
    model <- model.dmc(p.map=list(A="1",t0="1",st0="1",s="1",
                                  B0="1", B_T = "1", wS='1', s_T = "1", v_T = "1", t0T = "1",
                                  SR="1", aV="1",
                                  V0="1", wV="1"),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(st0=0, s=1, s_T = 1, SR=-10,
                                    A=0, t0T = 0.05),
                       factors=list(S=c("s1")), 
                       responses=c("r1","r2"),
                       type="norm")
    # arw-RL-mag-SAT-B2 -----------------------------------------------------------------
    }
    else if(modelName == 'arw-RL-timing_B_Tmod') {
    load_model("RW-Timing", "arw-RL-timing_B_Tmod.R")
    model <- model.dmc(p.map=list(A="1",t0="1",st0="1",s="1",
                                  B0="1", B_T = "1", wS='1', s_T = "1", v_T = "1", t0T = "1",
                                  SR="1", aV="1", B_T_mod = "cue",
                                  V0="1", wV="1"),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(st0=0, s=1, s_T = 1, SR=-10,
                                   A=0, t0T = 0.05, B_T_mod.ACC=0),
                       factors=list(S=c("s1"), cue=c('ACC', 'SPD')), 
                       responses=c("r1","r2"),
                       type="norm")
    }
    else if(modelName == 'arw-RL-timing_B_Tmod_NoV0') {
      load_model("RW-Timing", "arw-RL-timing_B_Tmod.R")
      model <- model.dmc(p.map=list(A="1",t0="1",st0="1",s="1", V0 = "1",
                                    B0="1", B_T = "1", wS='1', s_T = "1", v_T = "1", t0T = "1",
                                    SR="1", aV="1", B_T_mod = "cue", wV="1"),
                         match.map=list(M=list(s1=1, s1=2)),
                         constants=c(st0=0, s=1, s_T = 1, SR=-10,
                                     A=0, t0T = 0.05, V0 = 0, B_T_mod.ACC=0),
                         factors=list(S=c("s1"), cue=c('ACC', 'SPD')), 
                         responses=c("r1","r2"),
                         type="norm")
    # arw-RL-mag-SAT-B2 -----------------------------------------------------------------
    
  }
  else if(modelName == 'arw-RL-timing_B_Tmod_NoV0_EduGuess') {
    load_model("RW-Timing-EduGuess", "arw-RL-timing_B_Tmod.R")
    model <- model.dmc(p.map=list(A="1",t0="1",st0="1",s="1", V0 = "1",
                                  B0="1", B_T = "1", wS='1', s_T = "1", v_T = "1", t0T = "1",
                                  SR="1", aV="1", B_T_mod = "cue", wV="1"),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(st0=0, s=1, s_T = 1, SR=-10,
                                   A=0, t0T = 0.05, V0 = 0, B_T_mod.ACC=0),
                       factors=list(S=c("s1"), cue=c('ACC', 'SPD')), 
                       responses=c("r1","r2"),
                       type="norm")
    # arw-RL-mag-SAT-B2 -----------------------------------------------------------------
  }
  else if(modelName == 'arw-RL-timing_v_Tmod') {
    load_model("RW-Timing", "arw-RL-timing_v_Tmod.R")
    model <- model.dmc(p.map=list(A="1",t0="1",st0="1",s="1",
                                  B0="1", B_T = "1", wS='1', s_T = "1", v_T = "1", t0T = "1",
                                  SR="1", aV="1", v_T_mod = "cue",
                                  V0="1", wV="1"),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(st0=0, s=1, s_T = 1, SR=-10,
                                   A=0, t0T = 0.05, v_T_mod.ACC=0),
                       factors=list(S=c("s1"), cue=c('ACC', 'SPD')), 
                       responses=c("r1","r2"),
                       type="norm")
  }
  else if(modelName == 'arw-RL-timing_v_Tmod_EduGuess') {
    load_model("RW-Timing-EduGuess", "arw-RL-timing_v_Tmod.R")
    model <- model.dmc(p.map=list(A="1",t0="1",st0="1",s="1",
                                  B0="1", B_T = "1", wS='1', s_T = "1", v_T = "1", t0T = "1",
                                  SR="1", aV="1", v_T_mod = "cue",
                                  V0="1", wV="1"),
                       match.map=list(M=list(s1=1, s1=2)),
                       constants=c(st0=0, s=1, s_T = 1, SR=-10,
                                   A=0, t0T = 0.05, v_T_mod.ACC=0),
                       factors=list(S=c("s1"), cue=c('ACC', 'SPD')), 
                       responses=c("r1","r2"),
                       type="norm")
  }
  return(model)
}
