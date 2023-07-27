rms_old = 47
spidx = 1.8
v_old = 0.15
v_new = 0.65

print("for alpha=", -spidx, "the rms at", v_new, "GHz is", ((v_old/v_new)**spidx)*rms_old*1.2, r"$ \mu $ Jy/bea " )