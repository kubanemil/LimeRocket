import matplotlib.pyplot as plt
import numpy as np
from engine import Engine
from random import shuffle

ps = [p for p in range(1, 30)]
mdots = [mdot/10 for mdot in range(1, 20, 2)]
ofs = [0.5, 1, 2]
shuffle(mdots); shuffle(ps)
ps = ps[:10]
mdots = mdots[:10]
engs = []
for of in ofs:
    for p, mdot in zip(ps, mdots):
        engs.append(Engine(ch_pressure=p, mdot=mdot, o_f=of))

for i in range(len(engs)):
    engs[i].plot_engine(ix=i, save=True)

