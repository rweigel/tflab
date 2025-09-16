At 

http://mag.gmu.edu/git-data/rweigel/tflab/EarthScope/ORF03/tfs-20070831T014836-20070904T014835/figures/figures.html

I show a comparison of my calculation of the MT impedance with that in the 
EDI file http://ds.iris.edu/spud/emtf/14866915.

Given that Ex and Ey are so similar, I expect Zxx ~= Zyx and Zxy ~= Zyy. I find this to be true for my calculations, but this is not the case for the EDI Z. See Figures 5-8 where Zxx and Zyx have a similar phase but different magnitudes and Zxy have very different phases and amplitudes. I suspect that I am not using the same data.


phi_yx at 7.314290 s ~ +19 degrees

However, based on

<Period value="7.314290e0" units="secs">
<Z type="complex" size="2 2" units="[mV/km]/[nT]">
<value name="Zxx" output="Ex" input="Hx">-4.567480e-1 1.202799e-1</value>
<value name="Zxy" output="Ex" input="Hy">1.774532e0 1.087520e0</value>
<value name="Zyx" output="Ey" input="Hx">-2.744468e0 -9.216799e-1</value>
<value name="Zyy" output="Ey" input="Hy">6.867481e-1 -2.226799e-1</value>
</Z>

in the XML at https://ds.iris.edu/spudservice/data/21636090, I expect

(180/pi)*atan2(-9.216799e-1, -2.744468e0) ~ -161 degrees

In addition, sometimes I see missing phase information because of the axis
limits, e.g., https://ds.iris.edu/spud/emtf/15014349