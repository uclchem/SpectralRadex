from spectralradex import radex

for data_file in ['hcn.dat', 'o-nh3.dat', 'p-h3o+.dat', 'hc3n.dat', 'catom.dat', 'sio.dat', 'ch2_h2_para.dat', 'hnc.dat', 'hcl.dat', 'ch2_h2_ortho.dat', 'co.dat', 'hco+.dat', 'oh2s.dat', 'hd.dat', 'oh.dat', 'oh@hfs.dat', 'oh2cs.dat', 'n+.dat', 'hcl@hfs.dat', 'hcn@hfs.dat', 'oatom.dat', 'SO-pH2.dat', 'o-h3o+.dat', 'ph2cs.dat', 'o-c3h2.dat', 'c+.dat', 'ph2s.dat', 'p-c3h2.dat', 'o-sic2.dat', 'so2@lowT.dat', 'o2.dat', 'p-nh3.dat', 'oh+.dat']:
    print(data_file)
    line_df=radex.get_transition_table(data_file)
    print(line_df.head())