import gctools_enc as gce
import files as f

skip = [1,1,1,1,1,1,0,1,1,1,1]
datafile = ['wide_4147.csv',
            'wide_5024.csv',
            'wide_5053.csv',
            'wide_5272.csv',
            'wide_5466.csv',
            'wide_5904.csv',
            'NGC_6205_all.csv', #'wide_6205.csv',
            'wide_6341.csv',
            'wide_7078.csv',
            'wide_7089.csv',
            'wide_Pal5.csv' 
           ]
center_ra = [182.525, 198.228, 199.109, 205.545, 211.36, 229.641, 250.423, 259.1680,
             322.493, 323.362, 229.013]
center_dec = [18.530, 18.164, 17.697, 28.376, 28.53, 2.083, 36.460, 43.1033,
              12.167, -0.826, -0.123]
inner_cut = [0.0170, 0.0800, 0.0500, 0.0, 0.057, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
outer_cut = [0.0700, 0.2500, 0.1400, 0.3000, 0.180, 0.2100, 0.2400, 0.1350,
             0.1900, 0.1400, 0.0870]
background_cut = [0.0800, 0.3000, 0.1700, 0.4000, 0.2500, 0.2800, 0.3300, 0.2000,
                  0.2500, 0.1800, 0.1200]
outer_back_cut = 10.0

for i in range(len(datafile)):
    if i == 7: continue
    if skip[i] == 1: continue
    data = f.read_csv(datafile[i])
    clus_data = gce.make_cut(data, center_ra[i], center_dec[i], inner_cut[i], outer_cut[i])
    back_data = gce.make_cut(data, center_ra[i], center_dec[i], background_cut[i], outer_back_cut)
    clus_out = datafile[i][:-4]+'_cluster.csv'
    back_out = datafile[i][:-4]+'_background.csv'
    if (f.write_csv(clus_data, clus_out)) == 1:
        print '#-cluster file', datafile[i], 'successfully cut and saved as', clus_out
    else:
        print '!!!AN ERROR OCCURED - FILE NOT CUT CORRECTLY!!!'
    if (f.write_csv(back_data, back_out)) == 1:
        print '#-data file', datafile[i], 'successfully cut and saved as', back_out
    else:
        print '!!!AN ERROR OCCURED - FILE NOT CUT CORRECTLY!!!'
print '#---All Done'