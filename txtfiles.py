def writetxt(ra, dec, dered_g, dered_r, points, filename ):
    #write a column-delimited ".txt" file with the data points in it
    fhandle = open(filename, 'w')
    for j in range(points):
        a, b, c, d = str(ra[j]), str(dec[j]), str(dered_g[j]), str(dered_r[j])
        fhandle.write( a + '\t' + b + '\t' + c + '\t' + d + '\n' )
    fhandle.close()
    return 1

#svn list file:///data2/svn/newrepos/newby/python/
#svn checkout svn+ssh://fornax.phys.rpi.edu/data2/svn/newrepos/newby/
